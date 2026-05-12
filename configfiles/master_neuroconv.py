import sys
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
from neuroconv.datainterfaces import IntanRecordingInterface

# Arguments passed from intan2NWB_neuroconv.m:
#   sys.argv[1]  absolute path to the INTAN header file (info.rhd or .rhs)
#   sys.argv[2]  absolute path for the output .nwb file
#
# Note: no 'NCfolder' argument needed -- neuroconv is a pip-installed package
# and is resolved through the Python environment, not via os.chdir().

file_path    = sys.argv[1]   # e.g. D:\projectname\ABC\12122025\raw\info.rhd
nwbfile_path = sys.argv[2]   # e.g. D:\projectname\ABC\12122025\preprocessing\12122025.nwb

# Create the recording interface for INTAN amplifier channels.
interface = IntanRecordingInterface(file_path=file_path, verbose=True)

# Extract what metadata we can from the source file.
metadata = interface.get_metadata()

# -- Metadata diagnostic ------------------------------------------------------
# Remove this block once the timestamp is confirmed to be populating correctly.
print("\n[master_neuroconv] -- Metadata extracted from .rhd --")
print("  file_path : {}".format(file_path))
nwb_meta = metadata.get("NWBFile", {})
for key, val in nwb_meta.items():
    print("  NWBFile.{} : {!r}".format(key, val))
rec = interface.recording_extractor
print("\n  recording : {} channels @ {} Hz, {} samples ({:.1f} s)".format(
    rec.get_num_channels(),
    rec.get_sampling_frequency(),
    rec.get_num_samples(),
    rec.get_num_samples() / rec.get_sampling_frequency()
))
print("[master_neuroconv] --------------------------------------------------\n")
# -----------------------------------------------------------------------------

# session_start_time is mandatory for NWB.
# The Intan .rhd format does NOT store an absolute recording timestamp in the
# file header, so get_metadata() will never populate this field from the file.
# Primary: use the extracted value if somehow present (future-proof).
# Fallback: parse the session date from the raw-data folder name (DDMMYYYY),
#           which is the IKN pipeline convention (e.g. '12122025' -> 2025-12-12).
extracted_time = metadata["NWBFile"].get("session_start_time")

if extracted_time is not None:
    metadata["NWBFile"]["session_start_time"] = extracted_time.replace(
        tzinfo=ZoneInfo("Europe/Berlin")
    )
    print("[master_neuroconv] session_start_time from header: {}".format(
        metadata["NWBFile"]["session_start_time"]
    ))
else:
    # Parse DDMMYYYY from the immediate parent folder of info.rhd.
    session_folder = Path(file_path).parent.name
    try:
        day   = int(session_folder[0:2])
        month = int(session_folder[2:4])
        year  = int(session_folder[4:8])
        # Time of day is unknown; midnight is the best approximation.
        fallback_time = datetime(year, month, day, 0, 0, 0,
                                 tzinfo=ZoneInfo("Europe/Berlin"))
        metadata["NWBFile"]["session_start_time"] = fallback_time
        print("[master_neuroconv] WARNING: session_start_time not in .rhd header "
              "(Intan does not store it). Parsed from folder name '{}': {}".format(
              session_folder, fallback_time))
    except (ValueError, IndexError) as e:
        raise ValueError(
            "session_start_time is absent from the .rhd header and the session "
            "folder name '{}' could not be parsed as DDMMYYYY. "
            "Error: {}".format(session_folder, e)
        )

# Run the conversion. overwrite=True prevents crashes if a partial .nwb
# file was left behind by a previous interrupted run.
interface.run_conversion(nwbfile_path=nwbfile_path, metadata=metadata, overwrite=True)
