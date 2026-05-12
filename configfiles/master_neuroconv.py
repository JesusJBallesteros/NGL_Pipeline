import sys
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
from neuroconv.datainterfaces import IntanRecordingInterface
from neuroconv.utils import load_dict_from_file, dict_deep_update

# Arguments passed from intan2NWB_neuroconv.m:
#   sys.argv[1]  absolute path to the INTAN header file (info.rhd or .rhs)
#   sys.argv[2]  absolute path for the output .nwb file
#   sys.argv[3]  absolute path to the project nwb_metadata.yaml in analysisCode\

file_path     = sys.argv[1]
nwbfile_path  = sys.argv[2]
metadata_yaml = sys.argv[3]

# Create the recording interface for INTAN amplifier channels.
interface = IntanRecordingInterface(file_path=file_path, verbose=True)

# Start with whatever metadata neuroconv can auto-extract from the source file.
metadata = interface.get_metadata()

# Load project-level metadata from the YAML and merge it in.
# dict_deep_update: YAML values override auto-extracted ones where keys overlap.
# Fields absent from the YAML (e.g. session_start_time) are left as-is.
if Path(metadata_yaml).is_file():
    metadata_from_yaml = load_dict_from_file(file_path=metadata_yaml)
    metadata = dict_deep_update(metadata, metadata_from_yaml)
    print("[master_neuroconv] Metadata YAML loaded: {}".format(metadata_yaml))
else:
    print("[master_neuroconv] WARNING: metadata YAML not found at: {}".format(metadata_yaml))
    print("[master_neuroconv]   Proceeding with auto-extracted metadata only.")

# -- Metadata diagnostic ------------------------------------------------------
# Remove this block once the full conversion is verified.
print("\n[master_neuroconv] -- Merged metadata (NWBFile section) --")
print("  file_path : {}".format(file_path))
for key, val in metadata.get("NWBFile", {}).items():
    print("  NWBFile.{} : {!r}".format(key, val))
rec = interface.recording_extractor
print("\n  recording : {} channels @ {} Hz, {} samples ({:.1f} s)".format(
    rec.get_num_channels(),
    rec.get_sampling_frequency(),
    rec.get_num_samples(),
    rec.get_num_samples() / rec.get_sampling_frequency()
))
print("[master_neuroconv] ---------------------------------------------------\n")
# -----------------------------------------------------------------------------

# session_start_time is mandatory for NWB.
# The Intan .rhd format does NOT store an absolute recording timestamp, so
# get_metadata() will never populate it. We parse the session date from the
# raw-data folder name (DDMMYYYY -- IKN pipeline convention).
extracted_time = metadata["NWBFile"].get("session_start_time")

if extracted_time is not None:
    # Should not normally happen with Intan, but handle it cleanly if it does.
    if extracted_time.tzinfo is None:
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
        print("[master_neuroconv] WARNING: session_start_time absent from .rhd header "
              "(Intan does not store it). Parsed from folder '{}': {}".format(
              session_folder, fallback_time))
    except (ValueError, IndexError) as e:
        raise ValueError(
            "session_start_time is absent from the .rhd header and the session "
            "folder name '{}' could not be parsed as DDMMYYYY. "
            "Error: {}".format(session_folder, e)
        )

# Run the conversion. overwrite=True handles partial files from interrupted runs.
# Future stages (spike sorting, behavior) will append to this file with
# overwrite=False using their respective neuroconv interfaces.
interface.run_conversion(nwbfile_path=nwbfile_path, metadata=metadata, overwrite=True)
