import sys
import os

# Take matlab arguments (NCdir,)
neuroconvdir = sys.argv[1]

# Go to NCdir
os.chdir(neuroconvdir)

from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
from neuroconv.datainterfaces import IntanRecordingInterface

file_path = sys.argv[2] # This can be .rhd or .rhs
interface = IntanRecordingInterface(file_path=file_path, verbose=True)

# Extract what metadata we can from the source files
metadata = interface.get_metadata()

# session_start_time is required for conversion. If it cannot be inferred
# automatically from the source files you must supply one.
session_start_time = datetime(2024, 6, 6, 10, 30, 0, tzinfo=ZoneInfo("Europe/Berlin"))
metadata["NWBFile"].update(session_start_time=session_start_time)

nwbfile_path = sys.argv[3]  # This should be something like: "./saved_file.nwb"

interface.run_conversion(nwbfile_path=nwbfile_path, metadata=metadata)