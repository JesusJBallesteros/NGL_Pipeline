from ConvertIntanToNWB import *
from datetime import datetime
""" Callback function for when the user begins conversion.

Parameters
----------
None

Returns
-------
None
"""
        
# Get NWB settings and pass these to NWB conversion
intan_filename="info.rhd"
nwb_filename="info.nwb"
session_description=None
blocks_per_chunk=2000
use_compression=False
compression_level=4
lowpass_description=None
highpass_description=None
merge_files=None
subject=None
manual_start_time=datetime.now()

convert_to_nwb(intan_filename=intan_filename,
               nwb_filename=nwb_filename,
               session_description=session_description,
               blocks_per_chunk=blocks_per_chunk,
               use_compression=use_compression,
               compression_level=compression_level,
               lowpass_description=lowpass_description,
               highpass_description=highpass_description,
               merge_files=merge_files,
               subject=subject,
               manual_start_time=manual_start_time)