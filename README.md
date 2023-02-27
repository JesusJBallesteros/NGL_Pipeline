# First edition of 'Ephys-data-pipeline'
Main Scripts, functions and tools to work with electrophysiological data at NGL.

# A first Script 'NGL01_Main' transforms data from INTAN and Deuteron formats. 
The high-pass data will be converted into a .bin file to feed into Kilosort software.
The low-pass data will be converted to a Matlab (.mat) single file with a FieldTrip-formatted structure.
Additionally, data can be converted into Neurodata Without Borders (.nwb) single file, a class of HDF5 file format.

IMPORTANT: The creation of a NWB file is made by 'intan2NWB_wrapper.m', which performs the main INTAN-NWB transformation. This wrapper NEEDS a working python installation. This is because INTAN's tool (INTANToNWB) requires so. See:
https://github.com/Intan-Technologies/IntanToNWB

Python will never show up, neither the user have to interact with it at all. Everything happens from MATLAB.
However, requires Python installed in the machine. This is explained inside the Script.

To date, MATLAB 2021b accepts interaction with Python up to v3.9. It needs the 64 bits version:
(https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)

To check that MATLAB can interact with the Python Modules, look if variable 'pe' is correctly populated when debugging.


Coming soon...
# A second Script 'NGL02_NWB_managing' with commands to manage .nwb (HDF5) files.
This includes how to read electrophysiological series, event codes, etc. from .NWB files.
Also, commands to add data and metadata to the .NWB file.

# A third...
