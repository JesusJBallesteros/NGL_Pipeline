# First edition of 'Ephys-data-pipeline'
Main Scripts, functions and tools to work with electrophysiological data at NGL.

# A first Script 'NGL01_Main' transforms data from (most) INTAN formats. 
Data can be converted into:
- Neurodata Without Borders (.nwb) single file. A HDF5 file format. 
- Matlab (.mat) single file. A 'struct' with FieldTrip format.

IMPORTANT: The creation of a NWB file is made by 'intan2NWB_wrapper.m', which performs the main INTAN-NWB transformation. This wrapper NEEDS a working python installation. This is because INTAN's tool (INTANToNWB) requires so. See:
https://github.com/Intan-Technologies/IntanToNWB

With this wrapper, Python will never show up, neither the user will have to interact with it AT ALL. Everything happens from MATLAB.
Requires Python installed in the machine. As explained in the Script:

To date, MATLAB 2021b accepts interaction with Python up to v3.9. Needs the 64 bits version:
(https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)

To check that MATLAB can interact with the Python Modules, look if variable 'pe' is correctly populated when debugging.

# A second Script 'NGL02_NWB_managing' contains information and instructions to manage .nwb (HDF5) files.


# A third...
