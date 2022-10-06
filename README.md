# First edition of 'Ephys-data-pipeline'
Main Scripts, functions and tools to work with electrophysiological data at NGL.

# A first Script 'NGL01_Main' transforms data from (most) INTAN formats. Data can be converted into:
- Neurodata Without Borders (.nwb) single file. A HDF5 file format. 
- Matlab (.mat) single file. A 'struct' with FieldTrip format.

IMPORTANT: To use the function 'intan2NWB_wrapper', which performs the main INTAN-NWB transformation, a working python installation is required.
Python will never show up, neither the user will have to interact because everything happens from MATLAB.
Requires Python installed in the machine. As explained in the Script:

To date, MATLAB 2021b accepts up to Python 3.9. Install the 64 bits version:
(https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
To check access to Python Modules from MATLAB, look that 'pe' is correctly populated when running the script.

# A second Script 'NGL02_NWB_managing' contains information and instructions to manage .nwb (HDF5) files.


# A third...
