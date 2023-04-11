# First edition of 'Ephys-data-pipeline'
Main Scripts, functions and tools to work with electrophysiological data at NGL.

# A Script 'NGL01_Main' that transforms data from INTAN and Deuteron formats into .bin files (for kilosort) and Fieldtrip .mat files.
Additionally, data can be converted into a Neurodata Without Borders (.nwb) single file, a class of HDF5 file format.

Your data SHOULD be stored as the IKN standard Harddisk data structure. See:
  gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure

  *ABOUT .NWB: The creation of a NWB file is made by 'intan2NWB_wrapper.m', which performs the main INTAN-NWB transformation. This wrapper NEEDS a working python installation. This is because INTAN's tool (INTANToNWB) requires so. See:
  https://github.com/Intan-Technologies/IntanToNWB
  With this wrapper, Python will never show up, neither the user will have to interact with it AT ALL. Everything happens from MATLAB.
  Requires Python installed in the machine. As explained in the Script:
  To date, MATLAB 2021b accepts interaction with Python up to v3.9. Needs the 64 bits version:
  (https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
  To check that MATLAB can interact with the Python Modules, look if variable 'pe' is correctly populated when debugging.

# Getting Started
For the very first time: Set the folder where 'NGL01_Main' script is as working folder.

Adding this folder permanently to Matlab's path could be useful, the rest will be temporarily added during the run.

Fill out the Inputs ('datadrive', 'studyName', 'subjects' and 'dates').
Fill out the Options at wish.

Hit F5.

# Script Description.
Pipeline process INTAN and Deuteron continous data.
Will read and process INTAN, DEUTERON (or ALLEGO) data, from selected sessions for a given animal.
The main pipeline will be: INTAN/DEUTERON raw formats to be located, then
converted to .bin files (spike sorting), and Fieldtrip .mat structures
(for LFP). Once sorted, spike data will be attached to the FieldTrip
structure. For Arena experiments, motion sensor data will be extracted and
interpreted. Data will be trial-parsed using EventCodes.

The hard disk data structure SHOULD fit the IKN standard published at:
gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure

DEPENDENCIES:
Requires that all pipeline dependencies are properly located. 
I suggest to include the 'mainfolder' in Matlab's permanent path system.
The function 'set_default' will take care of the rest of folders on each run.

INPUTS:
      input.datadrive, char array with the drive where data is located. As 'D:\'
      input.studyName, char array with the project name, matching the
                         folder name where all data will be stored. As 'studyName'
      input.subjects,  char array with either 'all' OR a single subject name e.g. 'DOE'
      input.dates,     char array with either 'all' OR a cell array of dates 
                         for a SINGLE subject e.g. {'YYYYMMDD' 'yyyymmdd' ...)

OPTIONS: is a struct with many possible fields. All should have a
corresponding default inside whatever function is being called. Main ones
are:     
    opt.bin,              Creation of .bin file, input to Kilosort 2/4.
    opt.FTfile,           Creation of .mat file with FieldTrip format.
    opt.RetrieveEvents,   Retrieve event log from Deuteron system.
    opt.GetMotionSensors, Retrieve data from motion sensors in Deuteron.
    opt.kilosort,         Asks to proceed with KS processing and waits to retrieve its results.
    opt.set_filter,       If Deuteron data was adquired with a wideband.
    opt.lowpass,          Lowpass band to extract LFP from wideband.
    opt.highpass,         Highpass band to extract spike activity.
    opt.h5,               Creation of .h5 file (not really used, to deprecate?).

OUTPUTS:
For one single session or for a batch of sessions, from one single animal:
       Fieldtrip (.mat), binary (.bin), HDF5 (.h5) and/or .nwb files from
           1. Deuteron .DT2 or .DF1 data.
           2. INTAN file-per-type and file-per-channel format data.
           3. (ALLEGO data?)
       EventRecord.mat file, from Deuteron session.
       MotionData.mat file, From Deuteron sensors.
       Plots snippets of time- and frequency-domain data, from FieldTrip
       
Last modified 05.04.2023 (Jesus)
