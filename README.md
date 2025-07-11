# 'Ephys-data-pipeline'
Scripts, functions and tools to work with electrophysiological data at NGL.

# Use 'NGL00_Prep' to create a new project folder system to start storing your raw data and set initial readme info.
Then, drop your raw data into individual subject and session folders under 'data/raw'.
Your data SHOULD be stored as the IKN standard Harddisk data structure. 
See: gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure

# 'NGL01_Main' will transform raw data from INTAN and Deuteron into .bin files (for kilosort) and .mat (for Fieldtrip) files.
A set of options let the user decide wich pipelines to follow, and to specify any filters, thresholds and so on for the related steps.

This Script will process high-pass data and proceed to Kilosort it with no GUI. 
Inmediately after, it will call Bombcell to 'pre-curate' and create an initial set of tags for the sorted clusters.
Then, the user needs to manually curate the results. There is no way around this.

For low-pass data, the downsampled time series will be stored into .mat files with the FieldTrip expected format. 
Events will be used to trial-parse the data (or let it be continous) and give proper format to allow the use of FT functions.

# 'NGL02_postPhy' will proceed with typical steps to transform the manually-curated spike data to NLG data format.
It will read and extract data from the python-based files into MATLAB, generating spike matices according to the lab format.
This can then be feeded into further functions to analyze, plot, etc.
It will also process the spike data to fit the FieldTrip structures together with the LFP data, and trial parsed if required.
This would allow for spike-field analysis, as well as the use of FT funtions on both domains.

# Getting Started
For the very first time: Set the main repo folder (where 'NGL01_Main.m' is) as working folder.
*Adding this folder permanently to Matlab's path is recommended, the rest will be taken care of during the run.

Fill out the A, B and C Sections as needed. A is mandatory.

Hit F5.

# Script Description. (in progress)
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
       
Last modified 21.12.2023 (Jesus)
