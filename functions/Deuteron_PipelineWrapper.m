function Deuteron_PipelineWrapper(sessions, varargin)
% Adaptation from the common pipeline for Deuteron. Wraps up the most common 
% processing lines necessary to get data from Deuteron raw files. This
% includes the Neural data and the motion sensors, so far. Could be
% expanded to extract audio as well.
%
% DEPENDENCIES
%   Deuteron_EventFileReaderDll: To extract Event Record from Deuteron Block format.
%   Deuteron2Kilosort: To compile recorded data in a single file per channel.
%                      Can also split the data based on event codes. 
%   Deuteron_GetMotionSensors: To extract data from motion sensors.
%   Deuteron_PlotMotionSensors: To process and visualize data from motion sensors.
%
% INPUTS:
%    sessions: struct. Variable containing info about sessions in process
%    ss:       int. Current session ordinal in the pipeline
%    opt:      struct. optional inputs to override the defaults:
%                   h5:     logic. Creation of .h5 file. Normally 'false'
%                   bin:    logic. Creation of .bin file. Normally 'true'
%                   FTfile: logic. Creation of Fieltrip-formatted .mat file.
%                   RetrieveEvents:   logic. Retrieve eventlog from Deuteron (and extract eventcodes and timestamps from it).
%                   GetMotionSensors: logic. Extraction and processing of motion sensor data.
%                   lowpass:    int array. lower and upper boundaries for lowpass filter. [  0  300]
%                   highpass:   int array. lower and upper boundaries for highpass filter. [300 7500]
%                   DllFolder:  string. Location of the .dll file to process events in Deuteron.
%                   set_filter: logic. Filtering (and downsampling) request.
%                   StpSz:      int. Number of samples to be written per chunck.
%
% OUTPUS:
%    EventRecord.mat file and compressed events file.
%    .bin file, as channels x samples. If requested.
%    .h5 file, as channels x sample. If requested.
%       (both with channels in increasing order as required for Kilosort.)
%    MotionData.mat file, with [Accelerometer, Gyroscope, Magnetometer] variables
%       containing timeseries for each sensor readings, in physical units. Plus
%       a 'rotators' variable, containing the quaternions to create the
%       rotation matrices and other transformations.
%
% Version 06.03.2023 Jesus

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Options 
if ~isfield(opt,'bin'),             opt.bin                 = true;         end
if ~isfield(opt,'FTfile'),          opt.FTfile              = true;         end
if ~isfield(opt,'RetrieveEvents'),  opt.RetrieveEvents      = true;         end
if ~isfield(opt,'GetMotionSensors'),opt.GetMotionSensors    = false;        end

if ~isfield(opt,'set_filter'),      opt.set_filter          = 1;            end
if ~isfield(opt,'lowpass'),         opt.lowpass             = [  0  400];   end
if ~isfield(opt,'highpass'),        opt.highpass            = [500 7500];   end
if ~isfield(opt,'StpSz'),           opt.StpSz               = 1000000;      end

% Hardcode the .dll file from Deuteron. Not really an option.
opt.ReaderDll = 'C:\Code\Scripts\ephys-data-pipeline\functions\dlls\Event_File_Reader_8_3.dll';

%% Parameters
% Collect parameters to proceed with file creation. List all files.
opt.myFiles = sessions.info.files;
opt.ext     = sessions.info.fileformat;

% Sample rate.
opt.sampleRate  = sessions.info.sampleRate;

% ChunkSize of HDF5 file (e.g., 5 minutes is, 300s at 30000Hz = 9600000 samples)
%  this chunk size works well. optimal? Once it is, this variable no longer requires user input.
opt.HDF5chunkSize = 300*opt.sampleRate; 

% Get number of channels.
opt.numChannels     = sessions.info.numChannels;
opt.channelOrder    = 1:1:opt.numChannels; 

% We need this parameters from Deuteron's log and documentation, to convert 
% to physical units. (At least for .DT2)
opt.numberOfAdcBits   = sessions.info.numADCBits;
opt.voltageResolution = sessions.info.voltageRes;
opt.offset            = 2^(opt.numberOfAdcBits-1);

%% Event data, using dll
if opt.RetrieveEvents && strcmp(sessions.info.fileformat, 'DF1')
    disp('Retrieving Events from Deuteron BLOCK format.')
   
    % Proceed to extract all events during session.
    [EventRecord, sessions.info.numChannels] = ...
        Deuteron_EventFileReaderDll(opt, sessions);
    
else
   disp('Event extraction not requested or session is FLAT format. Skipping...')
end

%% Neural Data Conversion.
if opt.bin
    disp('Converting Deuteron files into .h5 and .bin...');

    % Converts Deuteron DT2 and DF1 files into single files (.bin and .h5) 
    % to further use (i.e. with Kilosort)
    Deuteron2Kilosort(opt);
end

if opt.FTfile
    % Convert Deuteron files into a FieldTrip formatted .mat file
    disp('Converting Deuteron files into a pseudo-FT file.');
    [DEUTdata] = Deuteron2mat(opt);

    %% Conversion to Fieldtrip formatted data
    disp('Giving proper FieldTrip format.');
    mat2FieldTrip(DEUTdata, opt);
end

%% Motion Data to Matlab
if opt.GetMotionSensors
    disp('Extracting Motion Sensor data from Deuteron...')
    Deuteron_GetMotionSensors(opt);
    
end

end