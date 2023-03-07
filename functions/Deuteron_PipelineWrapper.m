function Deuteron_PipelineWrapper(sessions, ss, varargin)
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
%    in:        struct. optional inputs to override the defaults:
%                   RetrieveEvents: logic. possibility to load event codes to build a restriced matrix
%                                   (e.g., the matrix starts at the first 'itiON' and ends at 'end'
%                                    experiment, removing paradigm irrelevant periods) 
%                   StpSz: int. relative to HDF5file: chunks in which ...  
%
% GENERATES:
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

if nargin < 3, opt = struct();
elseif nargin == 3, opt = varargin{1};
end

%% Options 
if ~isfield(opt,'RetrieveEvents'),  opt.RetrieveEvents      = true;         end
if ~isfield(opt,'StpSz'),           opt.StpSz               = 1000000;      end
if ~isfield(opt,'h5'),              opt.h5                  = true;         end
if ~isfield(opt,'bin'),             opt.bin                 = true;         end
if ~isfield(opt,'GetMotionSensors'),opt.GetMotionSensors    = true;         end
if ~isfield(opt,'FTfile'),          opt.FTfile              = true;         end
if ~isfield(opt,'lowpass'),         opt.lowpass             = [  0  300];   end
if ~isfield(opt,'highpass'),        opt.highpass            = [300 7500];   end
if ~isfield(opt,'DllFolder'),       opt.DllFolder           = 'C:\Code\Scripts\ephys-data-pipeline\functions\dlls'; end
if ~isfield(opt,'ReaderDll'),       opt.ReaderDll           = fullfile(opt.DllFolder, 'Event_File_Reader_8_3.dll'); end

if ~isfield(opt,'set_filter'),      opt.set_filter          = 0;            end

%% Parameters
% Collect parameters to proceed with file creation. List all files.
opt.myFiles = sessions.info{ss}.files;
opt.ext     = sessions.info{ss}.fileformat;

% Sample rate.
opt.sampleRate  = sessions.info{ss}.sampleRate;

% ChunkSize of HDF5 file (e.g., 5 minutes is, 300s at 30000Hz = 9600000 samples)
%  this chunk size works well. optimal? Once it is, this variable no longer requires user input.
opt.HDF5chunkSize = 300*opt.sampleRate; 

% Kilosort rearranges the rows of the input matrix according to the a channel map (which is developed in another file).
% Therefore, the matrix should be compiled with the channels in an increasing order.
opt.numChannels     = sessions.info{ss}.numChannels;
opt.channelOrder    = 1:1:opt.numChannels; 

% We need this parameters, obtained from Deuteron log and from the
% documentation, to convert to physical units. (At least for .DT2)
opt.numberOfAdcBits   = sessions.info{ss}.numADCBits;
opt.voltageResolution = 1.95e-7;
opt.offset            = 2^(opt.numberOfAdcBits-1);

%% Event data, using dll
if opt.RetrieveEvents
    % Proceed to extract all events during session. Give some feedback.
    disp('Retrieving Events from Deuteron...')
    [EventRecord, sessions.info{ss}.numChannels] = ...
        Deuteron_EventFileReaderDll(opt, sessions, ss);
    disp(EventRecord);
    disp(['Found ', int2str(sessions.info{ss}.numChannels), ' channels']);
else
   disp('Event extraction not requested, skipping...')
end

%% Neural Data Conversion.
if opt.h5 || opt.bin
    % Converts Deuteron DT2 and DF1 files into single files (.bin and .h5) 
    % to further use (i.e. with Kilosort)
    disp('Converting Deuteron files into .h5 and .bin...');
    Deuteron2KilosortV2(opt);
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
    [Accelerometer, Gyroscope, Magnetometer] = ...
        Deuteron_GetMotionSensors(opt);
    
    disp('Processing and Plotting Motion Sensor data.')
    % Add (..., 1, 1) to input, if visualization and video recording are wanted.
    Deuteron_PlotMotionSensors(Accelerometer, Gyroscope, Magnetometer, [], [])
end

end