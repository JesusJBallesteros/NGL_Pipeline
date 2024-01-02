function [sessions] = Deuteron_PipelineWrapper(input, varargin)
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
%                   lowpass:    int array. lower and upper boundaries for lowpass filter. e.g. [  0  300]
%                   highpass:   int array. lower and upper boundaries for highpass filter. e.g. [300 7500]
%                   DllFolder:  string. Location of the .dll file to process events in Deuteron.
%                   set_filter: logic. Filtering (and downsampling) request.
%                   StpSz:      int. Number of samples to be written per chunck.
%
% OUTPUTS:
%    EventRecord.mat file. If requested.
%    .bin file, as channels x samples. If requested.
%    .h5 file, as channels x sample. If requested.
%       (both with channels in increasing order as required for Kilosort.)
%    MotionData.mat file, with [Accelerometer, Gyroscope, Magnetometer] variables
%       containing timeseries for each sensor readings, in physical units. Plus
%       a 'rotators' variable, containing the quaternions to create the
%       rotation matrices and other transformations.
%
% Version 02.01.2024 (Jesus)

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Default options.
if ~isfield(opt,'bin'),             opt.bin                 = true;         end
if ~isfield(opt,'FTfile'),          opt.FTfile              = true;         end
if ~isfield(opt,'RetrieveEvents'),  opt.RetrieveEvents      = true;         end
if ~isfield(opt,'useexe'),          opt.useexe              = true;         end
if ~isfield(opt,'GetMotionSensors'),opt.GetMotionSensors    = false;        end
if ~isfield(opt,'set_filter'),      opt.set_filter          = 1;            end
if ~isfield(opt,'lowpass'),         opt.lowpass             = [  0  150];   end
if ~isfield(opt,'highpass'),        opt.highpass            = [450 5000];   end
if ~isfield(opt,'StpSz'),           opt.StpSz               = 1000000;      end
if ~isfield(opt,'useexe'),          opt.useexe              = true;         end
if ~isfield(opt,'usepar'),          opt.usepar              = false;        end
if ~isfield(opt,'parsetrial'),      opt.parsetrial          = true;         end

%% Parameters
    % exe/dll locations 
    opt.exefile     = input.exefile;
    opt.ReaderDll   = input.ReaderDll;

    % Collect parameters to proceed with file creation. List all files.
    opt.myFiles = input.sessions(input.run(1)).info.files;
    opt.ext     = input.sessions(input.run(1)).info.fileformat;
    
    % Sample rate.
    opt.sampleRate  = input.sessions(input.run(1)).info.amplifier_sample_rate;
    
    % ChunkSize of HDF5 file (e.g., 5 minutes: 300s x 30000Hz = 9600000 samples)
    opt.HDF5chunkSize = 300*opt.sampleRate; 
    
    % Get number of channels.
    if ~isempty(input.sessions(input.run(1)).info.nChannels)
        opt.numChannels     = input.sessions(input.run(1)).info.nChannels;
        opt.channelOrder    = 1:1:opt.numChannels; 
    else
        opt.numChannels     = [];
        opt.channelOrder    = []; 
    end

    % We need these parameters (from Deuteron's log and documentation), to convert to physical units.
    opt.numberOfAdcBits   = input.sessions(input.run(1)).info.numADCBits;
    opt.voltageResolution = input.sessions(input.run(1)).info.voltageRes;
    opt.offset            = 2^(opt.numberOfAdcBits-1);

%% Event data retrieval
if opt.RetrieveEvents && opt.useexe
    % Proceed to extract all events during session.
    disp('Retrieving Events from Deuteron.')
    [EventRecord, opt] = Deuteron_EventFileReaderDll(opt);
    
    % Create trial definition using the proper eventcodes
    % TODO
    trialdef = [];

elseif opt.RetrieveEvents && opt.usepar
    trialdef = [];
    
    % Proceed to define trials based on .par file
    % To discontinue once above works.
    if opt.parsetrial
        opt.def = [];
            % Set EventCode meaning
            opt.def.itiOn       = 1;
            opt.def.stimOn      = 2;
            opt.def.rwd         = 5;
            opt.def.tutor       = 4;
            opt.def.trialEnd    = 7;
            
            % Set event for trial t=0
            opt.def.alignto = opt.def.itiOn;
        
        if opt.usepar
            [trialdef, opt.def] = trialdefFromPar(opt);
        else
           % trialdef = trialdefFromEventRecord(def, input, sessions);
        end
    end
else
   disp('Events not requested. Skipping...')
end

%% High-pass Neural Data Conversion to .bin
if opt.bin
    % Converts Deuteron DT2/DF1 files into .bin and/or .h5 files.
    disp('Converting Deuteron files to .bin format...');
    Deuteron2Kilosort(opt);
end

%% Low pass Neural Data Conversion to FT format
if opt.FTfile
    if ~isfile(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_continous_FT.mat')))
        % Convert Deuteron files into a FieldTrip formatted .mat file
        disp('Converting Deuteron files into a pseudo-FT file.');
        [FT_data] = Deuteron2mat(opt);
    else
        load(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_continous_FT.mat')))
    end
    
    % Give proper Fieldtrip format and trial parse data.
    disp('Giving proper FieldTrip format.');
    MAT2FieldTrip(FT_data, opt, []); %trialdef after fixed
end

%% Motion Data to Matlab
if opt.GetMotionSensors
    disp('Extracting Motion Sensor data from Deuteron...')
    Deuteron_GetMotionSensors(opt);
    
end

end