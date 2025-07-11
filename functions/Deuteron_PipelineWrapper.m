function [opt] = Deuteron_PipelineWrapper(input, varargin)
% Adaptation from the common pipeline for Deuteron. Wraps up the most common 
% processing lines necessary to get data from Deuteron raw files. This
% includes the Neural data and the motion sensors, so far. Could be
% expanded to extract audio as well.
%
% DEPENDENCIES
%   % Deuteron_EventFileReaderDll: To extract Event Record from Deuteron Block format.
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
%                   lowpass:    int array. upper boundary for lowpass filter. e.g. 200
%                   highpass:   int array. lower boundary for highpass filter. e.g. 450
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
% Version 12.06.2024 (Jesus)

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Default options.
if ~isfield(opt,'bin'),             opt.bin                 = true;         end
if ~isfield(opt,'FieldTrip'),       opt.FieldTrip           = true;         end
if ~isfield(opt,'RetrieveEvents'),  opt.RetrieveEvents      = true;         end
if ~isfield(opt,'GetMotionSensors'),opt.GetMotionSensors    = false;        end
if ~isfield(opt,'lowpass'),         opt.lowpass             = 150;          end
if ~isfield(opt,'highpass'),        opt.highpass            = 400;          end
if ~isfield(opt,'StpSz'),           opt.StpSz               = 1000000;      end
if ~isfield(opt,'parsetrial'),      opt.parsetrial          = false;        end
if ~isfield(opt,'CAR'),             opt.CAR                 = true;         end
if ~isfield(opt,'timebreak'),       opt.timebreak           = false;        end

%% Set local options.
% Collect parameters to proceed with file creation. List all files.
opt.myFiles = input.sessions(input.run(1)).info.files;
opt.ext     = input.sessions(input.run(1)).info.fileformat;

% Set Sample rate.
opt.sampleRate  = input.sessions(input.run(1)).info.amplifier_sample_rate;

% Few specific parameters from Deuteron's log and documentation, to convert bits to physical units.
opt.numberOfAdcBits   = input.sessions(input.run(1)).info.numADCBits;
opt.voltageResolution = input.sessions(input.run(1)).info.voltageRes;
opt.offset            = 2^(opt.numberOfAdcBits-1);

%% Event data retrieval and trial definition.
% 'trialdef' outputted for later feed into fieldtrip transf.
% An empty output means that data shall be treated as continuous.
[events, trialdef, EventRecord] = EventProcess(opt);

% Register a timebreak if detected
if isfield(EventRecord,'TimeBreak')
    if ~isempty(EventRecord.TimeBreak{1,2})
        opt.timebreak = true;
    end
end

%% High-pass Neural Data Conversion to .bin
if opt.bin
    % Converts Deuteron DT2/DF1 files into .bin and/or .h5 files.
    disp('Converting Deuteron files to .bin format.');
    Deuteron2Kilosort(opt);
end

%% Low-pass Neural Data Conversion to FT format.
if opt.FieldTrip
    % Convert Deuteron files into a FieldTrip formatted .mat file
    disp('Converting Deuteron files into a pseudo-FT file.');
    FT_data = []; % If left empty, proper FT formatting will be skipped.
    [FT_data] = Deuteron2Fieldtrip(opt);
    
    % Give proper Fieldtrip format and trial parse data.
    disp('Giving proper FieldTrip format.');
    MAT2FieldTrip(FT_data, opt, trialdef, true); 
end

%% Get Motion Data into Matlab
if opt.GetMotionSensors
    disp('Extracting Motion Sensor data from Deuteron...')
    Deuteron_GetMotionSensors(opt);
end

end