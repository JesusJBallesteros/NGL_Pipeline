function [input, opt] = Deuteron_PipelineWrapper(input, varargin)
% Deuteron_PipelineWrapper  Entry point for Deuteron raw data processing.
%
% PURPOSE:
%   Orchestrates the complete Deuteron processing chain for a single session:
%   event extraction → wideband .bin creation (Kilosort input) →
%   LFP pseudo-FieldTrip file → optional motion-sensor extraction.
%   Sets session-specific fields on opt (myFiles, ext, sampleRate, ADC
%   parameters) before delegating to specialised sub-functions.
%
% USAGE:
%   [input, opt] = Deuteron_PipelineWrapper(input)
%   [input, opt] = Deuteron_PipelineWrapper(input, opt)
%   Called from NGL01_Main inside the session loop.
%
% INPUTS:
%   input  - struct built by set_default and findSessions; must contain:
%              .sessions(x).info.files                  dir-struct of NEUR*.DF1 files
%              .sessions(x).info.fileformat             'DF1' (or legacy 'DT2')
%              .sessions(x).info.amplifier_sample_rate  Hz (typically 32000)
%              .sessions(x).info.numADCBits             ADC bit depth (16)
%              .sessions(x).info.voltageRes             V/bit scaling factor
%              .run                                     index into sessions being processed
%   opt    - (optional) options struct; all fields default via default_opt if omitted
%
% OUTPUTS:
%   input  - passed through unchanged
%   opt    - updated with session-specific fields:
%              .myFiles, .ext, .sampleRate, .numberOfAdcBits,
%              .voltageResolution, .offset, .channelOrder, .numChannels,
%              .eventdef, .timebreak
%
% CALLS:
%   EventProcess, Deuteron2Kilosort, Deuteron2Fieldtrip,
%   MAT2FieldTrip, GetMotionSensors
%
% Last modified 08.05.2026 (Jesus)

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

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
[events, trialdef, EventRecord, opt] = EventProcess(input, opt);

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
    Deuteron2Kilosort(opt)
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
    GetMotionSensors(opt, input);
end

end