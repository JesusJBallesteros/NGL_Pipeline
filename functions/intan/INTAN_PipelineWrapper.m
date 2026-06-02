function [input, opt] = INTAN_PipelineWrapper(input, varargin)
% INTAN_PipelineWrapper. Concentrates all INTAN-specific preprocessing stages.
%
% PURPOSE:
%   Called from NGL01_Main for sessions in 'fileperch', 'filepertype', or
%   'tradFormat'. Coordinates event extraction, optional NWB export, binary
%   file creation for Kilosort, FieldTrip LFP conversion, and motion-sensor
%   extraction. Each sub-stage is gated by the corresponding opt flag and
%   skips with a warning if its output already exists on disk.
%
% USAGE:
%   [input, opt] = INTAN_PipelineWrapper(input, opt)
%   [input, opt] = INTAN_PipelineWrapper(input) % opt defaults used. % Unlikely
%
% INPUTS:
%   input  - struct (from set_default + prepforsession), contains input.run,
%              input.sessions(x).info, all path fields
%   opt    - options struct (from set_default); may be omitted (empty struct assumed)
%
% OUTPUTS:
%   input  - passed through (unchanged by this wrapper; sub-functions may
%              update input.sessions fields)
%   opt    - updated by EventProcess, adds opt.eventdef and opt.newEvent
%
% STAGES:
%   01. EventProcess            - extract events, build trialdef, run conditions_script
%   02. intan2NWB_neuroconv     - NWB export (only if opt.doNWB)
%   03. Intan2Kilosort_wrapper  - create .bin file (only if opt.bin and not cached)
%   04. intan2MAT_wrapper +
%       MAT2FieldTrip           - create FieldTrip LFP .mat (only if opt.FieldTrip)
%   05. GetMotionSensors        - extract accelerometer data (only if opt.GetMotionSensors)
%
% CALLS:
%   EventProcess, intan2NWB_neuroconv, Intan2Kilosort_wrapper,
%   intan2MAT_wrapper, MAT2FieldTrip, GetMotionSensors
%
% Version 06.05.2026 (Jesus)

if nargin < 2
    opt = struct();
    warning('opt variable not passed to INTAN wrapper. All Defaults will be used.')
elseif nargin == 2, opt = varargin{1};
end

%% 01. Event data retrieval and trial definition. INTAN version
% if isfile(fullfile(opt.trialSorted, "trialdef.mat")) % Removed ward for saved files
%     load(fullfile(opt.trialSorted, "trialdef.mat"))
%     if ~exist("trialdef","var") && exist("trialDefinition","var")
%         trialdef = trialDefinition.trl; clear trialDefinition
%     end
% else
    % 'trialdef' outputted for later feed into fieldtrip transf.
    [~, trialdef, ~, opt] = EventProcess(input, opt);
% end

%% 02. Create NWB file
if opt.doNWB % We want a .NWB file.
    % Run NeuroConv python app for the INTAN to NWB conversion:               
    % This NEEDS A PYTHON installation in the corresponding Conda Enviroment!
    % Detailed explanation:
    % WHAT IT IS: function to convert data from INTAN to .NWB format.
    % WHAT IT DOES: Checks for Python engine in computer. Checks given input 
    %  path to .rhs or .rhd file. Converts to output path .nwb file.
    %
    % Requires Python installed in the machine. Requires Neuroconv installed (in proper py enviroment).
    % To check access to Python Modules from MATLAB, look that 'pe' is correctly populated when running the script.
    intan2NWB_neuroconv(input, opt);  
else
    warning('Skipping NWB file creation.')
end 

%% 03. Run wrapper for the INTAN to Kilosort. Creates .bin and .h5 files
if opt.bin && ~isfile(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']))
    Intan2Kilosort_wrapper(input.sessions(input.run(1)), opt);
else
    warning('Skipping .bin file creation.')
end

%% 04. Run functions to convert INTAN dat to FIELDTRIP structure.
if opt.FieldTrip
    % Includes a mix of INTAN funtions. CREATES and GIVES proper
    % FieldTrip format without trial-parsing.
    INTANdata = intan2MAT_wrapper(input.sessions(input.run(1)), opt);
    MAT2FieldTrip(INTANdata, opt, trialdef, 1); %(data, options, trialdefinitions, force continuous)
end

%% 05. Get Motion Data into Matlab
if opt.GetMotionSensors
    disp('Extracting Motion Sensor data ...')
    GetMotionSensors(opt, input);
end
end