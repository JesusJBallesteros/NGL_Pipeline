function set_default(input)
% 'set_default' adds the dependencies, included under the main folder.
%
% It reads the inputs, if any, and validates them.
%
% If not valid values are provided, it ask for the necessary ones and puts 
% them in correct format. If any necessary input continues to be invalid or
% empty, throws error.
% 
% For missing optionals, it uses defaults.
%
% Jesus. 04.04.2023

% Set default to Extract data withour NWB file creation.
% Due to a conflict at h5 python-matlab dlls, when the two following pipelines 
% are requested, the NWB will perform well but the data extraction will not. 
% It will crash for not completely known reason. It needs a Matlab restart between runs.
if ~isfield(input,'ExtractData') || isempty(input.ExtractData), input.ExtractData = true; end
if ~isfield(input,'useNWB') || isempty(input.useNWB),           input.useNWB      = false;end
% Meaning, do not run both 'true' (for now). 

%% Set default paths. IKN Standard recommended.
input.toolbox    = 'C:\Code\Scripts\ephys-data-pipeline'; % Default: 'C:\Code\Scripts\ephys-data-pipeline'
input.datafolder = fullfile(input.datadrive, input.studyName, '\data\raw\');              % Default: '\data\raw'
input.processed  = fullfile(input.datadrive, input.studyName, '\data\preprocessing\');    % Default: '\data\preprocessing'

%% Find requested subjects.
if ~isfield(input,'subjects') || isempty(input.subjects)
    input.subjects = 'all';
end

% Get subjects
cd(fullfile(input.datafolder))
subjects = dir('???*');

if strcmp(input.subjects, 'all')
    input.subjects = subjects;
else
    subjidx = ismember({subjects.name}, input.subjects);
    input.subjects = subjects(subjidx); 
end
input.nsubjects = length(input.subjects);

%% Optional Inputs    
% If NWB requested, Python-based toolbox needed.
if input.useNWB 
    if ~isfield(input,'pyfolder') || isempty(input.pyfolder)
        input.pyfolder = 'C:\Code\Scripts\ephys-data-pipeline\toolboxes\IntanToNWB';
    end
end

% Plots
if ~isfield(input,'plots'),   input.plots     = []; end
if ~isfield(input,'test_ch'), input.test_ch   = []; end

%% Set Dependencies
cd(input.toolbox)
addpath functions\
addpath toolboxes\Deuteron
addpath toolboxes\Intan
addpath toolboxes\fieldtrip_light
addpath toolboxes\Viewer

addpath(genpath('toolboxes\spikes'))
addpath(genpath('toolboxes\npy-matlab'))
addpath(genpath('toolboxes\multitaper_prerau'))
ft_defaults

%% Send input to base workspace
assignin('base','input', input);
end