function input = set_default(input)
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
% Jesus. 21.12.2023

% Set default to extract data without NWB file creation.
% Due to a conflict at h5 python-matlab dlls, when the two following pipelines 
% are requested, the NWB will perform well but the data extraction will not. 
% It will crash for not completely known reason. It needs a Matlab restart between runs.
if ~isfield(input,'useNWB') || isempty(input.useNWB),           input.useNWB      = false; end
% MEANING: do not set both 'true'.

% input.datadrive = [input.datadrive ':\'];

%% Set default paths. IKN Standard recommended.
input.datafolder    = fullfile(input.datadrive, input.studyName, '\data\raw\');              % Default: '\data\raw'
input.processed     = fullfile(input.datadrive, input.studyName, '\data\preprocessing\');    % Default: '\data\preprocessing'
input.analysisCode  = fullfile(input.datadrive, input.studyName, '\analysisCode\');
input.sorted        = fullfile(input.datadrive, input.studyName, '\data\spikeSorted\');
input.analysisData  = fullfile(input.datadrive, input.studyName, '\data\analysis\');
input.bhvfolder     = fullfile(input.datadrive, input.studyName, '\data\behavior\');

%% Find requested subjects.
% In case is left empty or deleted, default to 'all'
if ~isfield(input,'subjects') || isempty(input.subjects)
    input.subjects = 'all';
end

% Get subjects. Read all existing content under datafolder
cd(fullfile(input.datafolder))
subjects = dir('???*');

if strcmp(input.subjects, 'all') % request is 'all'
    input.subjects = subjects; % Add them all to the list
else % Request is subset
    subjidx = ismember({subjects.name}, input.subjects); % Index those requested
    input.subjects = subjects(subjidx); % Add the indexed members
end

% Get final number of subjects added
input.nsubjects = length(input.subjects);

%% Optional Inputs    
% If NWB requested, Python-based toolbox needed. 
if input.useNWB 
    if ~isfield(input,'pyfolder') || isempty(input.pyfolder)
        input.pyfolder = [input.toolbox '\toolboxes\IntanToNWB']; % Add it
    end
end

% Plots. Normally left empty.
if ~isfield(input,'plots'),   input.plots     = []; end
if ~isfield(input,'test_ch'), input.test_ch   = []; end

%% Set Dependencies. Critical to find toolboxes.
cd(input.toolbox)

% Add functions
addpath functions\

% Add toolboxes
addpath toolboxes\Intan
addpath toolboxes\fieldtrip_light
addpath toolboxes\Viewer
addpath(genpath('toolboxes\Deuteron'))
addpath(genpath('toolboxes\npy-matlab'))
addpath(genpath('toolboxes\bombcell'))
% %addpath(genpath('toolboxes\spikes'))
% % addpath(genpath('toolboxes\multitaper_prerau'))

% Add Kilosort (external)
addpath(genpath('C:\KiloSort_2.0\')) % path to kilosort toolbox (Assumes Sorting PC, not local)

% Initialize FT
ft_defaults

% Hardcode Deuteron's exe/dll files location to RetrieveEvents.
input.ReaderDll   = [input.toolbox, '\toolboxes\Deuteron\software\Event_File_Reader_9_0.dll'];
input.exefile     = [input.toolbox, '\functions\dlls\EventFileReader\Event_File_Reader_9_0.exe'];

end