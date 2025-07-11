function input = set_default(input, opt)
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
% Jesus. 28.05.2024

%% Set default to extract data without NWB file creation.
% Due to a conflict at h5 python-matlab dlls, when the two following pipelines 
% are requested, the NWB will perform well but the data extraction will not. 
% It will crash for not completely known reason. It needs a Matlab restart between runs.
if ~isfield(opt,'doNWB') || isempty(opt.doNWB),           opt.doNWB      = false; end

%% Fix drive letter if needed.
if ~contains(input.datadrive,':\')
    input.datadrive = [input.datadrive ':\'];
end

%% Find toolbox
cd(input.toolbox)

%% Set default paths. IKN Standard recommended.
input.analysisCode  = fullfile(input.datadrive, input.studyName, '\analysisCode\');
input.datafolder    = fullfile(input.datadrive, input.studyName, '\data\raw\');              % Default: '\data\raw'
input.analysis      = fullfile(input.datadrive, input.studyName, '\data\analysis\');
input.bhvfolder     = fullfile(input.datadrive, input.studyName, '\data\behaviour\'); % deprecate?
input.spikeSorted   = fullfile(input.datadrive, input.studyName, '\data\spikeSorted\');
input.trialSorted   = fullfile(input.datadrive, input.studyName, '\data\trialSorted\');
input.processed     = fullfile(input.datadrive, input.studyName, '\data\preprocessing\');    % Default: '\data\preprocessing'

% When preprocessing as already be done, this could be overwritten later.
input.jump2plot = false;

%% Find requested subjects.
% In case is left empty or deleted, default to 'all'
if ~isfield(input,'subjects') || isempty(input.subjects)
    input.subjects = 'all';
end

% Get available subjects. Read all existing content under datafolder
cd(fullfile(input.datafolder))
    if isfile("_findatserver"), cd(fullfile(input.processed)), end
subjects = dir('???*');

if strcmp(input.subjects, 'all') % request is 'all'
    input.subjects = subjects; % Add them all to the list
else % Request is subset or already processed list
    if iscell(input.subjects)
        subjidx = ismember({subjects.name}, input.subjects); % Index those requested
        input.subjects = subjects(subjidx); % Add the indexed members
    else
        warning('Requested subjects and sessions as ran in NGL01. Re-run options if you want a different subset.')
    end
end

% Get final number of subjects added
if ~isfield(input,'nsubjects') || isempty(input.nsubjects)
    input.nsubjects = length(input.subjects);
end

%% If NWB requested, Python-based toolbox needed. 
if opt.doNWB 
    if ~isfield(input,'NCfolder') || isempty(input.NCfolder)
        input.NCfolder = 'C:\Code\miniconda3\envs\neuroconv\python.exe'; % Path to the neuroconv git-code. It could change among PCs
%         input.NCfolder = [input.toolbox '\toolboxes\IntanToNWB']; % Add it
    end
end

%% Kilosort-related
if isfield(opt,'kilosort')
    if opt.kilosort == 2
       input.KSpath = 'C:\Kilosort_2.0'; % Absolute path to kilosort, hardcoded. It could change among PCs.
       
       % Add Kilosort (external)
       addpath(genpath(input.KSpath)) % path to kilosort toolbox
    
    elseif opt.kilosort == 4
       input.KSpyfolder = 'C:\Code\miniconda3\envs\kilosort\Lib\site-packages\kilosort'; % Path to the kilosort git-code. It could change among PCs
%        input.KSpyenv_NGL = 'C:\Users\ballejd4\AppData\Local\miniconda3\envs\kilosort'; % Kilosort enviroment under python installation. It could change among PCs
       
       % Copy NGL customized py files to kilosort enviroment's library
       orig_dir = pwd; % where we come from
       cd(input.analysisCode) % go to where customized py files are, \analysisCode
            copyfile("*.py", input.KSpyfolder); % copy the customized py files to the working directory
       cd(orig_dir) % back to previous folder
    
       % No need to add to matlab path
    end
end

%% Set Dependencies. Critical to find toolboxes.
cd(input.toolbox)

% Add functions
addpath functions\
addpath(input.analysisCode)

% Add toolboxes
addpath toolboxes\Intan
addpath toolboxes\fieldtrip_light
addpath toolboxes\Viewer
addpath toolboxes\BDPAT_NGL
% Add toolboxes and all their subfolders
addpath(genpath('toolboxes\Deuteron'))
addpath(genpath('toolboxes\npy-matlab'))
addpath(genpath('toolboxes\bombcell'))
addpath(genpath('toolboxes\prettify_matlab'))
addpath(genpath('toolboxes\spikes'))

% Initialize FT
ft_defaults

% Hardcode Deuteron's exe/dll files location to RetrieveEvents.
input.ReaderDll   = [input.toolbox, '\toolboxes\Deuteron\software\Event_File_Reader_9_0.dll']; % hopefully wont be necessary
input.exefile     = [input.toolbox, '\toolboxes\Deuteron\software\Event_File_Reader_9_0.exe'];

end