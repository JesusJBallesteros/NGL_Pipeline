function [input, opt] = set_default(input, opt)
% set_default  Merge options with defaults, validate inputs, set paths and dependencies.
%   [input, opt] = set_default(input, opt)
%
%   This is the single validation and preparation gateway for the NGL toolbox.
%   It should be called once, early in each top-level pipeline script.
%   After this call:
%       - 'opt'   is guaranteed to be complete, validated, and consistent.
%       - 'input' carries all resolved paths and the subject list.
%
%   Jesus. 27.03.2026

%% SECTION 1: Merge user-provided opt with canonical defaults
% Pull the full default set from the single source of truth, then overwrite
% only the fields the user actually specified. Fields left unset by the user
% silently receive their safe default. Unrecognised fields trigger a warning
% so typos surface immediately rather than being silently ignored.
defaults   = default_opt();
userFields = fieldnames(opt);

% Should exist saved under your analysisCode\ folder (together with NGL_SetAndRunMe.m)
% Should match the specifics of the machine your intend to use
cfg = NGL_machineConfig();

for i = 1:numel(userFields)
    f = userFields{i};
    if isfield(defaults, f)
        defaults.(f) = opt.(f);           % user value wins
    else
        warning('NGL:unknownOption', ...
            'Option ''%s'' is not recognised by default_opt() and will be ignored.', f);
    end
end
opt = defaults; % opt is now guaranteed to be complete

%% SECTION 2: Validate and normalise opt fields
% Catch type errors and common mistakes early, with clear messages,
% rather than letting them surface as cryptic errors deep in a pipeline.

% alignto: forgive the very common mistake of passing a plain char
if ischar(opt.alignto)
    opt.alignto = {opt.alignto};
end
assert(iscell(opt.alignto) && all(cellfun(@ischar, opt.alignto)), ...
    'NGL:invalidOption', ...
    'opt.alignto must be a char array or cell array of char arrays (e.g. {''itiOn'', ''rwd''}).');

% numChannels
assert(isnumeric(opt.numChannels) && isscalar(opt.numChannels) && opt.numChannels > 0, ...
    'NGL:invalidOption', 'opt.numChannels must be a positive scalar integer.');

% -- lowpass: must be a sensible frequency or empty (= off)
if ~isempty(opt.lowpass)
    assert(isnumeric(opt.lowpass) && isscalar(opt.lowpass) && ...
           opt.lowpass > 0 && opt.lowpass < 9500, ...
        'NGL:invalidOption', ...
        'opt.lowpass must be a scalar between 0 and 9500 Hz, or empty (= off).');
end

% Kilosort version must be 2 or 4
assert(ismember(opt.kilosort, [2, 4]), ...
    'NGL:invalidOption', 'opt.kilosort must be 2 or 4.');

%% SECTION 3: Resolve inter-option dependencies
% These are logical consistency checks across pairs of options.

% NWB and H5 pipelines conflict due to a DLL collision; warn clearly.
% A MATLAB restart is required between the two.
if opt.doNWB
    warning('NGL:nwbConflict', ...
        ['opt.doNWB=true: the NWB and H5 pipelines share conflicting DLLs. ', ...
         'If both are requested in sequence, restart MATLAB between runs.']);
end

% Manual Phy curation without a prior Bombcell pass loses semi-automation.
if opt.phy && ~opt.bombcell
    warning('NGL:phyWithoutBombcell', ...
        ['opt.phy=true but opt.bombcell=false. ', ...
         'Manual curation will proceed without Bombcell QC pre-filtering.']);
end

% Asking for multiple alignment events without extracting events is contradictory.
if ~opt.RetrieveEvents && numel(opt.alignto) > 1
    warning('NGL:eventsMismatch', ...
        'opt.RetrieveEvents=false but multiple alignto events were specified.');
end

%% SECTION 4: Validate and normalise input struct
% Drive letter: normalise to 'X:\' format
if ~contains(input.datadrive, ':\')
    input.datadrive = [upper(input.datadrive(1)) ':\'];
end
assert(ischar(input.datadrive) && length(input.datadrive) == 3, ...
    'NGL:invalidInput', ...
    'input.datadrive must be a single drive letter (e.g. ''D'').');

% Study name must be a non-empty char
assert(ischar(input.studyName) && ~isempty(strtrim(input.studyName)), ...
    'NGL:invalidInput', 'input.studyName must be a non-empty character array.');

%% SECTION 5: Build standard IKN path structure
base = fullfile(input.datadrive, input.studyName);

input.analysisCode  = fullfile(base, 'analysisCode');
input.datafolder    = fullfile(base, 'data', 'raw');
input.analysis      = fullfile(base, 'data', 'analysis');
input.bhvfolder     = fullfile(base, 'data', 'behaviour');   % candidate for deprecation
input.spikeSorted   = fullfile(base, 'data', 'spikeSorted');
input.trialSorted   = fullfile(base, 'data', 'trialSorted');
input.processed     = fullfile(base, 'data', 'preprocessing');

% Toolbox and external tool paths (from cfg)
input.toolbox     = cfg.toolbox;
input.KSpythonExe = cfg.KSpythonExe;
input.KSpyfolder  = fullfile(cfg.KSpythonExe, 'Lib', 'site-packages', 'kilosort');
input.PHYpythonExe = cfg.PHYpythonExe;
input.PHYpyfolder  = fullfile(cfg.PHYpythonExe, 'Lib', 'site-packages', 'phy');

% Deuteron reader binaries
input.ReaderDll = fullfile(input.toolbox, 'toolboxes', 'Deuteron', 'software', ...
    'Event_File_Reader_9_0.dll');
input.exefile   = fullfile(input.toolbox, 'toolboxes', 'Deuteron', 'software', ...
    'Event_File_Reader_9_0.exe');

% NWB/NeuroConv path (only needed if requested)
if opt.doNWB
    input.NCfolder = cfg.NCpythonExe;
end

%% SECTION 6: Resolve subject list
% Default to 'all' if subjects was left empty or the field is missing
if ~isfield(input, 'subjects') || isempty(input.subjects)
    input.subjects = 'all';
end

% Determine the search root: jump to 'preprocessing' if a server flag exists
searchRoot = input.datafolder;
if isfile(fullfile(searchRoot, '_findatserver'))
    searchRoot = input.processed;
end

% Read all available subject folders (3+ character names)
cd(searchRoot)
available = dir('???*');

if strcmp(input.subjects, 'all')
    input.subjects = available;
elseif iscell(input.subjects)
    idx = ismember({available.name}, input.subjects);
    input.subjects = available(idx);
else
    % Assume subjects was already a processed struct from a previous run
    warning('NGL:subjectFormat', ...
        ['input.subjects was neither ''all'' nor a cell array. ', ...
         'Assuming it was set from a previous NGL01 run. Re-assign to change the subset.']);
end

% Cache the count for loop bounds in pipeline scripts
input.nsubjects = numel(input.subjects);

%% SECTION 7: Add dependencies to MATLAB path
cd(input.toolbox)

addpath('functions')
addpath(input.analysisCode)

addpath(fullfile('toolboxes', 'Intan'))
addpath(fullfile('toolboxes', 'fieldtrip_light'))
addpath(fullfile('toolboxes', 'Viewer'))
addpath(fullfile('toolboxes', 'BDPAT_NGL'))
addpath(genpath(fullfile('toolboxes', 'Deuteron')))
addpath(genpath(fullfile('toolboxes', 'npy-matlab')))
addpath(genpath(fullfile('toolboxes', 'bombcell')))
addpath(genpath(fullfile('toolboxes', 'prettify_matlab')))
addpath(genpath(fullfile('toolboxes', 'spikes')))

ft_defaults   % initialise FieldTrip

end
