function [input, opt] = set_default(input, opt)
% set_default  Merge options with defaults, validate inputs, resolve paths.
%
% PURPOSE:
%   validation and preparation of options for the NGL toolbox.
%   Must be called once, at the top of every pipeline script (NGL01_Main,
%   NGL02_postPhy, etc.), before any session loop begins.
%   After this call, 'opt' should be complete and validated, and
%   'input' carries all resolved paths, the subject list, and Python env
%   paths needed by downstream functions.
%
% USAGE:
%   [input, opt] = set_default(input, opt)
%
% INPUTS:
%   input - struct with at minimum:
%             .datadrive  (char) drive letter, e.g. 'D'
%             .studyName  (char) project folder name
%             .subjects   (char 'all' | cell of IDs)
%             .dates      (char 'all' | cell of YYYYMMDD strings)
%           Optional multi-area field:
%             .Areas      {1xN} cell of area labels, one per kcoords group
%   opt   - struct of user-specified options (can be empty struct)
%
% OUTPUTS:
%   input - struct with resolved IKN path fields:
%             .analysisCode, .datafolder, .analysis, .bhvfolder,
%             .spikeSorted, .trialSorted, .processed, .toolbox,
%             .KSpythonExe, .KSpyfolder, .PHYpythonExe, .PHYpyfolder
%             .ReaderDll, .exefile
%             .subjects   (dir-struct array)
%             .nsubjects  (scalar count)
%           If input.Areas was set:
%             .areaMap    struct from buildAreaMap (see buildAreaMap.m)
%   opt   - complete options struct; all fields from default_opt present,
%           user values override defaults, unknown fields warned and dropped.
%
% SECTIONS:
%   1. Merge user opt with defaults (default_opt)
%   2. Validate and normalise fields
%   3. Resolve dependencies (warn on conflicts)
%   4. Validate and normalise input struct
%   4b. Build area map if input.Areas is defined (multi-area)
%   5. Build IKN standard path structure
%   6. Resolve subject list
%   7. Add all dependency folders to MATLAB path, call ft_defaults
%
% REQUIRES:
%   default_opt.m, NGL_machineConfig.m (in analysisCode/)
%
% Last modified 13.05.2026 (Jesus)

%% SECTION 1: Merge user-provided opt with canonical defaults
defaults   = default_opt();
userFields = fieldnames(opt);

for i = 1:numel(userFields)
    f = userFields{i};
    if isfield(defaults, f)
        defaults.(f) = opt.(f);
    else
        warning('NGL:unknownOption', 'Option ''%s'' is not recognised by default_opt() and will be ignored.', f);
    end
end
opt = defaults;

%% SECTION 2: Validate and normalise opt fields
if ischar(opt.alignto)
    opt.alignto = {opt.alignto};
end
assert(iscell(opt.alignto) && all(cellfun(@ischar, opt.alignto)), ...
    'NGL:invalidOption', ...
    'opt.alignto must be a char array or cell array of char arrays (e.g. {''itiOn'', ''rwd''}).');

assert(isnumeric(opt.numChannels) && isscalar(opt.numChannels) && opt.numChannels > 0, ...
    'NGL:invalidOption', 'opt.numChannels must be a positive scalar integer.');

if ~isempty(opt.lowpass)
    assert(isnumeric(opt.lowpass) && isscalar(opt.lowpass) && ...
           opt.lowpass > 0 && opt.lowpass < 9500, ...
        'NGL:invalidOption', ...
        'opt.lowpass must be a scalar between 0 and 9500 Hz, or empty (= off).');
end

assert(islogical(opt.kilosort) && isscalar(opt.lowpass), 'NGL:invalidOption', 'opt.kilosort must be true or false.');

%% SECTION 3: Resolve dependencies
if opt.phy && ~opt.bombcell
    warning('NGL:phyWithoutBombcell', ...
        ['opt.phy=true but opt.bombcell=false. ', ...
         'Manual curation will proceed without Bombcell screening.']);
end

if ~opt.RetrieveEvents && numel(opt.alignto) > 1
    warning('NGL:eventsMismatch', ...
        'opt.RetrieveEvents=false but multiple alignto events were specified.');
end

%% SECTION 4: Validate and normalise input struct
if ~contains(input.datadrive, ':\')
    input.datadrive = [upper(input.datadrive(1)) ':\'];
end
assert(ischar(input.datadrive) && length(input.datadrive) == 3, ...
    'NGL:invalidInput', ...
    'input.datadrive must be a single drive letter (e.g. ''D'').');

assert(ischar(input.studyName) && ~isempty(strtrim(input.studyName)), ...
    'NGL:invalidInput', 'input.studyName must be a non-empty character array.');

%% SECTION 5: Build standard path structure
base = fullfile(input.datadrive, input.studyName);

input.analysisCode  = fullfile(base, 'analysisCode');
input.datafolder    = fullfile(base, 'data', 'raw');
input.analysis      = fullfile(base, 'data', 'analysis');
input.bhvfolder     = fullfile(base, 'data', 'behaviour');
input.spikeSorted   = fullfile(base, 'data', 'spikeSorted');
input.trialSorted   = fullfile(base, 'data', 'trialSorted');
input.processed     = fullfile(base, 'data', 'preprocessing');

addpath(input.analysisCode)
cfg = NGL_machineConfig();

input.toolbox      = cfg.toolbox;
input.KSpythonExe  = cfg.KSpythonExe;
input.KSpyfolder   = fullfile(cfg.KSpythonExe, 'Lib', 'site-packages', 'kilosort');
input.PHYpythonExe = cfg.PHYpythonExe;
input.PHYpyfolder  = fullfile(cfg.PHYpythonExe, 'Lib', 'site-packages', 'phy');

input.ReaderDll = fullfile(input.toolbox, 'toolboxes', 'Deuteron', 'software', 'Event_File_Reader_9_0.dll');
input.exefile   = fullfile(input.toolbox, 'toolboxes', 'Deuteron', 'software', 'Event_File_Reader_9_0.exe');

if opt.doNWB
    input.NCfolder = cfg.NCpythonExe;
end

%% SECTION 6: Resolve subject list
if ~isfield(input, 'subjects') || isempty(input.subjects)
    input.subjects = 'all';
end

searchRoot = input.datafolder;
if isfile(fullfile(searchRoot, '_findatserver'))
    searchRoot = input.processed;
end

cd(searchRoot)
available = dir('???*');

if strcmp(input.subjects, 'all')
    input.subjects = available;
elseif iscell(input.subjects)
    idx = ismember({available.name}, input.subjects);
    input.subjects = available(idx);
else
    warning('NGL:subjectFormat', ...
        ['input.subjects was neither ''all'' nor a cell array. ', ...
         'Assuming it was set from a previous NGL01 run. Re-assign to change the subset.']);
end

input.nsubjects = numel(input.subjects);

%% SECTION 7: Add dependencies to MATLAB path
cd(input.toolbox)

% functions/ subfolders — one addpath per logical group.
% Keeping them explicit (rather than genpath) avoids pulling in _deprecated/.
addpath(fullfile('functions', 'pipeline'))
addpath(fullfile('functions', 'intan'))
addpath(fullfile('functions', 'deuteron'))
addpath(fullfile('functions', 'sorting'))
addpath(fullfile('functions', 'events'))
addpath(fullfile('functions', 'analysis'))
addpath(fullfile('functions', 'video'))
addpath(fullfile('functions', 'plotting'))
addpath(fullfile('functions', 'ethology'))
addpath(fullfile('functions', 'utils'))

addpath(fullfile('toolboxes', 'Intan'))
addpath(fullfile('toolboxes', 'fieldtrip_light'))
addpath(fullfile('toolboxes', 'Viewer'))
addpath(fullfile('toolboxes', 'BDPAT_NGL'))
addpath(genpath(fullfile('toolboxes', 'Deuteron')))
addpath(genpath(fullfile('toolboxes', 'npy-matlab')))
addpath(genpath(fullfile('toolboxes', 'matnwb')))
addpath(genpath(fullfile('toolboxes', 'bombcell')))
addpath(genpath(fullfile('toolboxes', 'prettify_matlab')))
addpath(genpath(fullfile('toolboxes', 'spikes')))

ft_defaults

% Validate that all required config files are present in analysisCode\.
% Raises an error listing every gap at once so nothing is silently missing.
checkAnalysisCode(input, opt)

disp('Defaults and User Options successfully merged.')

%% 8: Build area map if input.Areas is defined (multi-area)
% When the user supplies input.Areas (e.g. {'NCL','NCL','STR'}), each entry
% labels the kcoords group at that index in the chanMap. buildAreaMap derives
% per-area channel masks and writes per-area chanMap .mat files to analysisCode.
% If input.Areas is absent this section is skipped and the pipeline runs in
% backward-compatible single-area mode.
if isfield(input, 'Areas') && ~isempty(input.Areas)
    assert(~isempty(opt.KSchanMapFile), ...
        'NGL:areaMap', ...
        'input.Areas is set but opt.KSchanMapFile is empty. Provide a chanMap filename.');
    chanMapFullPath = fullfile(input.analysisCode, opt.KSchanMapFile);
    input.areaMap   = buildAreaMap(input.Areas, chanMapFullPath);
    fprintf('Multi-area mode: %d unique areas (%s)\n', ...
        numel(input.areaMap.uniqueAreas), strjoin(input.areaMap.uniqueAreas, ', '));
end
end
