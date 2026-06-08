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
% Last modified 11.05.2026 (Jesus)

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

assert(islogical(opt.kilosort), 'NGL:invalidOption', 'opt.kilosort must be true or false.');

% Waveform extraction (NGL02 / loadSpikes).
assert(islogical(opt.getwF) && isscalar(opt.getwF), ...
    'NGL:invalidOption', 'opt.getwF must be true or false.');
assert(isstruct(opt.gwfparams), ...
    'NGL:invalidOption', 'opt.gwfparams must be a struct.');
assert(isfield(opt.gwfparams,'wfWin') && isnumeric(opt.gwfparams.wfWin) && ...
       numel(opt.gwfparams.wfWin) == 2 && ...
       opt.gwfparams.wfWin(1) <= 0 && opt.gwfparams.wfWin(2) > 0, ...
    'NGL:invalidOption', ...
    'opt.gwfparams.wfWin must be a 2-element numeric [pre post] with pre<=0<post (samples).');
assert(isfield(opt.gwfparams,'nWf') && isnumeric(opt.gwfparams.nWf) && ...
       isscalar(opt.gwfparams.nWf) && opt.gwfparams.nWf > 0 && ...
       opt.gwfparams.nWf == floor(opt.gwfparams.nWf), ...
    'NGL:invalidOption', ...
    'opt.gwfparams.nWf must be a positive integer (max waveforms per cluster).');

% Cluster-loading and ISI-binning knobs (NGL02 / loadSpikes / calc_isihist).
assert(isstruct(opt.spparams) && ...
       isfield(opt.spparams,'excludeNoise') && islogical(opt.spparams.excludeNoise) && ...
       isfield(opt.spparams,'loadPCs')      && islogical(opt.spparams.loadPCs), ...
    'NGL:invalidOption', ...
    'opt.spparams must be a struct with logical fields .excludeNoise and .loadPCs.');
assert(isnumeric(opt.isibins) && isvector(opt.isibins) && numel(opt.isibins) >= 2 && ...
       all(diff(opt.isibins) > 0), ...
    'NGL:invalidOption', ...
    'opt.isibins must be a strictly increasing numeric vector of bin edges (ms).');

% Firing-rate binning (canonical values for popDyn smoothing).
assert(isnumeric(opt.binSize_ms) && isscalar(opt.binSize_ms) && opt.binSize_ms > 0, ...
    'NGL:invalidOption', 'opt.binSize_ms must be a positive scalar (ms).');
assert(isnumeric(opt.stepSz_ms) && isscalar(opt.stepSz_ms) && opt.stepSz_ms > 0, ...
    'NGL:invalidOption', 'opt.stepSz_ms must be a positive scalar (ms).');
assert(opt.stepSz_ms <= opt.binSize_ms, ...
    'NGL:invalidOption', ...
    'opt.stepSz_ms (%g) must be <= opt.binSize_ms (%g); otherwise bins do not overlap.', ...
    opt.stepSz_ms, opt.binSize_ms);

% Population-dynamics family (NGL02 / calculate_population_dynamics).
assert(isstruct(opt.popDyn), 'NGL:invalidOption', 'opt.popDyn must be a struct.');
for f = {'do','pca','jPCA','GPFA','trialEmbed','dropAborted'}
    assert(isfield(opt.popDyn, f{1}) && islogical(opt.popDyn.(f{1})) && isscalar(opt.popDyn.(f{1})), ...
        'NGL:invalidOption', 'opt.popDyn.%s must be a logical scalar.', f{1});
end
assert(isfield(opt.popDyn,'smoothSigma') && isnumeric(opt.popDyn.smoothSigma) && ...
       isscalar(opt.popDyn.smoothSigma) && opt.popDyn.smoothSigma >= 0, ...
    'NGL:invalidOption', 'opt.popDyn.smoothSigma must be a non-negative scalar (seconds).');
assert(isfield(opt.popDyn,'nComponents') && isnumeric(opt.popDyn.nComponents) && ...
       isscalar(opt.popDyn.nComponents) && opt.popDyn.nComponents >= 1 && ...
       opt.popDyn.nComponents == floor(opt.popDyn.nComponents), ...
    'NGL:invalidOption', 'opt.popDyn.nComponents must be a positive integer.');
assert(isfield(opt.popDyn,'alignIdx') && isnumeric(opt.popDyn.alignIdx) && ...
       isscalar(opt.popDyn.alignIdx) && opt.popDyn.alignIdx >= 1 && ...
       opt.popDyn.alignIdx == floor(opt.popDyn.alignIdx), ...
    'NGL:invalidOption', 'opt.popDyn.alignIdx must be a positive integer.');
assert(opt.popDyn.alignIdx <= numel(opt.alignto), ...
    'NGL:invalidOption', ...
    'opt.popDyn.alignIdx (%d) exceeds numel(opt.alignto) (%d).', ...
    opt.popDyn.alignIdx, numel(opt.alignto));
assert(isfield(opt.popDyn,'conditionVar') && ischar(opt.popDyn.conditionVar), ...
    'NGL:invalidOption', 'opt.popDyn.conditionVar must be a char (empty = no grouping).');
assert(isfield(opt.popDyn,'trialEmbedMethod') && ischar(opt.popDyn.trialEmbedMethod) && ...
       ismember(lower(opt.popDyn.trialEmbedMethod), {'pca','tsne','umap'}), ...
    'NGL:invalidOption', 'opt.popDyn.trialEmbedMethod must be ''PCA'', ''tSNE'', or ''UMAP''.');

% LFP artifact rejection (NGL02_LFP / artifact_detRej_lfp).
assert(isnumeric(opt.artZvalue) && isscalar(opt.artZvalue) && opt.artZvalue > 0, ...
    'NGL:invalidOption', 'opt.artZvalue must be a positive scalar (z-value cutoff).');
assert(ischar(opt.rejValue) || (isnumeric(opt.rejValue) && isscalar(opt.rejValue)), ...
    'NGL:invalidOption', ...
    'opt.rejValue must be a char (''zero''|''nan'') or a numeric scalar.');

% Project-specific gates.
for f = {'proj_chgDtctPCue','proj_socialLearning','proj_extintion','proj_FLIP'}
    assert(isfield(opt, f{1}) && islogical(opt.(f{1})) && isscalar(opt.(f{1})), ...
        'NGL:invalidOption', 'opt.%s must be a logical scalar.', f{1});
end

% Cross-subject FR PSTH plotter (NGL04_fireRate).
assert(isstruct(opt.fireRatePlot), 'NGL:invalidOption', 'opt.fireRatePlot must be a struct.');
fp = opt.fireRatePlot;
assert(isnumeric(fp.interval) && numel(fp.interval) == 2 && fp.interval(1) < fp.interval(2), ...
    'NGL:invalidOption', 'opt.fireRatePlot.interval must be [pre post] ms with pre<post.');
assert(isnumeric(fp.binSize_ms) && isscalar(fp.binSize_ms) && fp.binSize_ms > 0, ...
    'NGL:invalidOption', 'opt.fireRatePlot.binSize_ms must be a positive scalar (ms).');
assert(isnumeric(fp.stepSz_ms) && isscalar(fp.stepSz_ms) && fp.stepSz_ms > 0 && fp.stepSz_ms <= fp.binSize_ms, ...
    'NGL:invalidOption', 'opt.fireRatePlot.stepSz_ms must be a positive scalar (ms) <= binSize_ms.');
assert(islogical(fp.smoothPlot) && isscalar(fp.smoothPlot), ...
    'NGL:invalidOption', 'opt.fireRatePlot.smoothPlot must be a logical scalar.');
assert(isnumeric(fp.errAlpha) && isscalar(fp.errAlpha) && fp.errAlpha >= 0 && fp.errAlpha <= 1, ...
    'NGL:invalidOption', 'opt.fireRatePlot.errAlpha must be in [0,1].');
assert(iscell(fp.labelPriority) && ~isempty(fp.labelPriority) && all(cellfun(@ischar, fp.labelPriority)), ...
    'NGL:invalidOption', 'opt.fireRatePlot.labelPriority must be a non-empty cell of chars.');
assert(isnumeric(fp.busyWarnTraces) && isscalar(fp.busyWarnTraces) && fp.busyWarnTraces >= 1, ...
    'NGL:invalidOption', 'opt.fireRatePlot.busyWarnTraces must be a positive scalar.');
assert(ischar(fp.outDir), 'NGL:invalidOption', 'opt.fireRatePlot.outDir must be a char (empty for default).');

% Cross-session aggregation (NGL03_acrossSession).
assert(islogical(opt.aggregateSessions) && isscalar(opt.aggregateSessions), ...
    'NGL:invalidOption', 'opt.aggregateSessions must be a logical scalar.');
assert(islogical(opt.aggregateSubjects) && isscalar(opt.aggregateSubjects), ...
    'NGL:invalidOption', 'opt.aggregateSubjects must be a logical scalar.');
assert(~opt.aggregateSubjects || opt.aggregateSessions, ...
    'NGL:invalidOption', ...
    'opt.aggregateSubjects=true requires opt.aggregateSessions=true (subjects only aggregate after sessions).');

% Per-area context tag (NGL02 sets this per area before loadSpikes).
assert(ischar(opt.area) && ~isempty(opt.area), ...
    'NGL:invalidOption', ...
    'opt.area must be a non-empty char (defaults to ''all''; NGL02 overrides per area in multi-area mode).');

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
