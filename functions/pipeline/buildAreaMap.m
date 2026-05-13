function areaMap = buildAreaMap(Areas, chanMapPath)
% buildAreaMap  Map area labels to channel groups and build per-area KS masks.
%
% PURPOSE:
%   Called once during set_default when input.Areas is provided. Loads the
%   Kilosort chanMap .mat file, maps each area label to the kcoords group(s)
%   it represents, derives the per-channel index sets for each unique area,
%   and builds per-area 'connected' vectors. These vectors are saved as
%   individual .mat files alongside the source chanMap so master_kilosort4
%   can pass them directly to the Python/KS4 run for each area.
%
%   The connected mask logic is:
%     connected_area(i) = 1  if channel i has kcoords matching this area
%     connected_area(i) = 0  otherwise
%   This overrides the original connected field, enabling KS4 to sort only
%   the channels that physically belong to the current area.
%
% USAGE:
%   areaMap = buildAreaMap(input.Areas, fullfile(input.analysisCode, opt.KSchanMapFile))
%   Called from set_default (Section 4b) when input.Areas is defined.
%
% INPUTS:
%   Areas       - {1×N} cell array of area labels, one per kcoords group.
%                 Areas{i} labels the channels where kcoords == i.
%                 Repeated labels assign multiple shanks to the same area.
%                 Example: {'NCL','NCL','STR'} for a 3-group 64-channel probe.
%   chanMapPath - Full path to the Kilosort chanMap .mat file.
%                 Must contain at minimum the fields: kcoords, connected.
%
% OUTPUT:
%   areaMap - struct with fields:
%     .uniqueAreas       {1×nAreas} unique area labels (first-appearance order)
%     .kgroups_per_area  {1×nAreas} kcoords values that belong to each area
%     .chanIdx_per_area  {1×nAreas} 1-based channel index vectors per area
%     .nChans_per_area   [1×nAreas] number of channels per area
%     .connected_masks   {1×nAreas} logical connected vectors (full length)
%     .chanMapFiles      {1×nAreas} full paths to per-area chanMap .mat files
%     .chanMapPath       char, source chanMap file path
%
%   Per-area chanMap files are written as:
%     <chanMapDir>/<chanMapBasename>_<AreaLabel>.mat
%   Each file is a copy of the source chanMap with the 'connected' field
%   replaced by the per-area mask.
%
% NOTES:
%   - kcoords values in the chanMap must be contiguous integers 1..N, where N
%     equals numel(Areas). The user must ensure kcoords are laid out this way.
%   - chanMap must be saved in a location that is writable (analysisCode/).
%
% CALLS:
%   load, save (built-in MATLAB file I/O)
%
% Last modified 11.05.2026 (Jesus)

%% Input validation
assert(iscell(Areas) && ~isempty(Areas), ...
    'NGL:buildAreaMap', 'input.Areas must be a non-empty cell array of char labels.');
assert(all(cellfun(@(s) ischar(s) || isstring(s), Areas)), ...
    'NGL:buildAreaMap', 'All entries in input.Areas must be character arrays or strings.');
Areas = cellfun(@char, Areas, 'UniformOutput', false); % normalise to char

assert(isfile(chanMapPath), ...
    'NGL:buildAreaMap', 'chanMap file not found: %s', chanMapPath);

%% Load chanMap
m = load(chanMapPath);
assert(isfield(m, 'kcoords'), ...
    'NGL:buildAreaMap', 'chanMap file must contain a ''kcoords'' field.');
assert(isfield(m, 'connected'), ...
    'NGL:buildAreaMap', 'chanMap file must contain a ''connected'' field.');

kcoords = m.kcoords(:);            % ensure column vector [nChan × 1]
nChan   = numel(kcoords);

%% Validate Areas covers all kcoords groups
nGroups = numel(Areas);
kMax    = max(kcoords);
assert(nGroups >= kMax, ...
    'NGL:buildAreaMap', ...
    ['input.Areas has %d entries but chanMap kcoords go up to %d. ' ...
     'Provide one Areas label per kcoords group (1 through %d).'], ...
    nGroups, kMax, kMax);

%% Identify unique areas (preserving first-appearance order)
[uniqueAreas, firstIdx] = unique(Areas(1:kMax), 'stable');
nAreas = numel(uniqueAreas);

fprintf('buildAreaMap: %d kcoords groups -> %d unique areas (%s)\n', ...
    kMax, nAreas, strjoin(uniqueAreas, ', '));

%% Pre-allocate output
areaMap.uniqueAreas      = uniqueAreas;
areaMap.kgroups_per_area = cell(1, nAreas);
areaMap.chanIdx_per_area = cell(1, nAreas);
areaMap.nChans_per_area  = zeros(1, nAreas);
areaMap.connected_masks  = cell(1, nAreas);
areaMap.chanMapFiles     = cell(1, nAreas);
areaMap.chanMapPath      = chanMapPath;

%% Build per-area chanMap files
[chanMapDir, chanMapBase, ~] = fileparts(chanMapPath);

for a = 1:nAreas
    label   = uniqueAreas{a};

    % kcoords indices belonging to this area
    kgroups = find(strcmp(Areas(1:kMax), label));
    areaMap.kgroups_per_area{a} = kgroups;

    % Channel indices (1-based) where kcoords is in kgroups
    inArea  = ismember(kcoords, kgroups);
    chanIdx = find(inArea);
    areaMap.chanIdx_per_area{a} = chanIdx;
    areaMap.nChans_per_area(a)  = numel(chanIdx);

    % Per-area connected mask: 1 for this area's channels, 0 for all others
    conn = false(nChan, 1);
    conn(inArea) = true;
    areaMap.connected_masks{a} = conn;

    % Build the per-area chanMap struct (copy source, swap connected)
    mArea           = m;
    mArea.connected = conn;

    % Save to <chanMapDir>/<base>_<label>.mat
    outFile = fullfile(chanMapDir, sprintf('%s_%s.mat', chanMapBase, label));
    save(outFile, '-struct', 'mArea');
    areaMap.chanMapFiles{a} = outFile;

    fprintf('  Area %-6s: kcoords [%s], %d channels -> %s\n', ...
        label, num2str(kgroups), numel(chanIdx), outFile);
end
