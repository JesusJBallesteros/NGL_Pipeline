function opt = reduceChanMap(input,opt)
% reduceChanMap  Generate a reduced Kilosort channel map for partial recordings.
%
% PURPOSE:
%   Called from prepforsession when the number of amp*.dat files in the
%   session folder is fewer than opt.numChannels. This happens when some
%   electrode channels were disabled before the recording. Reads settings.xml
%   to identify which channels were enabled, removes disabled channels from
%   the original channel map, and saves the reduced map to both analysisCode/
%   and the Kilosort preprocessing folder.
%
% USAGE:
%   opt = reduceChanMap(input, opt)
%   opt.KSchanMapFile must point to the original (non-reduced) map file.
%
% INPUTS:
%   input  - struct with:
%              .analysisCode  path to analysisCode folder (map file lives here)
%   opt    - struct with:
%              .PathRaw       session raw-data folder (settings.xml is here)
%              .KSchanMapFile original channel map filename (.mat)
%              .KSfolder      Kilosort output folder (receives a copy)
%              .SavFileName   session name (tagged into the reduced map)
%
% OUTPUT:
%   opt    - updated:
%              .KSchanMapFile  new filename: '<original>_reduced.mat'
%
% SAVED FILES:
%   <analysisCode>/<map>_reduced.mat  — for Kilosort use on this session
%   <KSfolder>/<map>_reduced.mat      — for post-hoc analysis reference
%   Both files contain the original map fields with inactive channels removed,
%   plus 'chanID' (original probe channel IDs, 0- or 1-indexed as per source
%   map) and 'session' (recording session name for traceability).
%
% NOTES:
%   - chanMap is reindexed as a sort-order (1..nActive) after removing
%     inactive channels; Kilosort requires contiguous 1-based channel indices.
%   - Safe to re-run: if opt.KSchanMapFile already ends in '_reduced', the
%     '_reduced' suffix is stripped before loading the original.
%
% Version 09.10.2025 (Winston)

%% Load data
settings = xmlwrite(xmlread(fullfile(opt.PathRaw,'settings.xml')));          % Load text from settings.xml in raw file folder, which specifies which channels were enabled for the recording
status = regexp(settings,'"[A-Z]-\d\d\d" Enabled=("True"|"False")','match'); % Finds channels by INTAN naming convention, with their enabled value
idx = cellfun(@(x) strcmp(x,'True'),regexp(status,'(True|False)','match'));  % Logical index of present channels
mapname = erase(opt.KSchanMapFile,'_reduced');                               % Prevent recursion. If re-run, opt.KSchanMapFile may already point to the reduced map
map = load(fullfile(input.analysisCode, mapname));                           % load original channel map (not reduced)
%% Modify channel map by removing inactive channels
fnames = fieldnames(map);
for i = 1:length(fnames)
    if length(map.(fnames{i})) == length(status)                             % field matches total number of channels of probe
        map.(fnames{i}) = map.(fnames{i})(idx);                              % removes inactive channels by idx from struct field
    end
end
%% Kilosort requires chanMap as indices, not channel ID (chanMap(n) cannot be > nChannels)
map.chanID = map.chanMap;                                                    % Saving probe channel IDs, for future use when working with multiple recordings (Maintains 0-index or 1-index, be certain of your own format)
[~,map.chanMap] = sort(map.chanMap);                                         % Sort order is equivalent to channel order after excluded channels
map.session = opt.SavFileName;                                               % Attach session name to the reduced channel map for future reference if misplaced/moved
%% Saving new channel map
opt.KSchanMapFile = [erase(mapname,'.mat') '_reduced.mat'];                  % assign new name (X_reduced.mat)
save(fullfile(input.analysisCode, opt.KSchanMapFile),"-struct","map")        % save to analysisCode (for kilosort use)
if ~isfolder(opt.KSfolder)                                                   % sequentially, KS output folder is only created after kilosort (but required to save the reduced chanMap appropriately)
    mkdir(opt.KSfolder)
end
save(fullfile(opt.KSfolder, opt.KSchanMapFile),"-struct","map")              % save to preprocessing (for future analysis use / tagged to recording)
disp(sprintf('Channel map "%s" created in analysisCode and preprocessing folders, and will be used by Kilosort',opt.KSchanMapFile))
end