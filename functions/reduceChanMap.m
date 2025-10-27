function opt = reduceChanMap(input,opt)
% This function is used when the number of raw data files do not match the
% intended/expected number of channels specified in NGL_SetAndRunMe.m, and
% a reduced channel map needs to be created for Kilosort to correctly
% assign data to channels (coordinates, etc.)
%
% OUTPUT:
% A reduced version of the original channel map excluding inactive channels
% specified by "settings.xml" in the raw data folder,
% "[Channel_Map_Name]_reduced.mat" is created (or overwritten) in the same
% analysisCode folder specified in input.analysisCode for Kilosort. 
% Another copy is saved to the preprocessing folder for the current 
% data/recording specified in opt.KSfolder, for future use in analysis and
% plotting. Channel IDs (not indices) corresponding to location on
% probe are saved in the new chanID variable
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
save(fullfile(opt.KSfolder, opt.KSchanMapFile),"-struct","map")              % save to preprocessing (for future analysis use / tagged to recording)
disp(sprintf('Channel map "%s" created in analysisCode and preprocessing folders, and will be used by Kilosort',opt.KSchanMapFile))
end