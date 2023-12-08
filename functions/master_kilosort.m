function master_kilosort(sessions, input, varargin)
% Run Kilosort processing line programatically, without GUI. 
% Uses some info from the current session and searches for configuration
% and channel map files on '\analysisCode' folder.
%
% MAKE SURE your config file IS in that folder beforehand. Also, several
% chanMaps can co-exist, but then it will need to make explicit which one
% you want, modifying the opt.KSchanMapFile below.
%
% INPUT:    sessions: stores info about current session. Relevant here to
%                     know number of channels withlut hard coding it.
%                     *IMPORTANT: It is used by running the ConfigFile.
%           input:    stores general info about project, paths and so on. 
%           varargin: optional input (opt) that can be given or not.
%
% Winston's script and functions together with Sara's fixes.
%
% Version 15.06.2023 (Jesus)

if nargin < 3, opt = struct();
elseif nargin == 3, opt = varargin{1};
end

%% Defaults, if not given.
% Config and Channelmap files are to be found under '\analysisCode'
if ~isfield(opt,'KSConfigFile') || isempty(opt.KSConfigFile),           opt.KSConfigFile    = input.analysisCode;  end 
if ~isfield(opt,'KSchanMapFile') || isempty(opt.KSchanMapFile),         opt.KSchanMapFile   = ls(fullfile(input.analysisCode, 'chanMap*.mat')); end 
if ~isfield(opt,'spkTh') || isempty(opt.spkTh),                         opt.spkTh           = -4; end 

%% Find .bin files (raw and temp) % JESUS, changed the name and left only one.
% I assume it will be always in a SDD for processing.
rootfolder = opt.FolderProcDataMat; % the raw data binary file is in this folder (for current subject and session)

%% Set configuration. Will run 'kilosortConfig.m'
% Added all ops INSIDE config file.
ops = [];

% Non-existing config file in adequate folder.
if ~isfile(fullfile(opt.KSConfigFile, 'kilosortConfig.m'))
    warning('Config File not found under expected folder ''analysisCode''. Using a default version.')
    % If exists, use the standard one stored within the toolbox.
    if isfile(fullfile(input.toolbox, '\Instructions\kilosortConfig.m'))
        copyfile(fullfile(input.toolbox, '\Instructions\kilosortConfig.m'), input.analysisCode);
    else
        % It does not exist for some reason.
        error('Could not find the default configuration file for Kilosort. Skipped.')
    end
end

% Valid file found.
run(fullfile(opt.KSConfigFile, 'kilosortConfig.m'));

% Override the threshold if user opt are different from config file ops
if ops.spkTh ~= opt.spkTh
    ops.spkTh = opt.spkTh;
end

%% Check for Channel map file. Will run 'createChannelMapFile.m' if necessary.
if ~isfile(fullfile(opt.KSConfigFile, opt.KSchanMapFile))
    % No channel map file located in expected folder. Warn and create a basic linear one.
    warning('No channel map file found under expected folder ''\analysisCode''! Using a simple linear map.')
    run(fullfile(input.toolbox, '\functions\createChannelMapFile.m'));
end

%% Jesus. Included ops to test a check for chanMap-actual number of channels matching.
% It can happen that some channels are disabled or known dead. It will use the complete
% chanMap and find unmatching arrays.
% There is a logic variable within the map file named 'connected' which
% could be used to use (1) or not (0) that channel.
%
% ops.actual_channels = [sessions.info.INTAN_hdr.amplifier_channels.custom_order].';
% ops.actual_channels = ops.actual_channels + 1; % to match the 1-indexed map 
% [ops.chanMap, ~, ~, ~, ~] = loadChanMap(ops.Mapchan); % function to load channel map file
% if any(~ismember(ops.chanMap,ops.actual_channels))
%    ch = ops.chanMap(~ismember(ops.chanMap,ops.actual_channels));
%    ops.chanMap(~ismember(ops.chanMap,ops.actual_channels)) = [];
%    ops.chanMap(ops.chanMap > ch) = ops.chanMap(ops.chanMap > ch) - 1;
% end

%% This block runs all the steps of the algorithm
% 11.05 Jesus adding a way to resume after creation of .rez file, since the
% option is given.
fprintf('Looking for data inside %s \n', rootfolder)

if ~isfile(fullfile(rootfolder, 'rez.mat'))
    % Find the binary file
    fs          = dir(fullfile(rootfolder, '*.bin'));
    ops.fbinary = fullfile(rootfolder, fs(1).name);
    
    % Preprocess data to create temp_wh.dat
    rez = preprocessDataSub(ops);
    
    % Time-reordering as a function of drift
    rez = clusterSingleBatches(rez);

    % Saving here is a good idea, because the rest can be resumed after loading rez
    save(fullfile(rootfolder, 'rez.mat'), 'rez', '-v7.3');
else
    load(fullfile(rootfolder, 'rez.mat'), 'rez');
end

% Main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% Final merges
rez = find_merges(rez, 1);

% Final splits by SVD
rez = splitAllClusters(rez, 1);

% Final splits by amplitudes
rez = splitAllClusters(rez, 0);

% Decide on cutoff
rez = set_cutoff(rez);

fprintf('Found %d good units \n', sum(rez.good>0))

% Write to Phy
fprintf('Saving results to Phy \n')
rez2Phy(rez, rootfolder); % function has been modified to additionally output template_bestchannels.mat (Winston)

%% If you want to save the results to a Matlab file...
% % TODO. We can discuss if we want to go this way. We also need to find out
% exctly which file we need for extracting the data we want afterwards, and
% keep only those.
%
% % discard features in final rez file (too slow to save)
% rez.cProj = [];
% rez.cProjPC = [];
% 
% % save final results as rez2
% fprintf('Saving final results in rez2 \n')
% fname = fullfile(rootZ, 'rez2.mat');
% save(fname, 'rez', '-v7.3');

end