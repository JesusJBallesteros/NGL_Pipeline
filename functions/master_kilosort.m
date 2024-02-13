function master_kilosort(input, varargin)
% Run Kilosort processing line programatically, without GUI. 
% Uses some info from the current session and searches for configuration
% and channel map files on '\analysisCode' folder.
%
% MAKE SURE your config file IS in that folder beforehand. Also, several
% chanMaps can co-exist, but then it will need to make explicit which one
% you want, modifying the opt.KSchanMapFile below.
%
% INPUT:    input:    stores general info about project, paths and so on. 
%           varargin: optional input (opt) that can be given or not.
%
% Winston's script and functions together with Sara's fixes.
%
% Jesus. This function has been modified to allow for more than one run, with
% different spike thresholds, outputting more than one sortings
%
% Version 01.02.2024 (Jesus)

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Defaults, if not given.
% Config and Channelmap files are to be found under '\analysisCode'
if ~isfield(opt,'KSConfigFile') || isempty(opt.KSConfigFile),           opt.KSConfigFile    = input.analysisCode;  end 
if ~isfield(opt,'KSchanMapFile') || isempty(opt.KSchanMapFile),         opt.KSchanMapFile   = ls(fullfile(input.analysisCode, 'chanMap*.mat')); end 
if ~isfield(opt,'spkTh') || isempty(opt.spkTh),                         opt.spkTh           = -4; end 

% Raw data folder
rootfolder = opt.FolderProcDataMat; % the raw data binary file is in this folder (for current subject and session)    

%% Several runs requested (more than 1 values of threshold given)
% Prepare append for output folder
if length(opt.spkTh) > 1
    for i=1:length(opt.spkTh)
        repeattxt{i} = sprintf('%2.1f',opt.spkTh(i));
    end
end

%% Find .bin files (raw and temp).
% I assume it will be always in a SDD for processing.
if length(opt.spkTh) > 1
    for i=1:length(opt.spkTh)
        outfolder{i} = [rootfolder, '\' , repeattxt{i}]; % Subfolder for each threshold
        mkdir(outfolder{i})
    end
else
    outfolder{1} = rootfolder; % In this case, output goes directly at raw level
end

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

%% To override default config must be done after the config file is ran
for i=1:length(opt.spkTh)
        % Create temporal copies of ops
        ops_r{i} = ops;
        % The threshold, if user opt are different from defaults
        ops_r{i}.spkTh = opt.spkTh(i);
        % temp_wh folder (output folder)
        ops_r{i}.fproc = fullfile(outfolder{i}, 'temp_wh.dat'); % proc file on a fast SSD   
end

% Substitute the single original for the new cell struct
ops = ops_r;
clear ops_r % remove the temp

%% Check for Channel map file. Will run 'createChannelMapFile.m' if necessary.
if ~isfile(fullfile(opt.KSConfigFile, opt.KSchanMapFile))
    % No channel map file located in expected folder. Warn and create a basic linear one.
    warning('No channel map file found under expected folder ''\analysisCode''! Using a simple linear map.')
    run(fullfile(input.toolbox, '\functions\createChannelMapFile.m'));
end

%% Jesus. Included ops to test a check for chanMap-actual number of channels matching.
% It can happen that some channels are disabled, known dead, or just not being analyzed. 
% It will use the complete chanMap and find unmatching arrays.
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
for i=1:length(opt.spkTh)
    fprintf('Looking for data inside %s \n', rootfolder)
    fprintf('Kilosort outputs to %s \n', outfolder{i})

    if ~isfile(fullfile(outfolder{i}, 'rez.mat'))
        % Find the binary file
        fs          = dir(fullfile(rootfolder, '*.bin'));
        ops{i}.fbinary = fullfile(rootfolder, fs(1).name);
        
        % Preprocess data to create temp_wh.dat
        rez = preprocessDataSub(ops{i});
        
        % Time-reordering as a function of drift
        rez = clusterSingleBatches(rez);
    
        % Saving here is a good idea, because the rest can be resumed after loading rez
        save(fullfile(outfolder{i}, 'rez.mat'), 'rez', '-v7.3');
    else
        load(fullfile(outfolder{i}, 'rez.mat'), 'rez');
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
    fprintf('Saving results for Phy. \n')
    rez2Phy(rez, outfolder{i}); % function has been modified to additionally output template_bestchannels.mat (Winston)

    clear rez
end

end