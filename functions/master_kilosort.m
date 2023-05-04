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
%           input:    stores general info about project, paths and so on. 
%           varargin: optional input (opt) that can be given or not.
%
% Winston's script and functions and Sara's fixes.
%
% Version 27.04.2023 (Jesus)

%% Defaults, if not given as opt
if nargin < 3, opt = struct();
elseif nargin == 3, opt = varargin{1};
end

% Config and Channelmap files are defaulted to be found under '\analysisCode'
addpath(genpath('C:\KiloSort2_SpikeSorting')) % path to kilosort toolbox (Assumes Sorting PC, not local)
if ~isfield(opt,'KSConfigFile') || isempty(opt.KSConfigFile),       opt.KSConfigFile   = input.analysisCode;  end 
if ~isfield(opt,'KSchanMapFile') || isempty(opt.KSchanMapFile),     opt.KSchanMapFile  = ls(fullfile(input.analysisCode, 'chanMap*.mat')); end 

% Find .bin files (raw and temp) % JESUS, changed the name and left only one.
% I assume it will be always in a SDD for processing.
rootfolder = opt.FolderProcDataMat; % the raw data binary file is in this folder (for current subject and session)
% rootfolder = opt.FolderProcDataMat; % path to temporary binary file (same size as data, should be on fast SSD)

% Total time and channels to process
ops.trange      = [0 Inf]; % time range to sort (defaulted to the whole recording)
ops.NchanTOT    = sessions.info.nchannels; % total number of channels in your recording

% Set configuration, SSD and channel map
run(fullfile(opt.KSConfigFile, 'kilosortConfig.m'))
ops.fproc   = fullfile(rootfolder, 'temp_wh.dat'); % proc file on a fast SSD
ops.chanMap = fullfile(opt.KSConfigFile, opt.KSchanMapFile); %changed to find path 

%% This block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootfolder)

% % is there a channel map file in this folder? 
% JESUS removed this. we stick to chanmap under \analysisCode
% fs = dir(fullfile(rootfolder, 'chan*.mat'));
% if ~isempty(fs)
%     ops.chanMap = fullfile(rootfolder, fs(1).name);
% end

% find the binary file
fs          = dir(fullfile(rootfolder, '*.bin')); % JESUS, dir(.binfile) should work
ops.fbinary = fullfile(rootfolder, fs(1).name);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

% saving here is a good idea, because the rest can be resumed after loading rez
save(fullfile(rootfolder, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
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
% fprintf('Saving final results in rez2  \n')
% fname = fullfile(rootZ, 'rez2.mat');
% save(fname, 'rez', '-v7.3');

end