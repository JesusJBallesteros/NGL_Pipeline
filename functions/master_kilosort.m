function master_kilosort(sessions, varargin)
% Run Kilosort processing line programatically, without GUI. 
% Uses some info from the current session and searches for configuration
% and channel map files on '\analysisCode' folder.
%
% MAKE SURE your config file IS in that folder beforehand. Also, several
% chanMaps can co-exist, but then it will need to make explicit which one
% you want, modifying the opt.KSchanMapFile below.
%
% Based on Winston's script and functions.
%
% Version 27.04.2023 (Jesus)

addpath(genpath('C:\KiloSort2_SpikeSorting')) % path to kilosort folder and all its subfolder (Assumes Sorting PC, not local)

%% Defaults, if not given as opt
if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

% Config and Channelmap files are defaulted to be found under '\analysisCode'
if ~isfield(opt,'KSConfigFile') || isempty(opt.KSConfigFile),       opt.KSConfigFile   = input.analysisCode;  end
if ~isfield(opt,'KSchanMapFile') || isempty(opt.KSchanMapFile),     opt.KSchanMapFile  = ls('chanMap*.mat'); end

% Find .bin files (raw and temp)
% TODO  Why are both the same? can we get rid of one?
rootZ = opt.FolderProcDataMat; % the raw data binary file is in this folder (for current subject and session)
rootH = opt.FolderProcDataMat; % path to temporary binary file (same size as data, should be on fast SSD)

% Total time and channels to process
ops.trange      = [0 Inf]; % time range to sort (defaulted to the whole recording)
ops.NchanTOT    = sessions.info.nchannels; % total number of channels in your recording

% Set configuration, SSD and channel map
run(fullfile(opt.KSConfigFile, 'kilosortConfig.m'))
ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
ops.chanMap = fullfile(opt.KSConfigFile, opt.KSchanMapFile);

%% This block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
ops.fbinary = fullfile(rootZ, fs(1).name);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

% saving here is a good idea, because the rest can be resumed after loading rez
save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

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
rezToPhy(rez, rootZ); % function has been modified to additionally output template_bestchannels.mat (Winston)

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