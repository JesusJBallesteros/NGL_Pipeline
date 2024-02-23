%% Configuration file for master_kilosort.m.
% Place this .m file inside '...\projectName\analysisCode' folder.
% Originally named 'StandardConfig_MOVEME.m' or 'configFile384.m'. 
% Original version can be found at ...\KiloSort2\Kilosort2-master\configFiles
%
% Version at 22.02.2024 by Jesus.

%% Often modified options
ops.spkTh       = -4;   % spike threshold in standard deviations (-6).
ops.ntbuff      = 32;   % samples of symmetrical buffer for whitening and spike detection
ops.NT          = 32*1000 + ops.ntbuff; %  must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
ops.CAR         = 1; % Common average referencing (median)
ops.AUCsplit    = 0.85; % splitting a cluster requires at least this much isolation for each sub-cluster (max = 1)
ops.minfr_goodchannels = 0.5; % minimum firing rate on a "good" channel (0 to skip)
ops.trange      = [0 Inf]; % Total time and channels to process (defaulted to the whole recording)

% ops.criterionNoiseChannels = 1; % EXPERIMENTAL. allows a fraction of all templates to span more channel groups.

%% Less often modified options
ops.lam          = 10; % how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.minFR        = 0.02; % minimum spike rate (Hz), if a cluster falls below this for too long it gets removed. was 1/50
ops.Th           = [10 4]; % threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.momentum     = [20 400]; % number of samples to average over (annealed from first to second value) 
ops.Nfilt        = 300;  % max number of clusters
ops.nfilt_factor = 4;    % max number of clusters per good channel (even temporary ones)

%% Options not to play around with. For determining PCs
ops.fshigh          = 500; % frequency for high pass filtering that Kilosort will do. We have pre-filtered anyways. TO REMOVE?
ops.sigmaMask       = 50; % spatial constant in um for computing residual variance of spike
ops.ThPre           = 4; % threshold crossings for pre-clustering (in PCA projection space)
ops.whiteningRange  = 32;   % number of channels to use for whitening each channel
ops.reorder         = 1;    % whether to reorder batches for drift correction. 
ops.nskip           = 2;    % how many batches to skip for determining spike PCs
ops.nSkipCov        = 25;    % compute whitening matrix from every N-th batch % was 25
ops.nPCs            = 3;    % how many PCs to project the spikes into
ops.scaleproc       = 200;  % int16 scaling of whitened data

%% Not really options. Necessary things.
ops.fproc   = fullfile(rootfolder, 'temp_wh.dat');  % proc file on a fast SSD. temp_wh folder
ops.chanMap = fullfile(opt.KSConfigFile, opt.KSchanMapFile); % channel map folder
ops.fs      = input.sessions(input.run(1)).info.amplifier_sample_rate; % sample rate can be obtain from session info.

try     ops.NchanTOT  = input.sessions(input.run(1)).info.numChannels;  % total number of channels in your recording
catch,  ops.NchanTOT  = input.sessions(input.run(1)).info.nChannels;  end % total number of channels in your recording

if isempty(ops.NchanTOT), ops.NchanTOT  = opt.numChannels; end % From Deuteron will be empty. Until we use one specific headstage, we override it with a general option.

ops.useRAM              = 0;    % not yet available
ops.GPU                 = 1;    % has to be 1, no CPU version yet, sorry
