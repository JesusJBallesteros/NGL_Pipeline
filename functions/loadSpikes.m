function [spike] = loadSpikes(opt)
% Function that allows to extract data from phzton based files, obtained
% after KS-Phy. Will extract all remaining clusters, excluding the ones
% labeled as 'noise', if not explicitly asked for. By default, it will not
% load the PCs and it will extract a maximum of 2000 waveforms from the raw data.

%% Default parameters
if ~isfield(opt,'spparams'),    opt.spparams = struct('excludeNoise', true, 'loadPCs', false);  end % For cluster loading
if ~isfield(opt,'getwF'),       opt.getwF = false;                                              end % For waveform loading
if ~isfield(opt,'isibins'),     opt.isibins = 0:0.5:200;                                        end % For ISI binning, msec

% Waveform extraction option defaults
gwfparams = struct('dataType', 'int16',  ... % Data type of .dat file
                   'wfWin',    [-32 63], ... % Number of samples around spiketime to include in waveform.
                   'nWf',      opt.gwfparams.nWf,     ... % N waveforms to extrac tper unit.
                   'dataDir',  fullfile(opt.KSfolder), ... % KiloSort/Phy output folder
                   'fileName', fullfile(opt.FolderProcDataMat, [opt.SavFileName, '.bin']), ... % .dat file containing the raw 
                   'nCh',      [], ... % we don't know yet
                   'spikeTimes', [], 'spikeClusters', []); % Reserved for each cluster

%% Extract data from python files into a matlab friendly matrix
% if ~exist(fullfile(opt.spikeSorted, 'spike.mat'), "file")
    spikes = loadKSdir(opt.KSfolder, opt.spparams); % Helper function from Cortex-lab toolbox
    
    %% Get relevant info
    clusters        = sort(unique(spikes.cids)); % get and sort clusters by id
    nclust          = numel(clusters); % number of clusters
    
    % for waveforms
    gwfparams.nCh = spikes.n_channels_dat;

    % Allocate memory
    spike.label     = cell(1,nclust); % for clusters IDs
    spike.timestamp = cell(1,nclust); % spikes timestamps
    
    %% Proceed to extract timestamps for each cluster
    disp('Extracting curated clusters from Phy files. If many waveforms are requested, it may take a while.')
    for cl = 1:nclust
        spike.label{cl}         = num2str(clusters(cl));
        spike.timestamp{cl}     = spikes.st(spikes.clu==clusters(cl)); % in seconds
        spike.depth{cl}         = spikes.spikeDepths(spikes.clu==clusters(cl));
        spike.ampl{cl}          = spikes.spikeAmps(spikes.clu==clusters(cl));
        spike.templampl{cl}     = spikes.tempScalingAmps(spikes.clu==clusters(cl));

        %% Extract waveforms
        if opt.getwF
            % a few more params for 'getWaveForms' dep on cluster
            gwfparams.spikeTimes = ceil(spike.timestamp{cl}*spikes.sample_rate); % Vector of cluster spike times (in samples) same length as .spikeClusters
            gwfparams.spikeClusters = spikes.clu(spikes.clu==clusters(cl));
                       
            % Get waveforms
            wF = getWaveForms(gwfparams);
        
            % Refine
            wF.waveForms = squeeze(wF.waveForms);
            wF.waveFormsMean = squeeze(wF.waveFormsMean);
        
            % Find averaged max amplitude channel
            [wF.maxamplch, ~, ~]    = find(wF.waveFormsMean==min(min(wF.waveFormsMean)));
            spike.waveform{cl}      = permute(wF.waveForms(:,wF.maxamplch,:),[2,3,1]);
            spike.waveform{cl}      = squeeze(spike.waveform{cl});
            spike.waveFormsMean{cl} = wF.waveFormsMean;
         end
    
    end

% Calculate the ISI histogram for all clusters
[spike.isihist] = calc_isihist(spike, opt);

end