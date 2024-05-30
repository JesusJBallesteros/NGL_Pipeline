function [spike] = loadSpikes(opt)
% Based on Juan's script for his own piloting

if ~exist(fullfile(opt.spikeSorted, [opt.SavFileName, '.mat']), "file")
    % Extract data from python files into a matlab friendly matrix
    spikes = loadKSdir(opt.KSfolder);
    
    % get relevant info
    clusters        = sort(unique(spikes.cids)); % get and sort clusters by id
    nclust          = numel(clusters); % number of clusters
    
    % Alocate memory
    spike.label     = cell(1,nclust); % for clusters IDs
    spike.timestamp = cell(1,nclust); % spikes timestamps
    
    % Proceed to extract timestamps for each cluster
    disp('Extracting curated clusters from Phy files. If many waveforms are requested, it may take a while.')
    for cl = 1:nclust
        spike.label{cl}           = num2str(clusters(cl));
        spike.timestamp{cl}       = spikes.st(spikes.clu==clusters(cl))*1000; % from ms to sec
        
        if opt.getwF
            % a few more params for 'getWaveForms' dep on cluster
            gwfparams.spikeTimes = ceil(extract.st(extract.clu==clusters(cl))*32000); % Vector of cluster spike times (in samples) same length as .spikeClusters
            gwfparams.spikeClusters = extract.clu(extract.clu==clusters(cl));
            
            % Get waveform
            wF = getWaveForms(gwfparams);
        
            % Refine
            wF.waveForms = squeeze(wF.waveForms);
            wF.waveFormsMean = squeeze(wF.waveFormsMean);
        
            % Find averaged max amplitude channel
            [wF.maxamplch, ~, ~] = find(wF.waveFormsMean==min(min(wF.waveFormsMean)));
            spike.waveform{cl}   = permute(wF.waveForms(:,wF.maxamplch,:),[2,3,1]);
         end
    
    end
    
    save(fullfile(opt.spikeSorted, [opt.SavFileName, '.mat']), "spike", '-mat');
else
    disp('Already existing Spike-sorted data for this session. Loading instead.')
    spike = load(fullfile(opt.spikeSorted, [opt.SavFileName, '.mat']));
end
end