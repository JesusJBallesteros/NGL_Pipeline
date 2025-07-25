function spikeStruct = loadKSdir(ksDir, varargin)

if ~isempty(varargin)
    params = varargin{1};
else
    params = struct();
end

if ~isfield(params, 'excludeNoise') || isempty(params.excludeNoise),  params.excludeNoise = true;   end
if ~isfield(params, 'loadPCs') || isempty(params.loadPCs),            params.loadPCs      = false;  end

% Load Parameters from the KS/Phy run
spikeStruct = loadParamsPy(fullfile(ksDir, 'params.py'));

% Load the spike times
ss = readNPY(fullfile(ksDir, 'spike_times.npy')); % Spike times in samples
st = double(ss)/spikeStruct.sample_rate;    % Spike times in sec

% Load the spike templates index
spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % Zero-indexed

% Load the spike cluster index
if exist(fullfile(ksDir, 'spike_clusters.npy'), "file")
    clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

% Load the spike amplitudes
tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));

% Load the principal components features
if params.loadPCs
    pcFeat = readNPY(fullfile(ksDir,'pc_features.npy'));    % nSpikes x nFeatures x nLocalChannels
    pcFeatInd = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat = [];
    pcFeatInd = [];
end

% Load the clusters tags
cgsFile = '';
if exist(fullfile(ksDir, 'cluster_groups.csv'), "file") 
    cgsFile = fullfile(ksDir, 'cluster_groups.csv');
elseif exist(fullfile(ksDir, 'cluster_group.tsv'), "file") 
   cgsFile = fullfile(ksDir, 'cluster_group.tsv');
end 

if ~isempty(cgsFile)
    [cids, cgs] = readClusterGroupsCSV(cgsFile);

    % Remove the clusters tagges as Noise
    if params.excludeNoise
        noiseClusters = cids(cgs==0);
        st = st(~ismember(clu, noiseClusters));
        spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
        tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));        
        
        if params.loadPCs
            pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
            %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
        end
        
        clu = clu(~ismember(clu, noiseClusters));
        cgs = cgs(~ismember(cids, noiseClusters));
        cids = cids(~ismember(cids, noiseClusters));
    end
    
else
    clu = spikeTemplates;
    cids = unique(spikeTemplates);
    cgs = 3*ones(size(cids));
end
    
% Load channel mapping
coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
ycoords = coords(:,2); xcoords = coords(:,1);

% Load templates
temps = readNPY(fullfile(ksDir, 'templates.npy'));

% Load whitening matrix 
winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

% Extract and process further info
[spikeAmps, spikeDepths, ~, ~, ~, templateDuration, ~] = ...
    templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps);

% Collect for output structure
spikeStruct.st = st;
spikeStruct.spikeTemplates = spikeTemplates;
spikeStruct.clu = clu;
spikeStruct.tempScalingAmps = tempScalingAmps;
spikeStruct.cgs = cgs;
spikeStruct.cids = cids;
spikeStruct.xcoords = xcoords;
spikeStruct.ycoords = ycoords;
spikeStruct.temps = temps;
spikeStruct.winv = winv;
spikeStruct.pcFeat = pcFeat;
spikeStruct.pcFeatInd = pcFeatInd;

spikeStruct.spikeAmps = spikeAmps;
spikeStruct.spikeDepths = spikeDepths;
spikeStruct.templateDuration = templateDuration;

end