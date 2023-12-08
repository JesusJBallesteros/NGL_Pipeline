function [spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, ...
    pcFeatures, pcFeatureIdx, channelPositions, goodChannels] = bc_loadEphysData(path)
% JF, Load ephys data (1-indexed)
% ------
% Inputs
% ------
% ephys_path.ephysKilosortPath: character array defining the path.ephysKilosortPath to your kilosorted output files 
% ------
% Outputs
% ------
% spikeTimes_samples: nSpikes × 1 uint64 vector giving each spike time in samples (*not* seconds)
% spikeTemplates: nSpikes × 1 uint32 vector giving the identity of each
%   spike's matched template
% templateWaveforms: nTemplates × nTimePoints × nChannels single matrix of
%   template waveforms for each template and channel
% templateAmplitudes: nSpikes × 1 double vector of the amplitude scaling factor
%   that was applied to the template when extracting that spike
% pcFeatures: nSpikes × nFeaturesPerChannel × nPCFeatures  single
%   matrix giving the PC values for each spike
% pcFeatureIdx: nTemplates × nPCFeatures uint32  matrix specifying which
%   channels contribute to each entry in dim 3 of the pc_features matrix
% channelPositions: nChannels x 2 double matrix, each row gives the x and y 
%   coordinates of each channel
% goodChannels: nChannels x 1 uint32 vector defining the channels used by
%   kilosort (some are dropped during the spike sorting process)
%
% Modified by Jesus 29/11/23

% Read spike templates
spike_templates_0idx = readNPY([path.ephysKilosortPath filesep 'spike_templates.npy']);

% Index them for matlab
spikeTemplates = spike_templates_0idx + 1;

% Read Spike times
if exist(fullfile(path.ephysKilosortPath,'spike_times_corrected.npy'), "file") % When running pyKS stitched you need the 'aligned / corrected' spike times
    spikeTimes_samples = double(readNPY([path.ephysKilosortPath filesep  'spike_times_corrected.npy']));
    %spikeTimes_datasets = double(readNPY([ephys_path.ephysKilosortPath filesep  'spike_datasets.npy'])) + 1; %  which dataset? (zero-indexed so +1)
else
    spikeTimes_samples = double(readNPY([path.ephysKilosortPath filesep 'spike_times.npy']));
    %spikeTimes_datasets = ones(size(spikeTimes_samples));
end
templateAmplitudes = readNPY([path.ephysKilosortPath filesep 'amplitudes.npy']);

% Load and unwhiten templates
templateWaveforms_whitened = readNPY([path.ephysKilosortPath filesep 'templates.npy']);
winv = readNPY([path.ephysKilosortPath filesep 'whitening_mat_inv.npy']);
templateWaveforms = zeros(size(templateWaveforms_whitened));
for t = 1:size(templateWaveforms,1)
    templateWaveforms(t,:,:) = squeeze(templateWaveforms_whitened(t,:,:))*winv;
end

% Read principal component features
pcFeatures = readNPY([path.ephysKilosortPath filesep  'pc_features.npy']);
pcFeatureIdx = readNPY([path.ephysKilosortPath filesep  'pc_feature_ind.npy']) + 1;

% Get channel topography and valid count
channelPositions = readNPY([path.ephysKilosortPath filesep  'channel_positions.npy']); 
goodChannels = readNPY([path.ephysKilosortPath filesep  'channel_map.npy']) + 1;

end