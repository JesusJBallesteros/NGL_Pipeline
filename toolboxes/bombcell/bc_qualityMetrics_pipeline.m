%% Adapted Bombcell pipeline (From JF-example pipeline)
% Set the paths here and the parameters in 'bc_qualityParamValues'
% This pipeline will:
%   (1) load your kilosorted data, 
%   (2) run bombcell on it, save the output and
%   (3) bring up summary plots.
% The first time, this pipeline will be significantly slower (10-20' more)
% than after because it extracts raw waveforms. Subsequent times these
% pre-extracted waveforms are simply loaded in.
% We recommend running this pipeline on a few datasets and deciding on
% quality metric thresholds depending on the summary plots (histograms 
% of the distributions of quality metrics for each unit) and GUI. 

%% Dependencies
addpath(genpath('C:\Code\bombcell'))
addpath(genpath('C:\Code\ephys-data-pipeline\toolboxes\npy-matlab'))

% Faster compute
cd('C:\Code\bombcell\ephysProperties\helpers');
mex -O CCGHeart.c 

%% Set paths - EDIT THESE
path.ephysKilosortPath  = 'F:\Pilot_SocialLearning\data\preprocessing\257\20231108\';% path to your kilosort output files 
path.ephysRawDir        = dir([path.ephysKilosortPath, '*.*bin']); % your raw .bin data
path.saveLocation       = path.ephysKilosortPath; % where you want to save the quality metrics 
path.savePath           = path.saveLocation;
path.decompressDataLocal = fullfile(path.ephysKilosortPath, 'decompressedData'); % where to save raw decompressed ephys data 
path.ephysMetaDir       = ''; % path to your meta file
% Detect whether data is compressed. Decompress locally, if necessary.
path.rawFile = [path.ephysRawDir.folder, filesep, path.ephysRawDir.name]; % Ours is never .cbin, so far.

%% Defaults
param = struct; % initialize structure 
     param.rerun        = 1; % To re-run and overwrite previous analisys
     param.nRawSpikesToExtract = 1000; % how many raw spikes to extract for each unit 
     
     % Recording system specific:
     param.ephys_sample_rate = 32000; % samples per second. 32KHz Deuteron, 30KHz Intan
     param.gain_to_uV   = 0.195; % Same for Deuteron and Intan. (vs their openephys stuff)

%Set the rest of quality metric parameters, config file:
param = bc_qualityParamValues(param, path);

%% Load data from Kilosort outputs.
[spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions] ...
     = bc_loadEphysData(path);

%% Compute quality metrics 
% First run
if ~param.qMetricsExist || param.rerun
    [qMetric, unitType] = bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
                                                  templateWaveforms, templateAmplitudes, pcFeatures, ...
                                                  pcFeatureIdx, channelPositions, path);
% If previous metrics exist
else
    [param, qMetric] = bc_loadSavedMetrics(path);
    unitType = bc_getQualityUnitType(param, qMetric, savePath);
end
