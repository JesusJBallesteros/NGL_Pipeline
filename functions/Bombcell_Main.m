function Bombcell_Main(input, opt)
% Adapted Bombcell pipeline 
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
%
% Jesus 29.02.2023

param = struct; % initialize bombcell param structure.
path = struct; % initialize bombcell path structure.

%% Non-existing config file in adequate folder.
if ~isfile(fullfile(input.analysisCode, 'bombcellConfig.m'))
    warning('Config File not found under expected folder ''analysisCode''. Using a default version.')
    % If exists, use the standard one stored within the toolbox.
    if isfile(fullfile(input.toolbox, '\Instructions\bombcellConfig.m'))
        copyfile(fullfile(input.toolbox, '\Instructions\bombcellConfig.m'), input.analysisCode);
    else
        % It does not exist for some reason.
        error('Could not find the default configuration file for Kilosort. Skipped.')
    end
end

% Valid file found. Run it.
run(fullfile(input.analysisCode, 'bombcellConfig.m'));

%% Override params based on opts
if isfield(opt,'rerun'),                param.rerun = opt.rerun;  end 
if isfield(opt,'nRawSpikesToExtract'),  param.nRawSpikesToExtract = opt.nRawSpikesToExtract; end 
if isfield(opt,'ephys_sample_rate'),    param.ephys_sample_rate = opt.ephys_sample_rate; end 
if isfield(opt,'gain_to_uV'),           param.gain_to_uV = opt.gain_to_uV; end 

% %% Set the rest of quality metric parameters, config function.
% [param, path] = bombcellConfig(param, opt);

%% Faster compute. Compile .mex file only if not done yet
if ~isfile('C:\Code\ephys-data-pipeline\toolboxes\bombcell\ephysProperties\helpers\CCGHeart.mexw64')
    orig = pwd;
    cd('C:\Code\ephys-data-pipeline\toolboxes\bombcell\ephysProperties\helpers');
    mex -O CCGHeart.c 
    cd(orig); clear orig
end

%% Load data from Kilosort outputs.
[spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions] ...
     = bc_loadEphysData(path);

%% Compute quality metrics 
% First run
if ~param.qMetricsExist || param.rerun
    [qMetric, unitType] = bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
                                                  templateWaveforms, templateAmplitudes, pcFeatures, ...
                                                  pcFeatureIdx, channelPositions, path);
else % If previous metrics exist
    [param, qMetric] = bc_loadSavedMetrics(path);
    unitType = bc_getQualityUnitType(param, qMetric, savePath);
end
