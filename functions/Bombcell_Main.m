function Bombcell_Main(varargin)
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
% Jesus 21.12.2023

if nargin < 1, opt = struct();
elseif nargin == 1, opt = varargin{1};
end

%% Config
if ~isfield(opt,'rerun') || isempty(opt.rerun),                                     opt.rerun = 1;  end 
if ~isfield(opt,'nRawSpikesToExtract') || isempty(opt.nRawSpikesToExtract),         opt.nRawSpikesToExtract = 1000; end 
if ~isfield(opt,'ephys_sample_rate') || isempty(opt.ephys_sample_rate),             opt.ephys_sample_rate = 32000; end 
if ~isfield(opt,'gain_to_uV') || isempty(opt.gain_to_uV),                           opt.gain_to_uV = 0.195; end 

%% Defaults, if not given.
param = struct; % initialize bombcell param structure. Get opts 
    param.rerun = opt.rerun;
    param.nRawSpikesToExtract = opt.nRawSpikesToExtract; % how many raw spikes to extract for each unit 
    param.ephys_sample_rate = opt.ephys_sample_rate; % samples per second. 32KHz Deuteron, 30KHz Intan
    param.gain_to_uV = opt.gain_to_uV; % Same for Deuteron and Intan. (vs their openephys stuff)

%% Faster compute. Compile .mex file only if not done yet
if ~isfile('C:\Code\ephys-data-pipeline\toolboxes\bombcell\ephysProperties\helpers\CCGHeart.mexw64')
    orig = pwd;
    cd('C:\Code\ephys-data-pipeline\toolboxes\bombcell\ephysProperties\helpers');
    mex -O CCGHeart.c 
    cd(orig); clear orig
end

%% Set paths - EDIT THESE
% Find .bin files. I assume it will be always in a SDD for processing.
path.ephysKilosortPath  = opt.FolderProcDataMat; % the raw data binary file is in this folder (for current subject and session)
path.ephysRawDir        = dir([opt.FolderProcDataMat, '\*.*bin']); % your raw .bin data
path.savePath           = opt.FolderProcDataMat; % where you want to save the quality metrics

% We don't need this. Deprecating
path.decompressDataLocal = fullfile(path.ephysKilosortPath, 'decompressedData'); % where to save raw decompressed ephys data 
path.ephysMetaDir       = ''; % path to your meta file

% Detect whether data is compressed. Decompress locally, if necessary.
path.rawFile = [path.ephysRawDir.folder, filesep, path.ephysRawDir.name]; % Ours is never .cbin, so far.

% Set the rest of quality metric parameters, config file.
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
