function param = bombcellConfig(param, path)
% JF, Load a parameter structure defining extraction and
% classification parameters
% 
% Inputs
% ephysMetaDir: dir() structure of the path to your .meta or .oebin meta
%   file
% rawFile: character array defining the path where your uncompressed raw
%   ephys data is
% 
% Outputs
% param: matlab structure defining extraction and
% classification parameters (see bc_qualityParamValues for required fields
% and suggested starting values)
% 
% Modified by Jesus 29/11/2023

% pre-Existing metrics?
param.qMetricsExist = ~isempty(dir(fullfile(path.savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(path.savePath, 'templates._bc_qMetrics.parquet')));

%% calculating quality metrics parameters 
param.plotDetails = 0; % lot of plots to check, debug or for a presentation

% plotting parameters 
param.plotGlobal    = 1; % plot summary of quality metrics 
param.verbose       = 1; % update user on progress
param.reextractRaw  = 1; % re extract raw waveforms or not 
param.ephysMetaFile = 'NaN';

% saving parameters 
param.saveAsTSV     = 1; % additionally save outputs in .tsv file - this is 
    % useful if you want to use phy after bombcell: each quality metric value
    % will appear as a column in the Cluster view
param.unitType_for_phy = 1; % whether to save the output of unitType in .tsv file for phy
param.saveMatFileForGUI = 1; % save certain outputs at .mat file - useful for GUI

% duplicate spikes parameters 
param.removeDuplicateSpikes = 1;
param.duplicateSpikeWindow_s = 0.00001; % in seconds 
param.saveSpikes_withoutDuplicates = 1;
param.recomputeDuplicateSpikes = 0;

% amplitude / raw waveform parameters
param.detrendWaveform = 1; % If this is set to 1, each raw extracted spike is
    % detrended (we remove the best straight-fit line from the spike)
    % using MATLAB's builtin function detrend.
param.saveMultipleRaw = 0; % If you wish to save the nRawSpikesToExtract 
param.decompressData = 0; % whether to decompress .cbin ephys data 
param.spikeWidth = 64; % width in samples. WAS 82
param.extractRaw = 1; % whether to extract raw waveforms or not 
param.probeType = []; % if you are using spikeGLX and your meta file does 
    % not contain information about your probe type for some reason
    % specify it here: '1' for 1.0 (3Bs) and '2' for 2.0 (single or 4-shanks)
    % For additional probe types, make a pull request with more
    % information.  If your spikeGLX meta file contains information about your probe
    % type, or if you are using open ephys, this paramater wil be ignored.

% signal to noise ratio
param.waveformBaselineNoiseWindow = 20; % time in samples at beginning of times
    % extracted to computer the mean raw waveform - this needs to be before the
    % waveform starts 

% refractory period parameters
param.tauR_valuesMin = 0.002; % refractory period time (s), usually 0.0020. 
    % If this value is different than param.tauR_valuesMax, bombcell will
    % estimate the tauR value taking possible values between :
    % param.tauR_valuesMin:param.tauR_valuesStep:param.tauR_valuesMax
param.tauR_valuesStep = 0.0005; %0.5/1000; % refractory period time (s) steps. Only 
    % used if param.tauR_valuesMin is different from param.tauR_valuesMax
param.tauR_valuesMax = 0.002; % refractory period time (s), usually 0.0020
param.tauC = 0.0001; % censored period time (s)

% percentage spikes missing parameters 
param.computeTimeChunks = 1; % compute fraction refractory period violations 
    % and percent spikes missing for different time chunks 
param.deltaTimeChunk = 360; %time in seconds 

% presence ratio 
param.presenceRatioBinSize = 60; % in seconds 

% drift estimate
param.driftBinSize = 60; % in seconds
param.computeDrift = 1; % whether to compute each units drift. this is a 
    % critically slow step that takes around 2seconds per unit 

% waveform parameters
param.waveformBaselineWindowStart = 20;
param.waveformBaselineWindowStop = 30; % in samples 
param.minThreshDetectPeaksTroughs = 0.2; % this is multiplied by the max value 
    % in a units waveform to give the minimum prominence to detect peaks using
    % matlab's findpeaks function.

% recording parameters
param.nChannels = 32; %number of recorded channels recorded in the raw data.
param.nSyncChannels = 1;

% distance metric parameters
param.computeDistanceMetrics = 1; % whether to compute distance metrics - this can be time consuming 
param.nChannelsIsoDist = 4; % number of nearby channels to use in distance metric computation 

%% classifying units into good/mua/noise parameters 
% whether to classify non-somatic units 
param.splitGoodAndMua_NonSomatic = 1;

% waveform 
param.maxNPeaks = 2; % maximum number of peaks
param.maxNTroughs = 1; % maximum number of troughs
param.somatic = 1; % keep only somatic units, and reject non-somatic ones
param.minWvDuration = 100; % in us
param.maxWvDuration = 1500; % in us
param.minSpatialDecaySlope = -0.003; % in a.u./um
param.maxWvBaselineFraction = 0.3; % maximum absolute value in waveform baseline
    % should not exceed this fraction of the waveform's abolute peak value

% distance metrics
param.isoDmin = 20; % minimum isolation distance value
param.lratioMax = 0.1; % maximum l-ratio value
param.ssMin = NaN; % minimum silhouette score 

% other classification params
param.minAmplitude = 20; % in uV
param.maxRPVviolations = 0.1; % fraction
param.RPV_tauR_estimate = NaN; 
param.maxPercSpikesMissing = 20; % in percentage
param.minNumSpikes = 1000; % number of spikes
param.maxDrift = 100;
param.minPresenceRatio = 0.7;
param.minSNR = 0.1;

end