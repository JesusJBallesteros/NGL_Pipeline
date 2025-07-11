% Creates parameter ans paths structures, defining extraction and classification parameters
%  
% Modified by Jesus 05/04/2024

%% Paths and other
% Find .bin files.
path.ephysKilosortPath  = [opt.FolderProcDataMat, '\kilosort2']; % the raw data binary file is in this folder (for current subject and session)
path.ephysRawDir        = dir([opt.FolderProcDataMat, '\*.*bin']); % your raw .bin data
path.savePath           = [opt.FolderProcDataMat, '\kilosort2']; % where you want to save the quality metrics

% Detect whether data is compressed. Decompress locally, if necessary.
path.rawFile = [path.ephysRawDir.folder, filesep, path.ephysRawDir.name]; % Ours is never .cbin, so far.

% Pre-Existing metrics?
param.qMetricsExist = ~isempty(dir(fullfile(path.savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(path.savePath, 'templates._bc_qMetrics.parquet')));

%% Switches
    param.rerun         = false;
    param.verbose       = false; % update user on progress
    param.plotDetails   = false; % lot of plots to check, debug or for a presentation
    param.reextractRaw  = true; % re-extract raw waveforms or not 

    % plotting parameters
    param.plotGlobal    = true; % plot summary of quality metrics 

    % saving parameters 
    param.saveAsTSV         = true; % save outputs also as .tsv file. To use phy after bombcell
    param.unitType_for_phy  = true; % to save the output of unitType in .tsv file for phy
    param.saveMatFileForGUI = true; % alos save outputs as .mat file. For GUI

    % duplicate spikes parameters 
    param.removeDuplicateSpikes         = true;
    param.saveSpikes_withoutDuplicates  = true;
    param.recomputeDuplicateSpikes      = false;

    % amplitude / raw waveform parameters
    param.detrendWaveform   = true; % If true, each raw extracted spike is detrended
    param.saveMultipleRaw   = true; % If you wish to save the nRawSpikesToExtract 
    param.decompressData    = false; % whether to decompress .cbin ephys data 
    param.extractRaw        = true; % whether to extract raw waveforms or not 
    
    param.computeDrift      = false; % whether to compute each units drift. this is critically slow step that takes around 2seconds per unit 
    param.computeTimeChunks = false; % compute fraction refractory period violations and percent spikes missing for different time chunks 
    param.somatic           = false; % keep only somatic units, and reject non-somatic ones
    param.computeDistanceMetrics = false; % whether to compute distance metrics - this can be time consuming 

%% Values
    % recording parameters
    param.nChannels         = 32; %number of channels recorded in the raw data.
    param.nSyncChannels     = 0;
    param.ephys_sample_rate = 32000;
    param.gain_to_uV        = 0.195;
    param.nRawSpikesToExtract = 5000;
    
    % duplicate spikes parameters 
    param.duplicateSpikeWindow_s = 0.00001; % in seconds
    
    % plotting parameters 
    param.ephysMetaFile = 'NaN';

    % amplitude / raw waveform parameters
    param.spikeWidth = 82; % width in samples.
    param.probeType = []; % For additional probe types. Not valid yet.

    % signal to noise ratio
    param.waveformBaselineNoiseWindow = 10; % samples at beginning extracted to compute the mean raw waveform. Needs to be before the waveform starts 

    % refractory period parameters
    param.tauR_valuesMin    = 0.002; % refractory period time (s)
    param.tauR_valuesStep   = 0.0005; % refractory period time (s) steps.
    param.tauR_valuesMax    = 0.002; % refractory period time (s)
    param.tauC              = 0.0005; % censored period time (s)

    % percentage spikes missing parameters 
    param.deltaTimeChunk    = 300; % time in seconds 

    % presence ratio 
    param.presenceRatioBinSize = 5; % in seconds 

    % drift estimate
    param.driftBinSize = 120; % in seconds

    % waveform parameters
    param.waveformBaselineWindowStart   = 21; % in samples 
    param.waveformBaselineWindowStop    = 30; % in samples 
    param.minThreshDetectPeaksTroughs   = 0.2; % this is multiplied by the max value 
        % in a units waveform to give the minimum prominence to detect peaks using
        % matlab's findpeaks function.

    % distance metric parameters
    param.nChannelsIsoDist = 4; % number of nearby channels to use in distance metric computation 

%% Classifying units parameters 
    % whether to classify non-somatic units 
    param.splitGoodAndMua_NonSomatic = true;

    % waveform 
    param.maxNPeaks     = 2; % maximum number of peaks
    param.maxNTroughs   = 1; % maximum number of troughs
    param.minWvDuration = 100; % in us
    param.maxWvDuration = 800; % in us
    param.minSpatialDecaySlope  = 0; % in a.u./um (was -0.003)
    param.maxWvBaselineFraction = 0.5; % maximum absolute value in waveform baseline
        % should not exceed this fraction of the waveform's abolute peak value (was 0.3)

    % distance metrics
    param.minIsoD   = 4; % minimum isolation distance value
    param.lratioMax = 0.1; % maximum l-ratio value
    param.ssMax     = NaN; % minimum silhouette score 

    % other classification params
    param.minAmplitude      = 30; % in uV
    param.maxRPVviolations  = 0.2; % fraction
    param.RPV_tauR_estimate = NaN;
    param.maxPercSpikesMissing = 40; % in percentage
    param.minNumSpikes      = 500; % number of spikes
    param.maxDrift          = 100;
    param.minPresenceRatio  = 0.7;
    param.minSNR            = 5;
