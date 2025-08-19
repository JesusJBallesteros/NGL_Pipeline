% Creates parameter ans paths structures, defining extraction and classification parameters
%  
% These defaults were check by Jesus 19/08/2025

%% Paths and other
% Find .bin files.
path.ephysKilosortPath  = [opt.FolderProcDataMat, '\kilosort4']; % the raw data binary file is in this folder (for current subject and session)
path.ephysRawDir        = dir([opt.FolderProcDataMat, '\*.*bin']); % your raw .bin data
path.savePath           = [opt.FolderProcDataMat, '\kilosort4']; % where you want to save the quality metrics

% Detect whether data is compressed. Decompress locally, if necessary.
path.rawFile = [path.ephysRawDir.folder, filesep, path.ephysRawDir.name]; % Ours is never .cbin, so far.

% Pre-Existing metrics?
param.qMetricsExist = ~isempty(dir(fullfile(path.savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(path.savePath, 'templates._bc_qMetrics.parquet')));

%% Switches
    param.rerun         = true;
    param.verbose       = true; % update user on progress
    param.plotDetails   = false; % lot of plots to check, debug or for a presentation
    param.reextractRaw  = false; % re-extract raw waveforms or not 

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
    param.computeDistanceMetrics = true; % whether to compute distance metrics - this can be time consuming 

%% Values
    % recording parameters
    param.nChannels         = 64; %number of channels recorded in the raw data.
    param.nSyncChannels     = 0;
    param.ephys_sample_rate = input.sessions.info.amplifier_sample_rate;
    param.nRawSpikesToExtract = 5000;
    
    % duplicate spikes parameters 
    param.duplicateSpikeWindow_s = 0.00001; % in seconds
    
    % plotting parameters 
    param.ephysMetaFile = 'NaN';

    % amplitude / raw waveform parameters
    if param.ephys_sample_rate == 30000 % INTAN
        param.gain_to_uV = 0.195;
        param.spikeWidth = 61; 
    elseif param.ephys_sample_rate == 32000 % DEUTERON
        param.gain_to_uV = 1;
        param.spikeWidth = 65; 
    end
    param.probeType = []; % For additional probe types. Not valid yet.

    % refractory period parameters
    param.tauR_valuesMin    = 1.5/1000; % refractory period time (s)
    param.tauR_valuesStep   = 0.1/1000; % refractory period time (s)
    param.tauR_valuesMax    = 2.4/1000; % refractory period time (s)
    param.tauC              = 0.0005; % censored period time (s)

    % percentage spikes missing parameters 
    param.deltaTimeChunk    = 360; % time in seconds 

    % presence ratio 
    param.presenceRatioBinSize = 60; % in seconds 

    % drift estimate
    param.driftBinSize = 300; % in seconds

    % Now calculate all spike-width-dependent parameters
    % Signal to noise ratio - baseline noise window
    if param.spikeWidth <= 70  % Shorter waveforms (like KS4)
        param.waveformBaselineNoiseWindow = round(param.spikeWidth * 10/61); % Scale from standard 10 samples at 61 width
    else  % Longer waveforms
        param.waveformBaselineNoiseWindow = round(param.spikeWidth * 20/82); % Scale from standard 20 samples at 82 width
    end
    param.waveformBaselineNoiseWindow = max(5, param.waveformBaselineNoiseWindow); % Ensure at least 5 samples
    
    % Waveform baseline windows
    if param.spikeWidth <= 70  % Shorter waveforms (like KS4)
        param.waveformBaselineWindowStart = max(1, round(param.spikeWidth * 1/61));
        param.waveformBaselineWindowStop = max(5, round(param.spikeWidth * 11/61)); % in samples 
    else  % Longer waveforms
        param.waveformBaselineWindowStart = max(1, round(param.spikeWidth * 20/82));
        param.waveformBaselineWindowStop = max(10, round(param.spikeWidth * 30/82)); % in samples 
    end

    % waveform parameters
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
    param.maxWvDuration = 950; % in us
    param.minSpatialDecaySlope  = 0.005; % in a.u./um (was -0.003)
    param.maxWvBaselineFraction = 0.35; % maximum absolute value in waveform baseline
        % should not exceed this fraction of the waveform's abolute peak value (was 0.3)

    % distance metrics
    param.minIsoD   = 4; % minimum isolation distance value
    param.lratioMax = 0.1; % maximum l-ratio value
    param.ssMax     = NaN; % minimum silhouette score 

    % other classification params
    param.minAmplitude      = 20; % in uV
    param.maxRPVviolations  = 0.2; % fraction
    param.RPV_tauR_estimate = NaN;
    param.maxPercSpikesMissing = 25; % in percentage
    param.minNumSpikes      = 500; % number of spikes
    param.maxDrift          = 100;
    param.minPresenceRatio  = 0.8;
    param.minSNR            = 9;
