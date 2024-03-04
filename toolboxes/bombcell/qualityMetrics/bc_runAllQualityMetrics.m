function [qMetric, unitType] = bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, path)
% Inputs
% param: parameter structure. See bc_qualityParamValues for all fields
%   anf information about them.
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
% channelPositions: nChannels x 2 double matrix corresponding to the x and
%   z locations of each channel on the probe, in um
% savePath: sting defining the path where to save bombcell's output
%
% Outputs
% qMetric: structure with fields:
%   percentageSpikesMissing : a gaussian is fit to the spike amplitudes with a
%       'cutoff' parameter below which there are no spikes to estimate the
%       percentage of spikes below the spike-sorting detection threshold - will
%       slightly underestimate in the case of 'bursty' cells with burst
%       adaptation (eg see Fig 5B of Harris/Buzsaki 2000 DOI: 10.1152/jn.2000.84.1.401)
%   fractionRefractoryPeriodViolations: percentage of false positives, ie spikes within the refractory period
%       defined by param.tauR of another spike. This also excludes
%       duplicated spikes that occur within param.tauC of another spike.
%   useTheseTimes : param.computeTimeChunks, this defines the time chunks
%       (deivding the recording in time of chunks of param.deltaTimeChunk size)
%       where the percentage of spike missing and percentage of false positives
%       is below param.maxPercSpikesMissing and param.maxRPVviolations
%   nSpikes : number of spikes for each unit
%   nPeaks : number of detected peaks in each units template waveform
%   nTroughs : number of detected troughs in each units template waveform
%   isSomatic : a unit is defined as Somatic of its trough precedes its main
%       peak (see Deligkaris/Frey DOI: 10.3389/fnins.2016.00421)
%   rawAmplitude : amplitude in uV of the units mean raw waveform at its peak
%       channel. The peak channel is defined by the template waveform.
%   spatialDecay : gets the minumum amplitude for each unit 5 channels from
%       the peak channel and calculates the slope of this decrease in amplitude.
%   isoD : isolation distance, a measure of how well a units spikes are seperate from
%       other nearby units spikes
%   Lratio : l-ratio, a similar measure to isolation distance. see
%       Schmitzer-Torbert/Redish 2005  DOI: 10.1016/j.neuroscience.2004.09.066
%       for a comparison of l-ratio/isolation distance
%   silhouetteScore : another measure similar ti isolation distance and
%       l-ratio. See Rousseeuw 1987 DOI: 10.1016/0377-0427(87)90125-7)
% unitType: nUnits x 1 vector indicating whether each unit met the
%   threshold criterion to be classified as a single unit (1), noise
%   (0) or multi-unit (2)
%
% Modified by Jesus 29/11/2023

% Initialize structures
qMetric = struct;
forGUI = struct;

%% Check parameter values
param = bc_checkParameterFields(param);

%% Get unit max channels
maxChannels = bc_getWaveformMaxChannel(templateWaveforms);

%% Extract and save or load in raw waveforms
[rawWaveformsFull, rawWaveformsPeakChan, signalToNoiseRatio] = ...
    bc_extractRawWaveformsFast(param, spikeTimes_samples, spikeTemplates, path); % it may take time

%% Remove duplicate spikes
[uniqueTemplates, ~, spikeTimes_samples, spikeTemplates, templateAmplitudes, pcFeatures, ...
 rawWaveformsFull, rawWaveformsPeakChan, signalToNoiseRatio, qMetric.maxChannels] = ...
    bc_removeDuplicateSpikes(spikeTimes_samples, spikeTemplates, templateAmplitudes, pcFeatures, rawWaveformsFull, ...
                             rawWaveformsPeakChan, signalToNoiseRatio, maxChannels, param, path);

%% Divide recording into time chunks
spikeTimes_seconds = spikeTimes_samples ./ param.ephys_sample_rate; % convert to seconds after using sample indices to extract raw waveforms

if param.computeTimeChunks
    timeChunks = [min(spikeTimes_seconds):param.deltaTimeChunk:max(spikeTimes_seconds), max(spikeTimes_seconds)];
else
    timeChunks = [min(spikeTimes_seconds), max(spikeTimes_seconds)];
end

%% Calculate quality metrics for spikes
fprintf('\n Extracting quality metrics from %s ... \n', path.rawFile)
for iUnit = 1:size(uniqueTemplates, 1)
    clearvars thisUnit theseSpikeTimes theseAmplis theseSpikeTemplates

    % Get this unit's attributes
    thisUnit = uniqueTemplates(iUnit);
    qMetric.phy_clusterID(iUnit) = thisUnit - 1; % this is the cluster ID as it appears in phy
    qMetric.clusterID(iUnit) = thisUnit; % this is the cluster ID as it appears in phy, 1-indexed (adding 1)

    theseSpikeTimes = spikeTimes_seconds(spikeTemplates == thisUnit);
    theseAmplis = templateAmplitudes(spikeTemplates == thisUnit);

    % Percentage spikes missing (false negatives)
    [percentageSpikesMissing_gaussian, percentageSpikesMissing_symmetric, ksTest_pValue, ~, ~, ~] = ...
        bc_percSpikesMissing(theseAmplis, theseSpikeTimes, timeChunks, param);

    if param.plotDetails
        mkdir(fullfile(path.savePath,"singleplots"))
        txt = [int2str(thisUnit) + "perc_spk_miss"];
        saveas(gcf, fullfile(path.savePath,"singleplots", txt), 'png');
        close(gcf);
    end

    % Fraction contamination (false positives)
    tauR_window = param.tauR_valuesMin:param.tauR_valuesStep:param.tauR_valuesMax;

    [fractionRPVs, ~, ~] = ...
        bc_fractionRPviolations(theseSpikeTimes, theseAmplis, tauR_window, timeChunks, param);

    if param.plotDetails
        txt = [int2str(thisUnit) + "_fract_RP_viol"];
        saveas(gcf, fullfile(path.savePath,"singleplots", txt), 'png');
        close(gcf);
    end

    % Define timechunks to keep: keep times with low percentage spikes missing and low fraction contamination
    [theseSpikeTimes, theseAmplis, theseSpikeTemplates, qMetric.useTheseTimesStart(iUnit), qMetric.useTheseTimesStop(iUnit), ...
     qMetric.RPV_tauR_estimate(iUnit)] = ...
        bc_defineTimechunksToKeep(percentageSpikesMissing_gaussian, fractionRPVs, theseAmplis, ...
                                  theseSpikeTimes, spikeTemplates, timeChunks, param); %QQ add kstest thing, symmetric ect

    % Re-compute percentage spikes missing and fraction contamination on timechunks
    thisUnits_timesToUse = [qMetric.useTheseTimesStart(iUnit), qMetric.useTheseTimesStop(iUnit)];

    [qMetric.percentageSpikesMissing_gaussian(iUnit), qMetric.percentageSpikesMissing_symmetric(iUnit), ...
     qMetric.ksTest_pValue(iUnit), forGUI.ampliBinCenters{iUnit}, forGUI.ampliBinCounts{iUnit}, ...
     forGUI.ampliGaussianFit{iUnit}] = ...
        bc_percSpikesMissing(theseAmplis, theseSpikeTimes, thisUnits_timesToUse, param);

    if param.plotDetails
        txt = [int2str(thisUnit) + "miss_spk"];
        saveas(gcf, fullfile(path.savePath,"singleplots", txt), 'png');
        close(gcf);
    end

    [qMetric.fractionRPVs(iUnit, :), ~, ~] = ...
        bc_fractionRPviolations(theseSpikeTimes, theseAmplis, tauR_window, thisUnits_timesToUse, param); %qMetric.RPV_tauR_estimate(iUnit),

    if param.plotDetails
        txt = [int2str(thisUnit) + "_fract_viol"];
        saveas(gcf, fullfile(path.savePath,"singleplots", txt), 'png');
        close(gcf);
    end

    % Presence ratio (potential false negatives)
    [qMetric.presenceRatio(iUnit)] = ...
        bc_presenceRatio(theseSpikeTimes, theseAmplis, qMetric.useTheseTimesStart(iUnit), ...
                         qMetric.useTheseTimesStop(iUnit), param);

    if param.plotDetails
        txt = [int2str(thisUnit) + "_pres_ratio"];
        saveas(gcf, fullfile(path.savePath,"singleplots", txt), 'png');
        close(gcf);
    end    

    % Maximum cumulative drift estimate
    [qMetric.maxDriftEstimate(iUnit), qMetric.cumDriftEstimate(iUnit)] = ...
        bc_maxDriftEstimate(pcFeatures, pcFeatureIdx, theseSpikeTemplates, theseSpikeTimes, ...
                            channelPositions(:, 2), thisUnit, param);

    % Number of spikes
    qMetric.nSpikes(iUnit) = bc_numberSpikes(theseSpikeTimes);

    % Waveform
    [qMetric.nPeaks(iUnit), qMetric.nTroughs(iUnit), qMetric.isSomatic(iUnit), forGUI.peakLocs{iUnit}, ...
     forGUI.troughLocs{iUnit}, qMetric.waveformDuration_peakTrough(iUnit), forGUI.spatialDecayPoints(iUnit, :), ...
     qMetric.spatialDecaySlope(iUnit), qMetric.waveformBaselineFlatness(iUnit), forGUI.tempWv(iUnit, :)] = ...
        bc_waveformShape(templateWaveforms, thisUnit, qMetric.maxChannels(thisUnit), channelPositions, param);

    if param.plotDetails
        txt = [int2str(thisUnit) + "_wavf_shape"];
        saveas(gcf, fullfile(path.savePath,"singleplots", txt), 'png');
        close(gcf);
    end
    
    % Amplitude
    if param.extractRaw
        rawWaveformsThisUnit = rawWaveformsFull(iUnit, rawWaveformsPeakChan(iUnit), :);
        qMetric.rawAmplitude(iUnit) = bc_getRawAmplitude(rawWaveformsThisUnit, param);

    else
        qMetric.rawAmplitude(iUnit) = NaN;
        qMetric.signalToNoiseRatio(iUnit) = NaN;
    end

    % Distance metrics
    if param.computeDistanceMetrics
        [qMetric.isoD(iUnit), qMetric.Lratio(iUnit), qMetric.silhouetteScore(iUnit), ...
         forGUI.d2_mahal{iUnit}, forGUI.mahalobnis_Xplot{iUnit}, forGUI.mahalobnis_Yplot{iUnit}] = ...
            bc_getDistanceMetrics(pcFeatures, pcFeatureIdx, thisUnit, sum(spikeTemplates == thisUnit), spikeTemplates == thisUnit, ...
                                  theseSpikeTemplates, param);            
    end

    if param.plotDetails
        txt = [int2str(thisUnit) + "_dist_metr"];
        saveas(gcf, fullfile(path.savePath,"singleplots", txt), 'png');
        close(gcf);
    end    

    %% Progress
    if ((mod(iUnit, 50) == 0) || iUnit == length(uniqueTemplates)) && param.verbose
        fprintf(['\n   Finished ', num2str(iUnit), ' / ', num2str(length(uniqueTemplates)), ' units.']);
    end

end

%% Tag unit types and save data
qMetric.maxChannels = qMetric.maxChannels(uniqueTemplates)';
if param.extractRaw
    qMetric.signalToNoiseRatio = signalToNoiseRatio';
end

fprintf('\n Finished extracting quality metrics')

qMetric = bc_saveQMetrics(param, qMetric, forGUI, path);
fprintf('\n Saved quality metrics to %s \n', path.savePath)

unitType = bc_getQualityUnitType(param, qMetric, path);

% Plot quality summary
bc_plotGlobalQualityMetric(qMetric, param, unitType, uniqueTemplates, forGUI.tempWv, path);

end
