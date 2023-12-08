function [param, qMetric, fractionRPVs_allTauR] = bc_loadSavedMetrics(path, saveTSV)
% JF, Load saved quality metrics
% Inputs
% 
% Outputs

qMetric = parquetread(fullfile(path.savepath, 'templates._bc_qMetrics.parquet'));
fractionRPVs_allTauR = parquetread([fullfile(path.savepath, 'templates._bc_fractionRefractoryPeriodViolationsPerTauR.parquet')]);
param = parquetread([fullfile(path.savepath, '_bc_parameters._bc_qMetrics.parquet')]);

% optionally, also save output as a .tsv file
if nargin < 2
    saveTSV = 1; 
end

if saveTSV
    fieldsToSave = {'percentageSpikesMissing_gaussian', ...
        'presenceRatio', 'maxDriftEstimate', 'nPeaks', 'nTroughs', 'isSomatic','waveformDuration_peakTrough', 'spatialDecaySlope','waveformBaselineFlatness',...
        'signalToNoiseRatio','fractionRPVs_estimatedTauR' };
    fieldsRename = {'%_spikes_missing', 'presence_ratio', 'max_drift', 'n_peaks', 'n_troughs', ...
        'is_somatic','waveform_dur', 'spatial_decay_slope','wv_baseline_flatness',...
        'SNR','frac_RPVs' };
     
    cluster_id_vector = qMetric.clusterID - 1; % from bombcell to phy nomenclature 
    if isfield(param,'ephysKilosortPath')
        saveTSV_path = param.ephysKilosortPath;
    else
        saveTSV_path = path.savepath;
    end
    for fid = 1:length(fieldsToSave)
        cluster_table = table(cluster_id_vector, qMetric.(fieldsToSave{fid}), 'VariableNames', {'cluster_id', fieldsRename{fid}});
        writetable(cluster_table,[saveTSV_path filesep 'cluster_' fieldsRename{fid} '.tsv'],'FileType', 'text','Delimiter','\t');  
    end
    
end

end