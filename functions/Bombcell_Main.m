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
% Jesus 08.04.2023

param = struct; % initialize bombcell param structure.
path = struct; % initialize bombcell path structure.

%% Non-existing config file in adequate folder.
if ~isfile(fullfile(input.analysisCode, 'bombcellConfig.m'))
    warning('Config File not found under expected folder ''analysisCode''. Using a default version.')
    % If exists, use the standard one stored within the toolbox.
    if isfile(fullfile(input.toolbox, '\Instructions\bombcellConfig.m'))
        copyfile(fullfile(input.toolbox, '\Instructions\bombcellConfig.m'), input.analysisCode);
        copyfile(fullfile(input.toolbox, '\Instructions\bombcellConfig_KS4.m'), input.analysisCode);
    else
        % It does not exist for some reason.
        error('Could not find the default configuration file for Kilosort. Skipped.')
    end
end

% Valid file found. Run it. Set the quality metric parameters, config function.
if  opt.bombcell == 2
    run(fullfile(input.analysisCode, 'bombcellConfig.m'));
elseif  opt.bombcell == 4
    run(fullfile(input.analysisCode, 'bombcellConfig_KS4.m'));
end

%% Override params based on opts
if isfield(opt,'rerun'),                param.rerun = opt.rerun;  end 

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

%% view units + quality metrics in GUI 
% % load data for GUI
% loadRawTraces = 0; % default: don't load in raw data (this makes the GUI significantly faster)
% bc_loadMetricsForGUI;
% 
% % GUI guide: 
% % left/right arrow: toggle between units 
% % g : go to next good unit 
% % m : go to next multi-unit 
% % n : go to next noise unit
% % up/down arrow: toggle between time chunks in the raw data
% % u: brings up a input dialog to enter the unit you want to go to
% 
% % currently this GUI works best with a screen in portrait mode - we are
% % working to get it to handle screens in landscape mode better. 
% unitQualityGuiHandle = bc_unitQualityGUI(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
%     param, probeLocation, unitType, loadRawTraces);

%% example: get the quality metrics for one unit
% this is an example to get the quality metric for the unit with the
% original kilosort and phy label of xx (0-indexed), which corresponds to
% the unit with qMetric.clusterID == xx + 1, and to
% qMetric.phy_clusterID == xx . This is *NOT NECESSARILY* the
% (xx + 1)th row of the structure qMetric - some of the  clusters that kilosort
% outputs are empty, because they were dropped in the last stages of the
% algorithm. These empty clusters are not included in the qMetric structure
% there are two ways to do this: 
% % 1:
% original_id_we_want_to_load = 0;
% id_we_want_to_load_1_indexed = original_id_we_want_to_load + 1; 
% number_of_spikes_for_this_cluster = qMetric.nSpikes(qMetric.clusterID == id_we_want_to_load_1_indexed);
% % or 2:
% original_id_we_want_to_load = 0;
% number_of_spikes_for_this_cluster = qMetric.nSpikes(qMetric.phy_clusterID == original_id_we_want_to_load);

%% example: get unit labels 
% % the output of `unitType = bc_getQualityUnitType(param, qMetric);` gives
% % the unitType in a number format. 1 indicates good units, 2 indicates mua units, 3
% % indicates non-somatic units and 0 indciates noise units (see below) 
%  
% goodUnits = unitType == 1;
% muaUnits = unitType == 2;
% noiseUnits = unitType == 0;
% nonSomaticUnits = unitType == 3; 
% 
% % example: get all good units number of spikes
% all_good_units_number_of_spikes = qMetric.nSpikes(goodUnits);
% 
% % (for use with another language: output a .tsv file of labels. You can then simply load this) 
% label_table = table(unitType);
% writetable(label_table,[savePath filesep 'templates._bc_unit_labels.tsv'],'FileType', 'text','Delimiter','\t');  
      
%% optional: additionally compute ephys properties for each unit and classify cell types 
% rerunEP = 0;
% region = ''; % options include 'Striatum' and 'Cortex'
% [ephysProperties, unitClassif] = bc_ephysPropertiesPipeline(ephysKilosortPath, savePath, rerunEP, region);
% 
% % example: get good MSN units 
% goodMSNs = strcmp(unitClassif, 'MSN') & unitType == 1; 