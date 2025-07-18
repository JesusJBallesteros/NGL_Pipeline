function Bombcell_Main(~, opt)
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
% Jesus 21.05.2025

%% Paths and other
% Find .bin files.
ephysKilosortPath  = fullfile(opt.FolderProcDataMat, 'kilosort4'); % the raw data binary file is in this folder (for current subject and session)
ephysRawDir        = dir([opt.FolderProcDataMat, '\*.*bin']); % your raw .bin data
savePath           = fullfile(ephysKilosortPath, 'bombcell'); % where you want to save the quality metrics
ephysMetaDir       = []; % path to your .meta or .oebin meta file
% Detect whether data is compressed. Decompress locally, if necessary.
ephysRawFile    = [ephysRawDir.folder, filesep, ephysRawDir.name]; % Ours is never .cbin, so far.

% Version
kilosortVersion = 4; % if using kilosort4 
gain_to_uV      = 0.195; % for DEUTERON, make sure you have this modified in your config file

%% Load default parameters
param = bc.qm.qualityParamValues(ephysMetaDir, ephysRawFile, ephysKilosortPath, gain_to_uV, kilosortVersion);

%% Override params based on opts
run("bombcellConfig.m"); % THIS OVERRIDES THE PREVIOUS 'param' CALL. Try Defaults first
opt.callBcGUI = 1;

%% Faster compute. Compile .mex file only if not done yet
if ~isfile('C:\Code\ephys-data-pipeline\toolboxes\bombcell\matlab\+bc\+ep\+helpers\CCGHeart.mexw64')
    orig = pwd;
    cd('C:\Code\ephys-data-pipeline\toolboxes\bombcell\matlab\+bc\+ep\+helpers');
    mex -O CCGHeart.c 
    cd(orig); clear orig
end

%% Load data from Kilosort outputs.
[spikeTimes_samples, spikeClusters, templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions] ...
     = bc.load.loadEphysData(ephysKilosortPath, savePath);

%% Compute quality metrics 
[qMetric, unitType] = bc.qm.runAllQualityMetrics(param, spikeTimes_samples, spikeClusters, ...
                    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);

%% view units + quality metrics in GUI 
if opt.callBcGUI
    % load data for GUI
    loadRawTraces = 0; % default: don't load in raw data (this makes the GUI significantly faster)
    bc.load.loadMetricsForGUI;
    
    % GUI guide: 
    % left/right arrow: toggle between units 
    % g : go to next good unit 
    % m : go to next multi-unit 
    % n : go to next noise unit
    % up/down arrow: toggle between time chunks in the raw data
    % u: brings up a input dialog to enter the unit you want to go to
    unitQualityGuiHandle = bc.viz.unitQualityGUI_synced(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
        param, probeLocation, unitType, loadRawTraces);
    
    GUIdlg = warndlg('Close this dialog to continue.', 'BombCell GUI Dialog.');
    waitfor(GUIdlg);
    disp('GUI Dialog closed.');
end

%% example: get the quality metrics for one unit
% % this is an example to get the quality metric for the unit with the
% % original kilosort and phy label of xx (0-indexed), which corresponds to
% % the unit with qMetric.clusterID == xx + 1, and to
% % qMetric.phy_clusterID == xx . This is *NOT NECESSARILY* the
% % (xx + 1)th row of the structure qMetric - some of the  clusters that kilosort
% % outputs are empty, because they were dropped in the last stages of the
% % algorithm. These empty clusters are not included in the qMetric structure
% % there are two ways to do this: 
% % 1:
% % original_id_we_want_to_load = 0;
% % id_we_want_to_load_1_indexed = original_id_we_want_to_load + 1; 
% % number_of_spikes_for_this_cluster = qMetric.nSpikes(qMetric.clusterID == id_we_want_to_load_1_indexed);
% 
% % or 2:
% % original_id_we_want_to_load = 0;
% % number_of_spikes_for_this_cluster = qMetric.nSpikes(qMetric.phy_clusterID == original_id_we_want_to_load);

%% example: get unit labels 
% % % the output of `unitType = bc_getQualityUnitType(param, qMetric);` gives
% % % the unitType in a number format. 1 indicates good units, 2 indicates mua units, 3
% % % indicates non-somatic units and 0 indciates noise units (see below) 
% goodUnits = unitType == 1;
% muaUnits = unitType == 2;
% noiseUnits = unitType == 0;
% nonSomaticUnits = unitType == 3; 
 
%% example: get all good units number of spikes
% all_good_units_number_of_spikes = qMetric.nSpikes(goodUnits);
% 
% % (for use with another language: output a .tsv file of labels. You can then simply load this) 
% label_table = table(unitType);
% writetable(label_table,[savePath filesep 'templates._bc_unit_labels.tsv'],'FileType', 'text','Delimiter','\t');  

%% Optionally get ephys properties for your cell. Bombcell will also attempt to classify your data if it is (a) from the cortex or striatum and (b) you specify this in the "region" variable.
rerunEP = 0;
region = ''; % options include 'Striatum' and 'Cortex'
[ephysProperties, unitClassif] = bc.ep.runAllEphysProperties(ephysKilosortPath, savePath, rerunEP, region);

save(fullfile(savePath, 'ephysProperties.mat'), 'ephysProperties', 'unitClassif');