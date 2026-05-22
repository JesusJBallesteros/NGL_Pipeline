function Bombcell_Main(input, opt)
% Bombcell_Main  Run Bombcell automatic quality metrics on Kilosort output.
%
% PURPOSE:
%   Called from NGL01_Main stage 05 (when opt.bombcell = true). Loads the
%   Kilosort 4 output and the raw .bin file, runs bc.qm.runAllQualityMetrics
%   (from the Bombcell toolbox), saves QC results to the bombcell/ subfolder,
%   and optionally shows the Bombcell GUI for interactive review.
%
% USAGE:
%   Bombcell_Main(input, opt)
%   Do not call directly; gated by opt.bombcell in NGL01_Main.
%
% INPUTS:
%   input  - struct from set_default; relevant fields:
%              .toolbox             root toolbox path (for mex compilation check)
%   opt    - options struct; relevant fields:
%              .FolderProcDataMat   preprocessing folder (contains .bin file)
%              .KSfolder            Kilosort output folder for this run
%                                   (single-area: kilosort\4\; multi-area: <Area>\)
%
% OUTPUTS:
%   <KSfolder>/bombcell/            Bombcell QC results directory
%     unitType.npy                  unit classification (good/MUA/noise)
%     qMetrics.mat                  full quality metrics struct
%     (+ Bombcell standard output files)
%
% PARAMETERS:
%   Bombcell quality metric thresholds are set in bombcellConfig.m (stored
%   in analysisCode/). Copy the template from the toolbox and adjust
%   per-dataset after reviewing the QC histograms and GUI.
%   On first run, raw waveforms are extracted from the .bin file (~10–20 min
%   extra). On re-runs, pre-extracted waveforms are loaded from disk.
%
% CITE:
%   Bombcell: https://github.com/Julie-Fabre/bombcell
%
% Jesus 12.05.2026

%% Paths and other
% KS output folder: use opt.KSfolder so single-area and multi-area runs both
% resolve correctly. In single-area mode this equals kilosort\4\; in
% multi-area mode it is the per-area sub-folder (e.g. preprocessing\NCL\).
ephysKilosortPath  = opt.KSfolder;
ephysRawDir        = dir([opt.FolderProcDataMat, '\*.*bin']); % raw .bin lives in the preprocessing folder (shared across areas)
savePath           = fullfile(ephysKilosortPath, 'bombcell'); % quality metrics saved alongside the KS output for this area
ephysMetaDir       = []; % path to your .meta or .oebin meta file
% Detect whether data is compressed. Decompress locally, if necessary.
ephysRawFile    = [ephysRawDir.folder, filesep, ephysRawDir.name]; % Ours is never .cbin, so far.

% Version
kilosortVersion = 4; % if using kilosort4 
gain_to_uV      = 0.195; % for DEUTERON, make sure you have this modified in your config file

if exist(fullfile(savePath,'ephysProperties.mat'), 'file')
    disp('Bombcell metrics already existing. To re-run BC, remove the previous data.')
    return
end

%% Load default parameters
param = bc.qm.qualityParamValues(ephysMetaDir, ephysRawFile, ephysKilosortPath, gain_to_uV, kilosortVersion);

%% Override params based on opts
run("bombcellConfig.m"); % THIS OVERRIDES THE PREVIOUS 'param' CALL. Try Defaults first

%% Faster compute. Compile .mex file only if not done yet
if ~isfile([input.toolbox, '\toolboxes\bombcell\matlab\+bc\+ep\+helpers\CCGHeart.mexw64'])
    orig = pwd;
    cd([input.toolbox, '\toolboxes\bombcell\matlab\+bc\+ep\+helpers']);
    mex -O CCGHeart.c 
    cd(orig); clear orig
end

%% Load data from Kilosort outputs.
[spikeTimes_samples, spikeClusters, templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions] ...
     = bc.load.loadEphysData(ephysKilosortPath, savePath);

%% Compute quality metrics 
[qMetric, unitType] = bc.qm.runAllQualityMetrics(param, spikeTimes_samples, spikeClusters, ...
                    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);

%% Tag area of origin (multi-area mode only)
% Derive area label from the last path component of opt.KSfolder.
% opt.KSfolders being present is the reliable multi-area guard (same
% convention used in master_kilosort4 and NGL01_Main).
if isfield(opt, 'KSfolders')
    [~, areaLabel] = fileparts(opt.KSfolder);
else
    areaLabel = '';
end

if ~isempty(areaLabel)
    nUnits = numel(qMetric.clusterID);

    % 1. Add .area field to qMetric struct and re-save qMetrics.mat.
    %    Enables filtering / concatenation in MATLAB analysis (e.g.
    %    strcmp(qMetric.area, 'NCL')).
    qMetric.area = repmat({areaLabel}, nUnits, 1);
    save(fullfile(savePath, 'qMetrics.mat'), 'qMetric');

    % 2. Write cluster_area.tsv to the KS output folder.
    %    Phy auto-loads any cluster_*.tsv it finds there and displays it
    %    as a column in the cluster table — no Phy config changes needed.
    %    cluster_id must be 0-indexed (phy_clusterID).
    T = table(qMetric.phy_clusterID(:), repmat({areaLabel}, nUnits, 1), ...
              'VariableNames', {'cluster_id', 'area'});
    writetable(T, fullfile(opt.KSfolder, 'cluster_area.tsv'), ...
               'FileType', 'text', 'Delimiter', '\t');

    % 3. Write area_label.txt alongside the bombcell .npy files.
    %    Python analysis scripts that load unitType.npy, etc. can read
    %    this one-liner to know which area the arrays belong to.
    fid = fopen(fullfile(savePath, 'area_label.txt'), 'w');
    fprintf(fid, '%s\n', areaLabel);
    fclose(fid);

    fprintf('Area label ''%s'' written to qMetrics.mat, cluster_area.tsv, and area_label.txt (%d units).\n', ...
            areaLabel, nUnits);
end

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