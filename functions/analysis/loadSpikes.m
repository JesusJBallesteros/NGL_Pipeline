function [spike] = loadSpikes(opt)
% loadSpikes  Read KS/Phy-curated clusters into a NGL-standard spike struct.
%
% PURPOSE:
%   Reads the curated Kilosort + Phy output for a single area's KS folder
%   (opt.KSfolder), excludes 'noise' clusters by default, and assembles
%   per-cluster timestamps, channel, shank, ROI, and (optionally) raw
%   waveforms into the lab-standard 'spike' struct.
%
% USAGE:
%   spike = loadSpikes(opt)
%
% INPUTS:
%   opt - resolved options struct (post-set_default). Required fields:
%           .KSfolder           Kilosort/Phy output folder for THIS area
%           .FolderProcDataMat  preprocessing folder containing the .bin
%           .SavFileName        session name (used to build .bin filename)
%           .spparams           cluster-loading flags (excludeNoise, loadPCs)
%           .isibins            ISI histogram bin edges, ms
%           .getwF              if true, also extract raw waveforms
%           .gwfparams          waveform-extraction config (wfWin, nWf, ...)
%           .area               area label to stamp on every cluster's
%                               spike.roi (set per-area by NGL02; 'all' in
%                               single-area mode)
%
% OUTPUT:
%   spike - struct with one entry per cluster in 1xNclust cell arrays:
%             .label         cluster ID (char, e.g. '42')
%             .timestamp     per-cluster spike times in SECONDS
%             .depth, .ampl, .templampl
%             .ch            max-amplitude channel from Phy
%             .shank         shank index derived from channel
%             .roi           area tag (opt.area; 'all' in single-area mode)
%             .KSLabel       Kilosort auto-classification ('good'/'mua'/'noise')
%             .bc_unitType   Bombcell classification ('GOOD'/'MUA'/'NOISE'/...);
%                            '' when Bombcell did not run for this session
%             .HumanLabel    Phy-curated 'group' label ('good'/'mua'/'noise'),
%                            i.e. the call the human made during Phy curation
%             .phyLabel      optional free-text annotation set in Phy
%                            (separate from .label, which is the cluster ID);
%                            '' when the cluster_info.tsv has no 'label' column
%             .waveform, .waveFormsMean   only when opt.getwF=true
%             .isihist       ISI histogram (always)
%
% NOTES:
%   - All option fields above are guaranteed present by set_default; this
%     function no longer carries inline defaults.
%   - Multi-area: NGL02 sets opt.KSfolder and opt.area per area, then
%     calls loadSpikes once per area. Every cluster returned from a given
%     call is tagged with opt.area.
%
% Jesus, 29.05.2026

% Waveform extraction option defaults
gwfparams = struct('dataType', 'int16',  ... % Data type of .dat file
                   'wfWin',    opt.gwfparams.wfWin, ... %[-32 63], ... % Number of samples around spiketime to include in waveform.
                   'nWf',      opt.gwfparams.nWf,     ... % N waveforms to extrac tper unit.
                   'dataDir',  fullfile(opt.KSfolder), ... % KiloSort/Phy output folder
                   'fileName', fullfile(opt.FolderProcDataMat, [opt.SavFileName, '.bin']), ... % .dat file containing the raw 
                   'nCh',      [], ... % we don't know yet
                   'spikeTimes', [], 'spikeClusters', []); % Reserved for each cluster

%% Extract data from python files into a matlab friendly matrix
% if ~exist(fullfile(opt.spikeSorted, 'spike.mat'), "file")
    spikes = loadKSdir(opt.KSfolder, opt.spparams); % Helper function from Cortex-lab toolbox
    
    if any(spikes.st <= -(opt.gwfparams.wfWin(1))/spikes.sample_rate)
       idx = spikes.st <= -(opt.gwfparams.wfWin(1))/spikes.sample_rate;
        spikes.st(idx)              = [];
        spikes.spikeTemplates(idx)  = [];
        spikes.clu(idx)             = [];
        spikes.tempScalingAmps(idx) = [];
        spikes.spikeAmps(idx)       = [];
        spikes.spikeDepths(idx)     = [];
    end

    %% Get relevant info
    % Phy2 must have been run before, so this file exists
    try phy_table = readtable(fullfile(opt.KSfolder, 'cluster_info.tsv'), "FileType", "text", 'Delimiter', '\t'); 
    catch ME
        % No 'clusters_info.tsv' file found. Phy needs to be run and clusters ACTUALLY human-labeled.
        error('NGL02:loadSpikes', 'Could not find clusters_info.tsv file: %s', ME.message);
    end
    shanksmap       = [spikes.chshanks, spikes.chmap];
    clusters        = sort(unique(spikes.cids)); % get and sort clusters by id
    nclust          = numel(clusters); % number of clusters

    % Which optional columns does this cluster_info.tsv carry? KSLabel and
    % group are produced by Kilosort/Phy and should always be present;
    % bc_unitType is added by Bombcell (may be absent on legacy sessions);
    % label is an optional Phy free-text annotation. Defensive lookup so
    % older sessions don't error.
    phy_cols    = phy_table.Properties.VariableNames;
    hasKSLabel  = ismember('KSLabel',     phy_cols);
    hasBcType   = ismember('bc_unitType', phy_cols);
    hasGroup    = ismember('group',       phy_cols);
    hasPhyLabel = ismember('label',       phy_cols);

    % for waveforms
    gwfparams.nCh = spikes.n_channels_dat;

    % Allocate memory
    spike.label     = cell(1,nclust); % for clusters IDs
    spike.timestamp = cell(1,nclust); % spikes timestamps
    
    %% Proceed to extract timestamps for each cluster
    disp('Extracting curated clusters from Phy files. If many waveforms are requested, it may take a while.')
    for cl = 1:nclust
        spike.label{cl}      = num2str(clusters(cl));
        spike.timestamp{cl}     = spikes.st(spikes.clu==clusters(cl)); % in seconds
        spike.depth{cl}         = spikes.spikeDepths(spikes.clu==clusters(cl));
        spike.ampl{cl}          = spikes.spikeAmps(spikes.clu==clusters(cl));
        spike.templampl{cl}     = spikes.tempScalingAmps(spikes.clu==clusters(cl));
        
        % find maxChannel using the cluster index of the bombcell table
        spike.ch{cl}            = phy_table.ch(find(phy_table.cluster_id==clusters(cl)));
        % Link it to the shank-ch equivalent
        spike.shank{cl}         = shanksmap(shanksmap(:,2) == spike.ch{cl}, 1);
        % Tag with the area being processed. NGL02 sets opt.area per area
        % in multi-area runs; single-area runs default to 'all'.
        % (Replaces the old opt.mapkey(spike.shank{cl}) lookup, which has
        % been retired in favour of input.Areas / input.areaMap.)
        spike.roi{cl}           = opt.area;

        % Curation labels: pulled from cluster_info.tsv. KSLabel = Kilosort
        % auto label; HumanLabel = the call the human made in Phy (the
        % 'group' column, renamed for clarity); bc_unitType = Bombcell's
        % classification; phyLabel = optional free-text annotation. Any
        % column not present in the TSV becomes ''.
        clusterRow              = find(phy_table.cluster_id == clusters(cl), 1, 'first');
        spike.KSLabel{cl}       = localCellChar(phy_table, 'KSLabel',     hasKSLabel,  clusterRow);
        spike.bc_unitType{cl}   = localCellChar(phy_table, 'bc_unitType', hasBcType,   clusterRow);
        spike.HumanLabel{cl}    = localCellChar(phy_table, 'group',       hasGroup,    clusterRow);
        spike.phyLabel{cl}      = localCellChar(phy_table, 'label',       hasPhyLabel, clusterRow);

        %% Extract waveforms
        if opt.getwF
            % a few more params for 'getWaveForms' dep on cluster
            gwfparams.spikeTimes = ceil(spike.timestamp{cl}*spikes.sample_rate); % Vector of cluster spike times (in samples) same length as .spikeClusters
            gwfparams.spikeClusters = spikes.clu(spikes.clu==clusters(cl));
                       
            % Get waveforms
            wF = getWaveForms(gwfparams);
        
            % Refine
            wF.waveForms = squeeze(wF.waveForms);
            wF.waveFormsMean = squeeze(wF.waveFormsMean);
        
            % Find averaged max amplitude channel
            [wF.maxamplch, ~, ~]    = find(wF.waveFormsMean==min(min(wF.waveFormsMean)));
            spike.waveform{cl}      = permute(wF.waveForms(:,wF.maxamplch,:),[2,3,1]);
            spike.waveform{cl}      = squeeze(spike.waveform{cl});
            spike.waveFormsMean{cl} = wF.waveFormsMean;
         end
    
    end

% Calculate the ISI histogram for all clusters
[spike.isihist] = calc_isihist(spike, opt);

end

function s = localCellChar(T, colName, hasCol, row)
% Safely fetch a single value from a readtable column and return it as a
% char vector. Handles the various shapes readtable can produce
% (categorical, string, cell-of-char, missing, numeric). Returns '' when
% the column does not exist, the row is missing, or the value is empty.
    if ~hasCol || isempty(row)
        s = '';
        return
    end
    val = T.(colName)(row);
    if iscell(val),        val = val{1};   end
    if iscategorical(val), val = char(val); end
    if isstring(val),      val = char(val); end
    if ismissing(val)
        s = '';
        return
    end
    if isnumeric(val) || islogical(val)
        if isempty(val) || all(isnan(double(val(:))))
            s = '';
        else
            s = num2str(val);
        end
        return
    end
    if ischar(val) && (strcmp(val,'<undefined>') || strcmp(val,'<missing>') || strcmp(val,'NaN'))
        s = '';
    elseif isempty(val)
        s = '';
    else
        s = val;
    end
end