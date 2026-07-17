function sfc = spikeFieldCoupling(FT_data, spike, opt)
% spikeFieldCoupling  Spike-LFP phase locking (PPC) + spike-field coherence.
%
% PURPOSE:
%   For every curated cluster in `spike` with enough events, compute:
%     (1) PPC per (cluster, frequency) via ft_spiketriggeredspectrum +
%         ft_spiketriggeredspectrum_stat. Vinck et al. 2010/2011 PPC
%         estimator; default 'ppc2' (robust to spike-rate dependencies).
%     (2) Spike-field coherence per (cluster, LFP channel, frequency)
%         via ft_freqanalysis on a joined spike+LFP FieldTrip struct,
%         then ft_connectivityanalysis with cfg.method = 'coh',
%         cfg.complex = 'absimag' to suppress volume-conduction bleed.
%
% USAGE:
%   sfc = spikeFieldCoupling(FT_data, spike, opt);
%
% INPUTS:
%   FT_data - continuous FieldTrip struct. Should carry .chanArea
%             (backfilled at load time in NGL07).
%   spike   - spike struct as produced by loadSpikes (fields .label,
%             .timeStamps or .spikeTimes; here we accept .spikeTimes
%             cell {nClusters x 1} of double vectors in SECONDS from
%             start of the recording, matching the FT time base).
%   opt     - resolved options. Consumes:
%               .lfp.bands              cell (default theta/beta/gamma)
%               .lfp.tfrAreaFilter      restrict LFP channels
%               .lfp.spikeField.minSpikes    default 50
%               .lfp.spikeField.ppcMethod    'ppc0'|'ppc1'|'ppc2'
%                                             (default 'ppc2')
%               .lfp.spikeField.timwin       [s] STA window around
%                                             each spike (default 0.5)
%               .lfp.spikeField.foi          freq vector for coherence
%                                             (default 2:2:100)
%
% OUTPUT (struct):
%   .clusters   - {nClustProcessed x 1} char cluster IDs
%   .nSpikes    - [nClustProcessed x 1] spike counts used
%   .channels   - {nCh x 1} LFP channel labels used
%   .chanArea   - {nCh x 1} LFP channel area tags
%   .foi        - [1 x nFreq] frequency vector
%   .ppc        - [nClust x nFreq] PPC per cluster per frequency
%                  (nan where the underlying stat could not be computed)
%   .coherence  - [nClust x nCh x nFreq] spike-field coherence
%                  (absolute imaginary coherence, [0 1])
%   .method     - {ppcMethod, timwin, foi}
%   .provenance - buildLFPProvenance snapshot
%
% NOTES:
%   * The two computations use different FieldTrip paths but the same
%     input. PPC is per-cluster + per-frequency (channel-marginalized
%     inside ft_spiketriggeredspectrum_stat when we pass all channels);
%     coherence keeps per-channel resolution so cross-area effects
%     (e.g. NCL-cluster ↔ STR-LFP) show up in the output cube.
%   * Cross-area analyses fall out for free from the chanArea tags -
%     downstream slicing by [ismember(chanArea, 'STR')] on the third
%     dim of .coherence gives NCL-cluster ↔ STR-LFP coherence, and
%     vice versa.
%
% SEE ALSO:
%   ft_spiketriggeredspectrum, ft_spiketriggeredspectrum_stat,
%   ft_freqanalysis, ft_connectivityanalysis, NGL07_LFPanalysis.
%
% Last modified 26.06.2026 (Jesus) - new (LFP Pass 3).

    minSpikes = localOptField(opt, {'lfp','spikeField','minSpikes'}, 50);
    ppcMethod = localOptField(opt, {'lfp','spikeField','ppcMethod'}, 'ppc2');
    timwin    = localOptField(opt, {'lfp','spikeField','timwin'},    0.5);
    foi       = localOptField(opt, {'lfp','spikeField','foi'},       2:2:100);

    %% Optional per-area LFP subselection.
    keepChanIdx = 1:numel(FT_data.label);
    if isfield(opt,'lfp') && isfield(opt.lfp,'tfrAreaFilter') ...
            && ~isempty(opt.lfp.tfrAreaFilter) && isfield(FT_data,'chanArea')
        keepChanIdx = find(ismember(FT_data.chanArea, cellstr(opt.lfp.tfrAreaFilter)));
        if isempty(keepChanIdx)
            error('NGL:spikeFieldCoupling:noChans', ...
                'opt.lfp.tfrAreaFilter matched 0 LFP channels.');
        end
    end
    selCfg = []; selCfg.channel = FT_data.label(keepChanIdx);
    FT_lfp = ft_selectdata(selCfg, FT_data);
    nCh    = numel(FT_lfp.label);

    %% Build a FieldTrip spike struct from `spike`.
    % NGL02_postPhy loadSpikes returns spike.spikeTimes as a cell of
    % vectors in seconds (matches loadSpikes.m; see its docstring).
    assert(isfield(spike, 'spikeTimes') && iscell(spike.spikeTimes), ...
        'NGL:spikeFieldCoupling:noSpikeTimes', ...
        'spike.spikeTimes cell (per-cluster spike-time vectors in seconds) is required.');
    nClust = numel(spike.spikeTimes);

    ftSpike           = struct();
    ftSpike.label     = cell(0, 1);
    ftSpike.timestamp = cell(0, 1);
    keepClust         = false(nClust, 1);
    for c = 1:nClust
        st = spike.spikeTimes{c};
        st = st(~isnan(st));
        if numel(st) < minSpikes, continue; end
        keepClust(c) = true;
        if isfield(spike, 'label') && numel(spike.label) >= c
            ftSpike.label{end+1, 1} = char(string(spike.label{c}));
        else
            ftSpike.label{end+1, 1} = sprintf('clust%03d', c);
        end
        ftSpike.timestamp{end+1, 1} = st(:) * FT_lfp.fsample;   % convert s -> samples
    end
    nProc = sum(keepClust);
    if nProc == 0
        warning('NGL:spikeFieldCoupling:noClusters', ...
            'No clusters with >= %d spikes; nothing to compute.', minSpikes);
        sfc = localEmptyOutput(FT_lfp, foi, minSpikes, ppcMethod, timwin, opt);
        return
    end
    ftSpike.trialtime = [FT_lfp.time{1}(1), FT_lfp.time{1}(end)];
    ftSpike.cfg       = struct();

    %% (1) PPC: per-spike phase spectrum, then Vinck estimator.
    cfg          = [];
    cfg.method   = 'convol';
    cfg.foi      = foi;
    cfg.timwin   = [-timwin timwin];
    cfg.channel  = FT_lfp.label;
    cfg.spikechannel = ftSpike.label;
    cfg.taper    = 'hanning';
    sts = ft_spiketriggeredspectrum(cfg, FT_lfp, ftSpike);

    ppc = nan(nProc, numel(foi));
    for k = 1:nProc
        cfgStat            = [];
        cfgStat.method     = ppcMethod;
        cfgStat.spikechannel = ftSpike.label{k};
        cfgStat.channel    = FT_lfp.label;
        cfgStat.avgoverchan = 'weighted';
        try
            stat = ft_spiketriggeredspectrum_stat(cfgStat, sts);
            % .ppc0 / .ppc1 / .ppc2 field per selected method
            fld = intersect({'ppc0','ppc1','ppc2'}, fieldnames(stat));
            if ~isempty(fld), ppc(k, :) = squeeze(stat.(fld{1})); end
        catch ME
            warning('NGL:spikeFieldCoupling:ppcFail', ...
                'PPC failed for cluster %s: %s', ftSpike.label{k}, ME.message);
        end
    end

    %% (2) Spike-field coherence: encode spikes as an LFP-rate 0/1 series
    %       joined with the LFP, then coh via ft_connectivityanalysis.
    coherence = nan(nProc, nCh, numel(foi));
    try
        cfgConv                 = [];
        cfgConv.fsample         = FT_lfp.fsample;
        spikeAsLFP              = ft_spike_convert2fieldtrip(cfgConv, ftSpike);
        % Concatenate spike channels + LFP channels into one struct.
        cfgApp                  = [];
        joined                  = ft_appenddata(cfgApp, spikeAsLFP, FT_lfp);
        cfgFreq                 = [];
        cfgFreq.method          = 'mtmfft';
        cfgFreq.output          = 'fourier';
        cfgFreq.foi             = foi;
        cfgFreq.taper           = 'dpss';
        cfgFreq.tapsmofrq       = 4;
        cfgFreq.keeptrials      = 'yes';
        cfgFreq.channel         = 'all';
        freq                    = ft_freqanalysis(cfgFreq, joined);
        cfgCoh                  = [];
        cfgCoh.method           = 'coh';
        cfgCoh.complex          = 'absimag';
        cfgCoh.channelcmb       = localMakeChanCombos(ftSpike.label, FT_lfp.label);
        c                       = ft_connectivityanalysis(cfgCoh, freq);
        coherence               = localReshapeCoherence(c, ftSpike.label, FT_lfp.label, foi);
    catch ME
        warning('NGL:spikeFieldCoupling:cohFail', ...
            'Coherence path failed: %s. PPC-only output.', ME.message);
    end

    %% Package.
    sfc.clusters   = ftSpike.label;
    sfc.nSpikes    = cellfun(@numel, ftSpike.timestamp);
    sfc.channels   = FT_lfp.label;
    if isfield(FT_lfp,'chanArea')
        sfc.chanArea = FT_lfp.chanArea;
    else
        sfc.chanArea = repmat({'main'}, nCh, 1);
    end
    sfc.foi        = foi;
    sfc.ppc        = ppc;
    sfc.coherence  = coherence;
    sfc.method     = struct( ...
        'ppcMethod',  ppcMethod, ...
        'timwin_s',   timwin, ...
        'foi',        foi, ...
        'minSpikes',  minSpikes);
    sfc.provenance = buildLFPProvenance('', opt, struct( ...
        'analysis',   'spikeFieldCoupling', ...
        'nClusters',  nProc, ...
        'nChans',     nCh));
end


% =======================================================================
function cmb = localMakeChanCombos(spikeLabels, lfpLabels)
% All (spike, lfp) channel pairs for ft_connectivityanalysis.
    [S, L] = ndgrid(1:numel(spikeLabels), 1:numel(lfpLabels));
    cmb = [spikeLabels(S(:)), lfpLabels(L(:))];
end

function cube = localReshapeCoherence(c, spikeLabels, lfpLabels, foi)
% Reshape ft_connectivityanalysis output (cohspctrm: [nCmb x nFreq]) into
% [nClust x nCh x nFreq]. Missing combos -> NaN.
    nClust = numel(spikeLabels);
    nCh    = numel(lfpLabels);
    cube   = nan(nClust, nCh, numel(foi));
    if ~isfield(c, 'cohspctrm') || ~isfield(c, 'labelcmb')
        return
    end
    for i = 1:size(c.labelcmb, 1)
        s = find(strcmp(c.labelcmb{i, 1}, spikeLabels), 1);
        l = find(strcmp(c.labelcmb{i, 2}, lfpLabels),   1);
        if ~isempty(s) && ~isempty(l)
            cube(s, l, :) = c.cohspctrm(i, :);
        end
    end
end

function sfc = localEmptyOutput(FT_lfp, foi, minSpikes, ppcMethod, timwin, opt)
    sfc = struct();
    sfc.clusters   = {};
    sfc.nSpikes    = [];
    sfc.channels   = FT_lfp.label;
    if isfield(FT_lfp,'chanArea')
        sfc.chanArea = FT_lfp.chanArea;
    else
        sfc.chanArea = repmat({'main'}, numel(FT_lfp.label), 1);
    end
    sfc.foi        = foi;
    sfc.ppc        = zeros(0, numel(foi));
    sfc.coherence  = zeros(0, numel(FT_lfp.label), numel(foi));
    sfc.method     = struct('ppcMethod', ppcMethod, 'timwin_s', timwin, ...
                            'foi', foi, 'minSpikes', minSpikes);
    sfc.provenance = buildLFPProvenance('', opt, struct('analysis','spikeFieldCoupling', 'empty', true));
end

function v = localOptField(opt, path, dflt)
    v = dflt;  cursor = opt;
    for k = 1:numel(path)
        if isstruct(cursor) && isfield(cursor, path{k})
            cursor = cursor.(path{k});
        else
            return
        end
    end
    v = cursor;
end
