function [stat, info] = lfpClusterStats(A, B, opt, varargin)
% lfpClusterStats  Cluster-permutation test on TFRs, one way for the whole pipeline.
%
% PURPOSE:
%   Time-frequency maps have thousands of bins, so an uncorrected test finds
%   "effects" everywhere. FieldTrip's cluster-permutation test is the standard
%   answer (Maris & Oostenveld 2007), but it takes a dozen cfg fields that have
%   to agree between analyses for their results to be comparable. This wraps it
%   once, so power contrasts, comodulograms and tagging tests are all corrected
%   the same way and say so in their provenance.
%
%   What it returns is a cluster-level claim: where a cluster is significant,
%   the data differ somewhere in it. The individual bins inside are not each
%   significant, and the cluster's edges are not a confidence interval.
%
% USAGE:
%   [stat, info] = lfpClusterStats(A, B, opt)                 % A vs B, per trial
%   [stat, info] = lfpClusterStats(A, [], opt, 'design', 'baseline')
%   [stat, info] = lfpClusterStats(A, B, opt, 'design', 'paired')
%
% INPUTS:
%   A, B - FieldTrip freq structures with trials kept (dimord starting 'rpt'),
%          same channels, frequencies and time axis. B may be [] for the
%          'baseline' design.
%   opt  - options struct; reads opt.lfp.stats.* (see below).
%   Name/value pairs override opt:
%     'design'      'trials' (default) | 'paired' | 'baseline'
%     'numrand'     permutations (default 1000)
%     'alpha'       cluster-level alpha, two-sided by default (0.05)
%     'clusteralpha' bin-level threshold for cluster formation (0.05)
%     'avgoverchan' collapse channels before testing (default true)
%     'latency'     [tmin tmax] to test; default 'all' (see NOTES)
%     'frequency'   [fmin fmax] to test; default 'all'
%     'neighbours'  FieldTrip neighbour structure, needed only to cluster
%                   across channels (default [] = no channel clustering)
%
%   Designs:
%     'trials'   A and B are different trials -> independent-samples t.
%                The session-level default: trials of two conditions.
%     'paired'   A and B are the same units in two conditions (sessions or
%                subjects at group level) -> dependent-samples t.
%     'baseline' A against its own baseline: each trial's baseline mean is
%                held against its time course, dependent-samples t. Use it on
%                RAW power, not on a baseline-normalised TFR - normalising
%                first makes the comparison circular.
%
% OUTPUT:
%   stat - ft_freqstatistics output: .stat, .prob, .mask, .posclusters, ...
%          .mask is the significance map to hand to plotTFRpanel.
%   info - what was run: design, numbers, counts of surviving clusters and
%          their p-values. Put it in the provenance of whatever you save.
%
% NOTES:
%   * Latency: leaving the baseline inside the tested window makes a
%     'baseline' design partly circular and weakens a 'trials' contrast by
%     adding bins where nothing can differ. Set 'latency' to the post-event
%     window when that is the question.
%   * With avgoverchan false and no neighbours, clusters cannot grow across
%     channels; each channel is clustered in time-frequency alone.
%   * numrand 1000 resolves p down to 0.001. Below ~500 the p-values are too
%     grainy to interpret near alpha.
%   * FieldTrip's "Not all replications are used for the computation of the
%     statistic" is expected on wavelet TFRs: the edge cone is NaN, so the
%     bins there have fewer trials than the design lists. FieldTrip drops
%     them per bin (nanmean/nanvar) and the test stays valid. Restricting
%     'latency' to the cone-free window silences it.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 0).

    p = inputParser;
    p.addParameter('design',       localOpt(opt, {'lfp','stats','design'}, 'trials'));
    p.addParameter('numrand',      localOpt(opt, {'lfp','stats','numrand'}, 1000));
    p.addParameter('alpha',        localOpt(opt, {'lfp','stats','alpha'}, 0.05));
    p.addParameter('clusteralpha', localOpt(opt, {'lfp','stats','clusteralpha'}, 0.05));
    p.addParameter('clusterstatistic', localOpt(opt, {'lfp','stats','clusterstatistic'}, 'maxsum'));
    p.addParameter('tail',         localOpt(opt, {'lfp','stats','tail'}, 0));
    p.addParameter('avgoverchan',  localOpt(opt, {'lfp','stats','avgoverchan'}, true));
    p.addParameter('minnbchan',    localOpt(opt, {'lfp','stats','minnbchan'}, 0));
    p.addParameter('latency',      localOpt(opt, {'lfp','stats','latency'}, 'all'));
    p.addParameter('frequency',    localOpt(opt, {'lfp','stats','frequency'}, 'all'));
    p.addParameter('neighbours',   []);
    p.parse(varargin{:});
    a = p.Results;
    design = lower(char(a.design));

    assert(isstruct(A) && isfield(A, 'powspctrm'), 'lfpClusterStats:input', ...
        'A must be a FieldTrip freq structure.');
    assert(contains(lower(localDimord(A)), 'rpt'), 'lfpClusterStats:noTrials', ...
        ['A has no trial dimension. Compute the TFR with cfg.keeptrials = ''yes'': ' ...
         'a permutation test needs the trials it permutes.']);

    if strcmp(design, 'baseline')
        assert(isempty(B), 'lfpClusterStats:baselineB', ...
            'the ''baseline'' design compares A with itself; pass B = [].');
        B = localBaselineClone(A, opt);
    end
    assert(isstruct(B) && isfield(B, 'powspctrm'), 'lfpClusterStats:inputB', ...
        'B must be a FieldTrip freq structure (or [] with design ''baseline'').');
    localAssertMatched(A, B);

    nA = size(A.powspctrm, 1);
    nB = size(B.powspctrm, 1);
    if ismember(design, {'paired', 'baseline'})
        assert(nA == nB, 'lfpClusterStats:unpaired', ...
            ['a paired design needs the same number of observations: A has %d, ' ...
             'B has %d. Use design ''trials'' for unequal, unpaired data.'], nA, nB);
        cfgDesign = [1:nA, 1:nB; ones(1, nA), 2 * ones(1, nB)];
        cfg.statistic = 'ft_statfun_depsamplesT';
        cfg.uvar = 1; cfg.ivar = 2;
    else
        cfgDesign = [ones(1, nA), 2 * ones(1, nB)];
        cfg.statistic = 'ft_statfun_indepsamplesT';
        cfg.ivar = 1;
    end

    cfg.method            = 'montecarlo';
    cfg.correctm          = 'cluster';
    cfg.clusteralpha      = a.clusteralpha;
    cfg.clusterstatistic  = a.clusterstatistic;
    cfg.minnbchan         = a.minnbchan;
    cfg.tail              = a.tail;
    cfg.clustertail       = a.tail;
    cfg.alpha             = a.alpha;
    cfg.correcttail       = 'prob';    % two-sided alpha, not alpha/2 by hand
    cfg.numrandomization  = a.numrand;
    cfg.design            = cfgDesign;
    cfg.latency           = a.latency;
    cfg.frequency         = a.frequency;
    cfg.avgoverchan       = localYesNo(a.avgoverchan);
    cfg.neighbours        = a.neighbours;
    if isempty(a.neighbours) && ~a.avgoverchan && numel(A.label) > 1
        warning('lfpClusterStats:noNeighbours', ...
            ['%d channels, avgoverchan off and no neighbours: clusters cannot ' ...
             'grow across channels, so each is tested in time-frequency alone.'], ...
            numel(A.label));
    end

    stat = ft_freqstatistics(cfg, A, B);

    info = struct('design', design, 'statistic', cfg.statistic, ...
                  'numrand', a.numrand, 'alpha', a.alpha, ...
                  'clusteralpha', a.clusteralpha, 'tail', a.tail, ...
                  'avgoverchan', logical(a.avgoverchan), ...
                  'latency', a.latency, 'frequency', a.frequency, ...
                  'nA', nA, 'nB', nB);
    [info.nPos, info.pPos] = localClusters(stat, 'posclusters', a.alpha);
    [info.nNeg, info.pNeg] = localClusters(stat, 'negclusters', a.alpha);
    info.anySignificant = info.nPos + info.nNeg > 0;
    % The smallest p a permutation test can return is 1/numrand: at the floor
    % the true p is only known to be below it.
    info.pFloor = 1 / a.numrand;
end

% ---------------- helpers ----------------
function B = localBaselineClone(A, opt)
% Each trial's baseline mean, spread over the full time axis: the "no change
% from baseline" null, with the same shape as A so FieldTrip can pair them.
    win = localOpt(opt, {'lfp','norm','baseline'}, ...
          localOpt(opt, {'lfp','plot','baseline'}, [-0.5 0]));
    idx = A.time >= win(1) & A.time <= win(2);
    assert(any(idx), 'lfpClusterStats:baselineEmpty', ...
        'baseline [%g %g] s contains no samples of the time axis (%g to %g s).', ...
        win(1), win(2), A.time(1), A.time(end));
    B = A;
    mu = mean(A.powspctrm(:, :, :, idx), 4, 'omitnan');
    B.powspctrm = repmat(mu, [1 1 1 numel(A.time)]);
end

function localAssertMatched(A, B)
    assert(isequal(numel(A.freq), numel(B.freq)) && isequal(numel(A.time), numel(B.time)), ...
        'lfpClusterStats:axes', ...
        ['A and B must share the frequency and time axes (A: %d freqs x %d times, ' ...
         'B: %d x %d). Compute both with the same cfg, or ft_selectdata them first.'], ...
        numel(A.freq), numel(A.time), numel(B.freq), numel(B.time));
    assert(isequal(sort(A.label), sort(B.label)), 'lfpClusterStats:labels', ...
        'A and B must contain the same channels.');
end

function [n, pvals] = localClusters(stat, field, alpha)
    n = 0; pvals = [];
    if isfield(stat, field) && ~isempty(stat.(field))
        pvals = [stat.(field).prob];
        n = sum(pvals < alpha);
    end
end

function s = localYesNo(tf)
    if ischar(tf), s = tf; elseif tf, s = 'yes'; else, s = 'no'; end
end

function d = localDimord(freq)
    if isfield(freq, 'dimord') && ~isempty(freq.dimord)
        d = freq.dimord;
    elseif ndims(freq.powspctrm) == 4
        d = 'rpt_chan_freq_time';
    else
        d = 'chan_freq_time';
    end
end

function v = localOpt(opt, path, default)
    v = default;
    s = opt;
    for k = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{k}), return; end
        s = s.(path{k});
    end
    if ~isempty(s) || ischar(s), v = s; end
end
