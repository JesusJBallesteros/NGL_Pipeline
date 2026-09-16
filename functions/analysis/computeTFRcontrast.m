function res = computeTFRcontrast(freq, spec, opt, varargin)
% computeTFRcontrast  Event-centered power, two conditions compared and tested.
%
% PURPOSE:
%   The core of the event-centered power question: does power differ between
%   two sets of trials, where, and is the difference more than the trial-to-
%   trial variability would give by chance? Normalisation, the test and the
%   maps come from the Phase 0 helpers, so this answer is expressed the same
%   way as every other analysis in the pipeline.
%
% USAGE:
%   res = computeTFRcontrast(freq, spec, opt)
%   res = computeTFRcontrast(freq, spec, opt, 'area', 'NCL', 'align', 'stim2')
%
% INPUTS:
%   freq - FieldTrip freq structure for ONE band and ONE alignment, trials
%          kept ('rpt_chan_freq_time'), channels already restricted to the
%          area of interest. Raw power: this function normalises.
%   spec - contrast from parseTrialContrast (.A / .B masks and labels).
%   opt  - options; reads opt.lfp.norm.* and opt.lfp.stats.* (Phase 0).
%   Name/value pairs:
%     'area' / 'align'  recorded in the result for naming and titles
%     'stats'           false to skip the permutation test (plots only)
%     'statArgs'        cell of name/value pairs passed to lfpClusterStats
%
% OUTPUT (struct):
%   .meanA/.meanB  [nFreq x nTime] normalised power, averaged over trials
%                  then channels
%   .diff          meanA - meanB
%   .time/.freq    axes
%   .stat          ft_freqstatistics output ([] when 'stats' false)
%   .mask          [nFreq x nTime] logical, the cluster mask on the full time
%                  axis (false outside the tested window)
%   .norm/.stats   info structs for provenance
%   .spec/.area/.align/.nA/.nB
%
% NOTES:
%   * The test runs on RAW power, the maps on normalised power. Baseline
%     normalisation rescales each trial by its own baseline, which changes
%     the variance the test weighs; the contrast itself does not need it,
%     since both sides share the baseline treatment.
%   * Averaging order is trials first, then channels: the per-trial dB values
%     are what the average should be taken over (see normalizeTFR).
%   * The mask is returned on the full time axis even when the test ran on a
%     shorter window, so it can be drawn over the whole map without the
%     caller re-deriving the indices.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 1).

    p = inputParser;
    p.addParameter('area',     '');
    p.addParameter('align',    '');
    p.addParameter('stats',    true);
    p.addParameter('statArgs', {});
    p.parse(varargin{:});
    a = p.Results;

    assert(isstruct(freq) && isfield(freq, 'powspctrm'), 'computeTFRcontrast:input', ...
        'freq must be a FieldTrip freq structure.');
    nTrials = size(freq.powspctrm, 1);
    assert(numel(spec.A) == nTrials, 'computeTFRcontrast:masks', ...
        ['the contrast masks cover %d trials but this TFR holds %d. Parse the ', ...
         'contrast against the same TFR you are testing.'], numel(spec.A), nTrials);

    A = ft_selectdata(struct('trials', find(spec.A)), freq);
    B = ft_selectdata(struct('trials', find(spec.B)), freq);

    [nA, normInfo] = normalizeTFR(A, opt);
    nB = normalizeTFR(B, opt);

    res = struct();
    res.meanA = localReduce(nA.powspctrm);
    res.meanB = localReduce(nB.powspctrm);
    res.diff  = res.meanA - res.meanB;
    res.time  = freq.time;
    res.freq  = freq.freq;
    res.norm  = normInfo;
    res.spec  = spec;
    res.area  = char(a.area);
    res.align = char(a.align);
    res.nA    = sum(spec.A);
    res.nB    = sum(spec.B);
    res.stat  = [];
    res.stats = struct();
    res.mask  = false(numel(freq.freq), numel(freq.time));

    if ~a.stats
        return
    end
    [res.stat, res.stats] = lfpClusterStats(A, B, opt, a.statArgs{:});
    res.mask = localFullMask(res.stat, freq);
end

% ---------------- helpers ----------------
function M = localReduce(P)
% [rpt x chan x freq x time] -> [freq x time], trials first then channels.
    M = squeeze(mean(mean(P, 1, 'omitnan'), 2, 'omitnan'));
end

function mask = localFullMask(stat, freq)
    mask = false(numel(freq.freq), numel(freq.time));
    if ~isfield(stat, 'mask') || isempty(stat.mask), return; end
    m = squeeze(logical(stat.mask));
    [~, fIdx] = ismember(round(stat.freq(:), 6), round(freq.freq(:), 6));
    [~, tIdx] = ismember(round(stat.time(:), 6), round(freq.time(:), 6));
    ok = fIdx > 0; okT = tIdx > 0;
    if ~all(ok) || ~all(okT)
        % ft_freqstatistics can return axes that do not land exactly on the
        % input grid (averaging, rounding). Fall back to nearest bins rather
        % than dropping the mask.
        fIdx = interp1(freq.freq, 1:numel(freq.freq), stat.freq, 'nearest', 'extrap');
        tIdx = interp1(freq.time, 1:numel(freq.time), stat.time, 'nearest', 'extrap');
    end
    mask(fIdx, tIdx) = m;
end
