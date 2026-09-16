function [freq, info] = normalizeTFR(freq, opt, varargin)
% normalizeTFR  Baseline-normalise a FieldTrip TFR, one way for the whole pipeline.
%
% PURPOSE:
%   Every power figure and every power statistic needs the same question
%   answered first: power relative to what? Left to each analysis, the answer
%   drifts (one plot in dB, one in % change, one raw), and the results stop
%   being comparable. This is the single place that answer is given, and the
%   choice it made travels with the data in `info` so provenance records it.
%
%   ft_freqbaseline does the same job, but only for 'absolute' / 'relative' /
%   'relchange' / 'db' on a non-rpt structure, and it silently averages
%   trials away. This keeps single trials (statistics need them), adds the
%   z-score against baseline variability, and refuses quietly-wrong input
%   instead of returning Inf.
%
% USAGE:
%   [freq, info] = normalizeTFR(freq, opt)
%   [freq, info] = normalizeTFR(freq, opt, 'method', 'db', 'baseline', [-0.4 0])
%
% INPUTS:
%   freq - FieldTrip freq structure from ft_freqanalysis, dimord
%          'rpt_chan_freq_time' or 'chan_freq_time'. Needs .powspctrm/.time.
%   opt  - options struct. Reads opt.lfp.norm.method / .baseline /
%          .singleTrial, falling back to opt.lfp.plot.baseline for the window
%          so existing configs keep working.
%   Name/value pairs override the opt values:
%     'method'      'none' | 'absolute' | 'relchange' | 'percent' | 'db' | 'z'
%     'baseline'    [tmin tmax] in seconds, on freq.time
%     'singleTrial' logical; true = each trial against its own baseline
%
% OUTPUT:
%   freq - same structure, .powspctrm normalised. Units change with method:
%          absolute = power difference, relchange = fraction, percent = %,
%          db = 10*log10(ratio), z = SD of the baseline.
%   info - what was done: .method .baseline .nSamples .singleTrial .units
%          Put it in the provenance of whatever you save.
%
% NOTES:
%   * 'db', 'relchange' and 'percent' divide by baseline power, so a channel
%     whose baseline is zero or negative yields NaN rather than Inf, with one
%     warning naming the channels. Negative power only happens on already
%     normalised input - normalising twice is the usual cause.
%   * 'z' needs variability, so it is computed across the baseline's time
%     samples (and trials, when not single-trial). With fewer than 2 baseline
%     samples it errors rather than dividing by zero.
%   * Trials are kept. Average after normalising, not before: with 'db' the
%     two orders differ, and the per-trial one is what statistics need.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 0).

    p = inputParser;
    p.addParameter('method',      localOpt(opt, {'lfp','norm','method'}, 'db'));
    p.addParameter('baseline',    localOpt(opt, {'lfp','norm','baseline'}, ...
                                  localOpt(opt, {'lfp','plot','baseline'}, [-0.5 0])));
    p.addParameter('singleTrial', localOpt(opt, {'lfp','norm','singleTrial'}, true));
    p.parse(varargin{:});
    method      = lower(char(p.Results.method));
    baseline    = p.Results.baseline;
    singleTrial = logical(p.Results.singleTrial);

    assert(isstruct(freq) && isfield(freq, 'powspctrm') && isfield(freq, 'time'), ...
        'normalizeTFR:input', 'freq must be a FieldTrip freq structure with powspctrm and time.');
    valid = {'none', 'absolute', 'relchange', 'percent', 'db', 'z'};
    assert(ismember(method, valid), 'normalizeTFR:method', ...
        'method must be one of: %s (got ''%s'').', strjoin(valid, ', '), method);

    info = struct('method', method, 'baseline', baseline, 'nSamples', 0, ...
                  'singleTrial', singleTrial, 'units', localUnits(method));
    if strcmp(method, 'none')
        return
    end

    assert(isnumeric(baseline) && numel(baseline) == 2 && baseline(1) < baseline(2), ...
        'normalizeTFR:baseline', 'baseline must be [tmin tmax] with tmin < tmax.');
    bIdx = freq.time >= baseline(1) & freq.time <= baseline(2);
    info.nSamples = sum(bIdx);
    assert(info.nSamples > 0, 'normalizeTFR:baselineEmpty', ...
        ['baseline [%g %g] s contains no samples of freq.time (%g to %g s). ' ...
         'Widen it, or check that this alignment has pre-event time.'], ...
        baseline(1), baseline(2), freq.time(1), freq.time(end));
    if strcmp(method, 'z')
        assert(info.nSamples > 1, 'normalizeTFR:baselineTooShort', ...
            'z-scoring needs at least 2 baseline samples; [%g %g] s gives 1.', ...
            baseline(1), baseline(2));
    end

    % Work in [rpt x chan x freq x time] so one code path covers both dimords.
    hasRpt = contains(lower(localDimord(freq)), 'rpt');
    P = freq.powspctrm;
    if ~hasRpt
        P = reshape(P, [1 size(P, 1) size(P, 2) size(P, 3)]);
        singleTrial = true;      % one "trial": the two paths coincide
    end

    base = P(:, :, :, bIdx);
    mu = mean(base, 4, 'omitnan');
    if ~singleTrial
        mu = repmat(mean(mu, 1, 'omitnan'), [size(P, 1) 1 1 1]);
    end

    if strcmp(method, 'z')
        if singleTrial
            sd = std(base, 0, 4, 'omitnan');
        else
            flat = reshape(permute(base, [1 4 2 3]), [], size(P, 2), size(P, 3));
            sd = repmat(reshape(std(flat, 0, 1, 'omitnan'), ...
                                [1 size(P, 2) size(P, 3)]), [size(P, 1) 1 1]);
        end
        sd(sd == 0) = NaN;                    % flat baseline: undefined, not Inf
        P = (P - mu) ./ sd;
    elseif strcmp(method, 'absolute')
        P = P - mu;
    else
        bad = mu <= 0 | ~isfinite(mu);
        if any(bad(:))
            chans = localBadChannels(freq, bad);
            warning('normalizeTFR:nonPositiveBaseline', ...
                ['%d of %d baseline values are zero, negative or non-finite, so ' ...
                 '''%s'' is undefined there and returns NaN (channels: %s). ' ...
                 'Normalising an already-normalised TFR is the usual cause.'], ...
                sum(bad(:)), numel(bad), method, chans);
            mu(bad) = NaN;
        end
        switch method
            case 'db',        P = 10 * log10(P ./ mu);
            case 'relchange', P = (P - mu) ./ mu;
            case 'percent',   P = 100 * (P - mu) ./ mu;
        end
    end

    if ~hasRpt
        P = reshape(P, size(P, 2), size(P, 3), size(P, 4));
    end
    freq.powspctrm = P;
    freq.normalize = info;      % survives ft_selectdata; also copied to provenance
end

% ---------------- helpers ----------------
function d = localDimord(freq)
    if isfield(freq, 'dimord'), d = freq.dimord; else, d = ''; end
    if isempty(d)
        % No dimord: infer from shape. 4-D can only be rpt_chan_freq_time.
        if ndims(freq.powspctrm) == 4, d = 'rpt_chan_freq_time';
        else,                          d = 'chan_freq_time';
        end
    end
end

function s = localUnits(method)
    switch method
        case 'db',        s = 'dB re baseline';
        case 'relchange', s = 'change (fraction of baseline)';
        case 'percent',   s = '% change from baseline';
        case 'z',         s = 'z (SD of baseline)';
        case 'absolute',  s = 'power - baseline';
        otherwise,        s = 'power (a.u.)';
    end
end

function txt = localBadChannels(freq, bad)
    idx = find(any(any(any(bad, 1), 3), 4));
    if isfield(freq, 'label') && ~isempty(freq.label)
        names = freq.label(idx(idx <= numel(freq.label)));
    else
        names = arrayfun(@(k) sprintf('#%d', k), idx, 'UniformOutput', false);
    end
    if numel(names) > 6, names = [names(1:6); {'...'}]; end
    txt = strjoin(reshape(names, 1, []), ', ');
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
