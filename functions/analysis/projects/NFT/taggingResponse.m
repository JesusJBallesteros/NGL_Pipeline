function resp = taggingResponse(spec, baseFreq, opt, varargin)
% taggingResponse  Is there a response at the tagging frequency, and how big?
%
% PURPOSE:
%   A tagged response is a peak at a frequency fixed by the experiment, so it
%   is tested against its own neighbourhood rather than against a baseline
%   period: the neighbouring bins hold the same noise at almost the same
%   frequency, measured at the same time. That is the standard frequency-
%   tagging logic (Retter & Rossion 2016, and the SSVEP literature before it).
%
%   Two numbers per frequency: SNR (amplitude / mean neighbour), which is 1
%   when nothing is there, and z ((amplitude - mean) / SD of neighbours),
%   which says how unlikely that is. A response usually spreads over the
%   harmonics of the base, so the summed baseline-corrected amplitude over
%   significant harmonics is the response's size.
%
% USAGE:
%   resp = taggingResponse(spec, 1.3, opt)
%   resp = taggingResponse(spec, 1.3, opt, 'exclude', [2.6 5.2])
%
% INPUTS:
%   spec     - from computeTaggingSpectrum.
%   baseFreq - tagging frequency of this block, Hz.
%   opt      - reads opt.nft.neighbours / .gap / .zThreshold / .maxHarmonic /
%              .lineFreq.
%   Name/value pairs:
%     'neighbours'  bins per side forming the baseline (default 12)
%     'gap'         bins next to the peak excluded from the baseline
%                   (default 1) - a tapered peak bleeds into its neighbours,
%                   and counting that bleed as baseline hides the peak
%     'zThreshold'  z above which a harmonic counts (default 1.64, one-sided
%                   p = 0.05; with many harmonics tested, treat it as a
%                   screen rather than a corrected test)
%     'maxHarmonic' highest harmonic considered (default 8)
%     'exclude'     frequencies whose bins must not enter any baseline or
%                   harmonic set (e.g. another block's tag)
%     'lineFreq'    mains frequency to exclude, Hz (default 50; 0 = ignore)
%
% OUTPUT (struct):
%   .harmonics   [1 x nH] frequencies actually used (n x base, minus excluded)
%   .snr/.z      [nChan x nH] at each harmonic
%   .amp         [nChan x nH] raw amplitude
%   .corrected   [nChan x nH] amplitude minus its neighbourhood mean
%   .sumAmp      [nChan x 1] summed corrected amplitude over SIGNIFICANT
%                harmonics (zero when none reach threshold)
%   .nSignificant[nChan x 1]
%   .snrSpectrum/.zSpectrum  [nChan x nFreq] for plotting
%   .info        the settings used
%
% NOTES:
%   * Summing only significant harmonics is a choice with a bias: it cannot
%     go below zero, so a channel with no response sums to exactly zero rather
%     than scattering around it. Compare channels by `sumAmp` only when both
%     have at least one significant harmonic; otherwise compare `z`.
%   * A harmonic landing within one bin of an excluded frequency is dropped,
%     not shifted; when blocks share harmonics (2.6 Hz is the 2nd harmonic of
%     1.3 Hz) that overlap is a property of the design, and hiding it by
%     nudging the bin would be worse than losing the harmonic.
%
% Last modified 16.09.2026 (Jesus) - new (NFT project).

    p = inputParser;
    p.addParameter('neighbours',  localOpt(opt, {'nft','neighbours'}, 12));
    p.addParameter('gap',         localOpt(opt, {'nft','gap'}, 1));
    p.addParameter('zThreshold',  localOpt(opt, {'nft','zThreshold'}, 1.64));
    p.addParameter('maxHarmonic', localOpt(opt, {'nft','maxHarmonic'}, 8));
    p.addParameter('lineFreq',    localOpt(opt, {'nft','lineFreq'}, 50));
    p.addParameter('exclude',     []);
    p.parse(varargin{:});
    a = p.Results;

    f = spec.freq(:)';
    A = spec.amp;
    if isvector(A), A = reshape(A, 1, []); end
    nChan = size(A, 1);
    df = median(diff(f));
    assert(df > 0, 'taggingResponse:freq', 'the spectrum has no usable frequency axis.');

    % --- SNR / z for every bin, so the whole spectrum can be plotted --------
    [snrSpec, zSpec] = localNeighbourStats(A, a.neighbours, a.gap);

    % --- which harmonics to use -------------------------------------------
    cand = baseFreq * (1:a.maxHarmonic);
    cand = cand(cand <= max(f) - a.neighbours * df);   % need a neighbourhood
    drop = false(size(cand));
    excl = a.exclude(:)';
    if a.lineFreq > 0
        excl = [excl, a.lineFreq * (1:floor(max(f) / a.lineFreq))];
    end
    for k = 1:numel(cand)
        if any(abs(cand(k) - excl) <= df)
            drop(k) = true;
        end
    end
    harm = cand(~drop);
    assert(~isempty(harm), 'taggingResponse:noHarmonics', ...
        ['no harmonic of %g Hz survives: all fall on excluded frequencies or ', ...
         'outside the spectrum.'], baseFreq);

    idx = arrayfun(@(x) localNearestBin(f, x), harm);
    offBy = abs(f(idx) - harm);
    if any(offBy > df / 2 + eps)
        warning('taggingResponse:offBin', ...
            ['a harmonic misses its bin by up to %.4f Hz (bin width %.4f Hz). ', ...
             'Epochs holding whole cycles put it exactly on one - see nftEpochs.'], ...
            max(offBy), df);
    end

    resp = struct();
    resp.harmonics   = f(idx);
    resp.amp         = A(:, idx);
    resp.snr         = snrSpec(:, idx);
    resp.z           = zSpec(:, idx);
    resp.corrected   = A(:, idx) - localNeighbourMean(A, a.neighbours, a.gap, idx);
    sig              = resp.z >= a.zThreshold;
    resp.significant = sig;
    resp.corrected(~isfinite(resp.corrected)) = 0;
    resp.sumAmp      = sum(resp.corrected .* sig, 2);
    resp.nSignificant = sum(sig, 2);
    resp.snrSpectrum = snrSpec;
    resp.zSpectrum   = zSpec;
    resp.label       = spec.label;
    resp.freq        = f;
    resp.info = struct('baseFreq', baseFreq, 'neighbours', a.neighbours, ...
                       'gap', a.gap, 'zThreshold', a.zThreshold, ...
                       'maxHarmonic', a.maxHarmonic, 'lineFreq', a.lineFreq, ...
                       'excluded', excl, 'resolution', df, 'nChan', nChan);
end

% ---------------- helpers ----------------
function [snr, z] = localNeighbourStats(A, nb, gap)
    [nChan, nF] = size(A);
    snr = nan(nChan, nF); z = nan(nChan, nF);
    for k = 1:nF
        idx = localNeighbourIdx(k, nb, gap, nF);
        if numel(idx) < 4, continue; end
        base = A(:, idx);
        mu = mean(base, 2, 'omitnan');
        sd = std(base, 0, 2, 'omitnan');
        snr(:, k) = A(:, k) ./ mu;
        sd(sd == 0) = NaN;                   % flat neighbourhood: undefined
        z(:, k) = (A(:, k) - mu) ./ sd;
    end
end

function mu = localNeighbourMean(A, nb, gap, bins)
    mu = nan(size(A, 1), numel(bins));
    for k = 1:numel(bins)
        idx = localNeighbourIdx(bins(k), nb, gap, size(A, 2));
        if isempty(idx), continue; end
        mu(:, k) = mean(A(:, idx), 2, 'omitnan');
    end
end

function idx = localNeighbourIdx(k, nb, gap, nF)
    lo = max(1, k - nb - gap) : max(1, k - gap - 1);
    hi = min(nF, k + gap + 1) : min(nF, k + nb + gap);
    idx = unique([lo, hi]);
    idx = idx(idx >= 1 & idx <= nF & abs(idx - k) > gap);
end

function b = localNearestBin(f, x)
    [~, b] = min(abs(f - x));
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
