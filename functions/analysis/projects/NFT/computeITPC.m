function itpc = computeITPC(spec, freqs, varargin)
% computeITPC  Phase consistency across epochs at the tagging frequencies.
%
% PURPOSE:
%   Amplitude says a rhythm is present; phase consistency says it is locked to
%   the stimulation. A response driven by the stimulus keeps the same phase in
%   every epoch, so the unit vectors add; ongoing activity at the same
%   frequency does not, so they cancel. That separation is what makes tagging
%   a test of driving rather than of power.
%
% USAGE:
%   itpc = computeITPC(spec, [1.3 2.6 3.9])
%   itpc = computeITPC(spec, resp.harmonics, 'alpha', 0.01)
%
% INPUTS:
%   spec  - from computeTaggingSpectrum (needs .fourier).
%   freqs - frequencies to evaluate, Hz; each snapped to its nearest bin.
%   Name/value pairs:
%     'alpha' significance level for the Rayleigh test (default 0.05)
%
% OUTPUT (struct):
%   .freq       [1 x nF] the bins actually used
%   .itpc       [nChan x nF] in [0 1]; 1 = identical phase every epoch
%   .rayleighZ  [nChan x nF] = n * itpc^2
%   .p          [nChan x nF] Rayleigh p, exp(-Z) with the standard
%               small-sample correction
%   .significant[nChan x nF] logical at `alpha`
%   .nEpochs    epochs the estimate rests on
%
% NOTES:
%   * ITPC is biased upward when epochs are few: with n epochs, even random
%     phases give an expected ITPC near sqrt(pi)/2/sqrt(n). Below ~10 epochs
%     the number is more bias than signal, and the function warns. Compare
%     ITPC only between conditions with the same epoch count, or use
%     rayleighZ, which accounts for n.
%   * Overlapping epochs (nftEpochs 'overlap' > 0) share samples, so their
%     phases are not independent and both the Rayleigh test and the bias
%     correction become optimistic. The warning fires on epoch count alone
%     and cannot see that - check how the epochs were cut.
%
% Last modified 16.09.2026 (Jesus) - new (NFT project).

    p = inputParser;
    p.addParameter('alpha', 0.05);
    p.parse(varargin{:});
    alpha = p.Results.alpha;

    assert(isstruct(spec) && isfield(spec, 'fourier') && ~isempty(spec.fourier), ...
        'computeITPC:input', ...
        ['spec must carry .fourier - compute the spectrum with ', ...
         'computeTaggingSpectrum, which keeps the complex values.']);
    F = spec.fourier;                        % [rpt x chan x freq]
    n = size(F, 1);
    assert(n > 1, 'computeITPC:oneEpoch', ...
        ['phase consistency needs more than one epoch; this spectrum has %d. ', ...
         'Cut the block with nftEpochs in ''split'' mode.'], n);
    if n < 10
        warning('computeITPC:fewEpochs', ...
            ['%d epochs: ITPC is biased upward at small n (random phases would ', ...
             'already give about %.2f). Prefer rayleighZ, which accounts for n.'], ...
            n, sqrt(pi) / 2 / sqrt(n));
    end

    bins = arrayfun(@(x) localNearestBin(spec.freq, x), freqs(:)');
    itpc = struct();
    itpc.freq = spec.freq(bins);
    itpc.nEpochs = n;

    Z = F(:, :, bins);                       % [rpt x chan x nF]
    unit = Z ./ abs(Z);
    unit(~isfinite(unit)) = NaN;             % a zero-amplitude bin has no phase
    R = squeeze(abs(mean(unit, 1, 'omitnan')));
    if size(F, 2) == 1, R = reshape(R, 1, []); end
    itpc.itpc = R;
    itpc.rayleighZ = n * R.^2;
    % Zar's small-sample form; exact enough from n ~ 10 and conservative below.
    itpc.p = exp(sqrt(1 + 4*n + 4*(n^2 - (n*R).^2)) - (1 + 2*n));
    itpc.significant = itpc.p < alpha;
    itpc.alpha = alpha;
    itpc.label = spec.label;
end

function b = localNearestBin(f, x)
    [~, b] = min(abs(f - x));
end
