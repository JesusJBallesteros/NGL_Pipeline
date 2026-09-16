function spec = computeTaggingSpectrum(epochs, opt, varargin)
% computeTaggingSpectrum  Amplitude spectrum of tagging epochs, with phase kept.
%
% PURPOSE:
%   The frequency-domain view a tagging experiment is read from: amplitude per
%   bin, averaged over epochs, plus the complex spectrum so phase consistency
%   across epochs can be measured from the same transform.
%
% USAGE:
%   spec = computeTaggingSpectrum(epochs, opt)
%   spec = computeTaggingSpectrum(epochs, opt, 'fmax', 30, 'taper', 'boxcar')
%
% INPUTS:
%   epochs - FieldTrip raw structure from nftEpochs (trials = epochs, all the
%            same length - unequal epochs would put the bins in different
%            places and the average would be meaningless).
%   opt    - options; reads opt.nft.taper / .fmax.
%   Name/value pairs: 'taper', 'fmax', 'pad'.
%
% OUTPUT (struct):
%   .freq      [1 x nFreq] bin centres, Hz
%   .amp       [nChan x nFreq] amplitude, averaged over epochs
%   .fourier   [nEpoch x nChan x nFreq] complex, for computeITPC
%   .label     channel labels
%   .nEpochs, .taper, .resolution
%
% NOTES:
%   * Amplitude, not power: the tagging literature reports amplitude (and
%     baseline-subtracted amplitude sums over harmonics), and amplitude keeps
%     the units of the signal. Square it for power if a comparison needs it.
%   * Averaging amplitude over epochs is deliberate: it keeps a response that
%     is present in every epoch but drifting in phase, which averaging the
%     complex spectra would cancel. Phase consistency is measured separately,
%     by computeITPC, where cancelling is the point.
%   * Taper: 'hanning' (default) is safe. With whole-cycle epochs 'boxcar'
%     gives the sharpest peak - no leakage to taper away - but any drift in
%     the stimulation rate then leaks badly. Prefer boxcar only when the
%     hardware timing is trusted.
%
% Last modified 16.09.2026 (Jesus) - new (NFT project).

    p = inputParser;
    p.addParameter('taper', localOpt(opt, {'nft','taper'}, 'hanning'));
    p.addParameter('fmax',  localOpt(opt, {'nft','fmax'},  40));
    p.addParameter('pad',   'maxperlen');
    p.parse(varargin{:});
    a = p.Results;

    assert(isstruct(epochs) && isfield(epochs, 'trial') && ~isempty(epochs.trial), ...
        'computeTaggingSpectrum:input', 'epochs must be a FieldTrip raw structure.');
    lens = cellfun(@(x) size(x, 2), epochs.trial);
    assert(isscalar(unique(lens)), 'computeTaggingSpectrum:unequal', ...
        ['epochs differ in length (%d to %d samples). Equal epochs are what puts ', ...
         'every FFT bin in the same place; use nftEpochs to cut them.'], ...
        min(lens), max(lens));

    cfg = struct('method', 'mtmfft', 'output', 'fourier', 'taper', a.taper, ...
                 'foilim', [0 a.fmax], 'pad', a.pad, 'keeptrials', 'yes');
    F = ft_freqanalysis(cfg, epochs);

    % fourierspctrm is [rpt x chan x freq]; amplitude is its modulus. The
    % factor 2/N would convert to physical amplitude, but every measure here
    % (SNR, z, ITPC) is a ratio within the same spectrum, so the scaling
    % cancels and is left alone rather than half-applied.
    spec = struct();
    spec.freq       = F.freq;
    spec.fourier    = F.fourierspctrm;
    spec.amp        = squeeze(mean(abs(F.fourierspctrm), 1, 'omitnan'));
    if isvector(spec.amp), spec.amp = reshape(spec.amp, 1, []); end
    spec.label      = F.label;
    spec.nEpochs    = size(F.fourierspctrm, 1);
    spec.taper      = char(a.taper);
    spec.resolution = median(diff(F.freq));
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
