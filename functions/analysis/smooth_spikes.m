function smoothed = smooth_spikes(rateTensor, sigma_s, binSize_s)
% smooth_spikes  Gaussian-kernel smoothing of a binned rate tensor.
%
% PURPOSE:
%   Convolve along the time axis of a [Nclust x Nbins x Ntrials] firing-rate
%   tensor with a Gaussian kernel. This is the standard pre-processing
%   step before any population dimensionality-reduction analysis (PCA,
%   jPCA, GPFA two-stage variants). Operates per-cluster, per-trial.
%
% USAGE:
%   smoothed = smooth_spikes(rateTensor, sigma_s, binSize_s)
%
% INPUTS:
%   rateTensor - [Nclust x Nbins x Ntrials] numeric. Typically built by
%                fireRate_to_tensor from fireRate.sps.
%   sigma_s    - Gaussian kernel SD in SECONDS.
%   binSize_s  - bin width of rateTensor in SECONDS.
%
% OUTPUT:
%   smoothed   - same shape as rateTensor, smoothed along axis 2 (time).
%
% NOTES:
%   - sigma_s = 0 returns the input unchanged.
%   - Edges are handled with 'same'-length zero padding (MATLAB conv default).
%     For most binSize values used here (~20 ms with ~50 ms sigma) the
%     first/last few bins are slightly attenuated; acceptable for
%     exploratory visualisation. If you need a less biased estimate at the
%     edges, use a method-of-images reflection or build kernel coverage
%     yourself.
%   - Reference: Park & Brockmeier (2013) tutorial, arxiv 1302.5964;
%     Yu et al. 2009, J. Neurophysiol.
%
% SEE ALSO:
%   fireRate_to_tensor, calculate_neural_pca
%
% Last modified 29.05.2026 (Jesus)

assert(isnumeric(rateTensor) && ndims(rateTensor) <= 3, ...
    'NGL:smooth_spikes', 'rateTensor must be a numeric 2-D or 3-D array.');
assert(isnumeric(sigma_s) && isscalar(sigma_s) && sigma_s >= 0, ...
    'NGL:smooth_spikes', 'sigma_s must be a non-negative scalar (seconds).');
assert(isnumeric(binSize_s) && isscalar(binSize_s) && binSize_s > 0, ...
    'NGL:smooth_spikes', 'binSize_s must be a positive scalar (seconds).');

if sigma_s == 0
    smoothed = rateTensor;
    return
end

%% Build Gaussian kernel in bin units, normalised to unit sum.
sigma_bins  = sigma_s / binSize_s;
halfWindow  = max(1, ceil(3 * sigma_bins));      % +/- 3 sigma
kernelIdx   = -halfWindow : halfWindow;
kernel      = exp(-(kernelIdx.^2) / (2 * sigma_bins^2));
kernel      = kernel / sum(kernel);

%% Convolve along time axis for every (cluster, trial) slice.
% Vectorised via reshape: stack into [Nbins x (Nclust*Ntrials)] columns.
sz = size(rateTensor);
if numel(sz) == 2,  sz(3) = 1; end
T  = sz(2);

stack    = permute(rateTensor, [2 1 3]);                 % [Nbins x Nclust x Ntrials]
stack    = reshape(stack, T, sz(1) * sz(3));             % [Nbins x (Nclust*Ntrials)]
smStack  = conv2(stack, kernel(:), 'same');              % column-wise convolution
smStack  = reshape(smStack, T, sz(1), sz(3));            % back to [Nbins x Nclust x Ntrials]
smoothed = permute(smStack, [2 1 3]);                    % back to [Nclust x Nbins x Ntrials]
end
