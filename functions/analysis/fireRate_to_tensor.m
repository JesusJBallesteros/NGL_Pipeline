function tensor = fireRate_to_tensor(fireRate, alignIdx)
% fireRate_to_tensor  Stack fireRate.sps cell-per-cluster into a 3-D tensor
%                     for a single alignment.
%
% PURPOSE:
%   The lab convention from calculate_fireRate_general stores per-cluster
%   firing-rate matrices in a cell array indexed by [cluster, alignment]
%   (#26): fireRate.sps{c, a}  ->  [Ntrials x Nbins]
%   Population analyses need a single contiguous tensor
%       [Nclust x Nbins x Ntrials]
%   for ONE alignment at a time. This helper picks the requested
%   alignment column and repacks it into a tensor.
%
% USAGE:
%   tensor = fireRate_to_tensor(fireRate)            % alignIdx = 1
%   tensor = fireRate_to_tensor(fireRate, alignIdx)  % explicit pick
%
% INPUTS:
%   fireRate - struct from calculate_fireRate_general. Required field:
%                .sps   {Nclust x Nalign} of [Ntrials x Nbins] matrices.
%              All cells in the picked column must share the same
%              (Ntrials, Nbins) shape.
%   alignIdx - positive integer column index into fireRate.sps (default 1).
%              Must be <= size(fireRate.sps, 2).
%
% OUTPUT:
%   tensor   - [Nclust x Nbins x Ntrials] double.
%
% NOTES:
%   - Clusters with shape mismatches are skipped with a warning, leaving
%     zeros in their slice. This shouldn't happen if all clusters come
%     from the same calculate_fireRate_general call.
%   - Pre-#26 fireRate.sps was {Nclust x 1}; with alignIdx defaulting to
%     1 this function still works on those caches transparently.
%
% SEE ALSO:
%   calculate_fireRate_general, calculate_neural_pca, smooth_spikes
%
% Last modified 02.06.2026 (Jesus) - alignIdx selector added (#26)

if nargin < 2 || isempty(alignIdx), alignIdx = 1; end

assert(isstruct(fireRate) && isfield(fireRate,'sps') && iscell(fireRate.sps), ...
    'NGL:fireRate_to_tensor', 'fireRate.sps must be a cell array.');

[Nclust, Nalign] = size(fireRate.sps);
assert(Nclust > 0, 'NGL:fireRate_to_tensor', 'fireRate.sps is empty.');
assert(isnumeric(alignIdx) && isscalar(alignIdx) && alignIdx >= 1 && ...
       alignIdx == floor(alignIdx) && alignIdx <= Nalign, ...
    'NGL:fireRate_to_tensor', ...
    'alignIdx must be a positive integer <= size(fireRate.sps, 2) = %d.', Nalign);

%% Determine the canonical (Ntrials, Nbins) from the first non-empty cell
%  in the requested alignment column.
ref = [];
for c = 1:Nclust
    cell_c = fireRate.sps{c, alignIdx};
    if ~isempty(cell_c) && isnumeric(cell_c)
        ref = size(cell_c);
        break
    end
end
assert(~isempty(ref), 'NGL:fireRate_to_tensor', ...
    'fireRate.sps(:, %d) contains no numeric matrices.', alignIdx);
Ntrials = ref(1);
Nbins   = ref(2);

tensor = zeros(Nclust, Nbins, Ntrials);
for c = 1:Nclust
    sps = fireRate.sps{c, alignIdx};
    if isempty(sps),                 continue; end
    if ~isnumeric(sps),              warning('NGL:fireRate_to_tensor', ...
        'fireRate.sps{%d, %d} is not numeric; skipping.', c, alignIdx); continue; end
    if ~isequal(size(sps), ref),     warning('NGL:fireRate_to_tensor', ...
        'fireRate.sps{%d, %d} shape %s mismatches reference %s; skipping.', ...
        c, alignIdx, mat2str(size(sps)), mat2str(ref)); continue; end
    tensor(c, :, :) = sps';  % [Nbins x Ntrials] -> placed into [1 x Nbins x Ntrials]
end
end
