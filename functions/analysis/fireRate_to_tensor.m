function tensor = fireRate_to_tensor(fireRate)
% fireRate_to_tensor  Stack fireRate.sps cell-per-cluster into a 3-D tensor.
%
% PURPOSE:
%   The lab convention from calculate_fireRate_general stores per-cluster
%   firing-rate matrices in a cell array:
%       fireRate.sps{c}   ->  [Ntrials x Nbins]
%   Population analyses need a single contiguous tensor
%       [Nclust x Nbins x Ntrials]
%   This helper performs that repacking.
%
% USAGE:
%   tensor = fireRate_to_tensor(fireRate)
%
% INPUT:
%   fireRate - struct from calculate_fireRate_general. Required field:
%                .sps   {Nclust x 1} of [Ntrials x Nbins] matrices.
%              All cells must share the same (Ntrials, Nbins) shape.
%
% OUTPUT:
%   tensor   - [Nclust x Nbins x Ntrials] double.
%
% NOTES:
%   - Clusters with shape mismatches are skipped with a warning, leaving
%     zeros in their slice. This shouldn't happen if all clusters come
%     from the same calculate_fireRate_general call.
%
% SEE ALSO:
%   calculate_fireRate_general, calculate_neural_pca, smooth_spikes
%
% Last modified 29.05.2026 (Jesus)

assert(isstruct(fireRate) && isfield(fireRate,'sps') && iscell(fireRate.sps), ...
    'NGL:fireRate_to_tensor', 'fireRate.sps must be a cell array.');

Nclust = numel(fireRate.sps);
assert(Nclust > 0, 'NGL:fireRate_to_tensor', 'fireRate.sps is empty.');

%% Determine the canonical (Ntrials, Nbins) from the first non-empty cell.
ref = [];
for c = 1:Nclust
    if ~isempty(fireRate.sps{c}) && isnumeric(fireRate.sps{c})
        ref = size(fireRate.sps{c});
        break
    end
end
assert(~isempty(ref), 'NGL:fireRate_to_tensor', 'fireRate.sps contains no numeric matrices.');
Ntrials = ref(1);
Nbins   = ref(2);

tensor = zeros(Nclust, Nbins, Ntrials);
for c = 1:Nclust
    sps = fireRate.sps{c};
    if isempty(sps),                 continue; end
    if ~isnumeric(sps),              warning('NGL:fireRate_to_tensor', ...
        'fireRate.sps{%d} is not numeric; skipping.', c); continue; end
    if ~isequal(size(sps), ref),     warning('NGL:fireRate_to_tensor', ...
        'fireRate.sps{%d} shape %s mismatches reference %s; skipping.', ...
        c, mat2str(size(sps)), mat2str(ref)); continue; end
    tensor(c, :, :) = sps';  % [Nbins x Ntrials] -> placed into [1 x Nbins x Ntrials]
end
end
