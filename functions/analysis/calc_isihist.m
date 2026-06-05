function [isihist] = calc_isihist(spike, opt)
% calc_isihist  Per-cluster inter-spike interval histogram.
%
% PURPOSE:
%   For every cluster in the spike struct, computes ISI = diff(timestamps)
%   in ms and bins it into opt.isibins. Called by loadSpikes at the end
%   of its per-cluster loop as a quick quality-inspection helper.
%
% USAGE:
%   isihist = calc_isihist(spike, opt)
%
% INPUTS:
%   spike - NGL spike struct from loadSpikes. Required: .label, .timestamp
%           (cell-per-cluster, SECONDS).
%   opt   - resolved options struct. Used field: opt.isibins (defaulted
%           by set_default; vector of bin edges in ms).
%
% OUTPUT:
%   isihist - {1 x nClust} cell of histogram counts (one row vector per
%             cluster), aligned with opt.isibins edges.
%
% Last modified 02.06.2026 (Jesus) - docstring (#10 Q)

%% Get relevant info
nclust = numel(spike.label); % number of clusters
bins   = opt.isibins;        % in milliseconds (always present after set_default)

%% Run per cluster
for cl = 1:nclust
    % ISI is diff ts(n)-ts(n-1)
    isi = [NaN diff((spike.timestamp{cl}(:)*1000)')]'; % timestamps from seconds to miliseconds
    isihist{cl}   = histcounts(isi,bins);    
end

end