function plot_fireRate_session(fireRate, condition, param, opt)
% plot_fireRate_session  Per-cluster + session-level firing-rate plots,
%                        with trial2plot selection applied here (#19).
%
% PURPOSE:
%   Owns the plotting that used to live inside calculate_fireRate_general
%   (audit item S). Since #19, calculate_fireRate_general no longer drops
%   rows for param.trial2plot — it returns the FULL trial axis. This
%   plotting helper is now where the trial2plot selection is applied,
%   via applyTrialFilter, before the per-cluster heatmaps and session-
%   level multi-cluster plot are drawn.
%
% USAGE:
%   plot_fireRate_session(fireRate, condition, param, opt)
%   Called by NGL02_postPhy after calculate_fireRate_general (per area
%   in multi-area mode).
%
% INPUTS:
%   fireRate  - struct from calculate_fireRate_general. Required:
%                 .sps      {Nclust x 1} of [Ntotal x Nbins] (raw)
%                 .Norm     {Nclust x 1} of [Ntotal x Nbins] (normalised)
%                 .meanNorm {Nclust x 1} of [1 x Nbins] (mean over ALL
%                           trials; this function recomputes the filtered
%                           mean for the multi-cluster plot).
%   condition - per-trial condition struct. Used by applyTrialFilter to
%               build the trial2plot mask. Pass an empty struct() to skip
%               filtering and draw every row.
%   param     - param struct. Used: .plot (gate), .trial2plot (filter rule),
%               .cl (set per-cluster), .ROI/.KSLabel/etc. (subtitle).
%   opt       - resolved options struct; forwarded to the plot helpers
%               (alignto, area, analysis, SavFileName, etc.).
%
% PARAM DEFAULTS:
%   .trial2plot 'allInitiated'  used by applyTrialFilter.
%   .plot       true            gate; false makes this function a no-op.
%
% SEE ALSO:
%   calculate_fireRate_general, applyTrialFilter, plot_single_fireRate,
%   plot_multi_fireRate.
%
% Last modified 02.06.2026 (Jesus) - trial2plot selection moved here (#19)

%% No-op gate.
if isfield(param,'plot') && ~param.plot, return; end
if ~isfield(param,'trial2plot') || isempty(param.trial2plot)
    param.trial2plot = 'allInitiated';
end

%% Build the row mask from condition. If condition lacks the needed
%  fields, fall back to "no filtering" (all rows).
try
    rowMask = applyTrialFilter(condition, param.trial2plot);
catch ME
    warning('NGL:plot_fireRate_session:filterFallback', ...
        'applyTrialFilter failed (%s); drawing all rows.', ME.message);
    if ~isempty(fireRate.sps) && ~isempty(fireRate.sps{1})
        rowMask = true(1, size(fireRate.sps{1}, 1));
    else
        rowMask = [];
    end
end

%% Per-cluster plots.
Nclust = numel(fireRate.sps);
for c = 1:Nclust
    param.cl = [1 c 1];  % alignment 1, cluster c, level 1
    sps  = fireRate.sps{c,1};
    nrm  = fireRate.Norm{c,1};
    if ~isempty(rowMask) && size(sps,1) == numel(rowMask)
        sps = sps(rowMask, :);
        nrm = nrm(rowMask, :);
    end
    plot_single_fireRate(sps, nrm, param, opt);
end

%% Session-level multi-cluster plot.
% Recompute meanNorm from the filtered .Norm so the multi-cluster heatmap
% reflects the same trial selection as the per-cluster plots.
if ~isempty(rowMask)
    filteredMeanNorm = cell(Nclust, 1);
    for c = 1:Nclust
        nrm = fireRate.Norm{c,1};
        if size(nrm,1) == numel(rowMask)
            filteredMeanNorm{c,1} = mean(nrm(rowMask,:), 1, 'omitnan');
        else
            filteredMeanNorm{c,1} = fireRate.meanNorm{c,1};
        end
    end
    plot_multi_fireRate(filteredMeanNorm, param, opt);
else
    plot_multi_fireRate(fireRate.meanNorm, param, opt);
end
end
