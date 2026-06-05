function plot_fireRate_session(fireRate, param, opt)
% plot_fireRate_session  Per-cluster + session-level firing-rate plots.
%
% PURPOSE:
%   Owns the plotting that used to live inside calculate_fireRate_general
%   (audit item S). Iterates clusters, calls plot_single_fireRate for
%   each, then calls plot_multi_fireRate once at the session level.
%   Decoupling the plot calls from the analysis function lets NGL02 (and
%   any other caller) compute firing rate without paying the plotting
%   cost when they only want the numbers.
%
% USAGE:
%   plot_fireRate_session(fireRate, param, opt)
%   Called by NGL02_postPhy after calculate_fireRate_general (per area in
%   multi-area mode).
%
% INPUTS:
%   fireRate - struct from calculate_fireRate_general. Required:
%                .sps      {Nclust x 1} of [Ntrials x Nbins] (raw)
%                .Norm     {Nclust x 1} of [Ntrials x Nbins] (normalised)
%                .meanNorm {Nclust x 1} of [1 x Nbins] (mean normalised)
%   param    - param struct from calculate_fireRate_general. Used:
%                .plot   gate; if false this function is a no-op.
%                .cl     [alignmentIdx clusterIdx levelIdx]; updated
%                        per-cluster inside this function.
%                .ROI, .KSLabel, .HumanLabel, .bc_unitType, .phyLabel
%                        per-cluster metadata forwarded to plot subtitles.
%   opt      - resolved options struct; forwarded to the plot helpers
%                (alignto, area, analysis, SavFileName, etc.).
%
% NOTES:
%   - This function does not recompute anything. It assumes the caller
%     already has the per-cluster cells in fireRate populated.
%   - For multi-alignment runs, today this only renders alignment index 1
%     (matches the previous internal behaviour, where param.cl(1)=a was
%     set inside the cluster loop but the function only iterated `a=1`
%     by convention). When alignment-aware plotting is required, the
%     caller should iterate alignments before calling.
%
% SEE ALSO:
%   calculate_fireRate_general, plot_single_fireRate, plot_multi_fireRate.
%
% Last modified 02.06.2026 (Jesus) - extracted from calculate_fireRate_general (#10 S)

%% No-op gate.
if isfield(param,'plot') && ~param.plot, return; end

%% Per-cluster plots.
Nclust = numel(fireRate.sps);
for c = 1:Nclust
    param.cl = [1 c 1];  % alignment 1, cluster c, level 1
    plot_single_fireRate(fireRate.sps{c,1}, fireRate.Norm{c,1}, param, opt);
end

%% Session-level multi-cluster plot.
plot_multi_fireRate(fireRate.meanNorm, param, opt);
end
