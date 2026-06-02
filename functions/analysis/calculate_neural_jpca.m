function result = calculate_neural_jpca(neurons, fireRate, condition, opt) %#ok<INUSD>
% calculate_neural_jpca  PLACEHOLDER for jPCA (Churchland 2012).
%
% PURPOSE (when implemented):
%   Identify 2-D planes in PC space with strong rotational dynamics in
%   trial-averaged neural population activity. Useful for motor planning,
%   perceptual decision-making, and other paradigms where preparatory or
%   evolving states rotate through state space.
%
% STATUS:
%   NOT YET IMPLEMENTED. This file is a placeholder so the
%   calculate_population_dynamics wrapper can route opt.popDyn.jPCA=true
%   to a stable call site that emits a clear warning. Activate by
%   completing roadmap task: see docs/audit_calculate_neural_dynamics.md §6.
%
% PLAN (for the eventual implementation):
%   1. Bring Churchland Lab's MATLAB jPCA codepack into toolboxes/jPCA/.
%      Reference port: https://github.com/bantin/jPCA (Python).
%   2. Build a [Nclust x Nbins x Ncond] trial-averaged smoothed-rate
%      tensor (reuse smooth_spikes + fireRate_to_tensor + condition
%      grouping from calculate_neural_pca).
%   3. Stack conditions into the format jPCA expects (Data array of
%      structs with .A = [Nbins x Nclust] per condition).
%   4. Call the jPCA core: returns rotational planes and projected
%      trajectories.
%   5. Plot rotational planes (PC1' vs PC2', PC3' vs PC4') per condition,
%      save to opt.analysis/plots/population_dynamics/.
%   6. Return result.method='jPCA' with the rotation matrix, projected
%      trajectories, eigenvalues of the skew-symmetric projection.
%
% USAGE (today):
%   result = calculate_neural_jpca(neurons, fireRate, condition, opt)
%     -> issues NGL:notImplemented warning, returns a stub struct
%        result.method='jPCA', result.status='not_implemented'.
%
% SEE ALSO:
%   calculate_population_dynamics, calculate_neural_pca, smooth_spikes.
%
% Last modified 29.05.2026 (Jesus)

warning('NGL:notImplemented', ...
    ['calculate_neural_jpca is a placeholder. jPCA is not yet integrated; ', ...
     'see docs/audit_calculate_neural_dynamics.md §6 and the roadmap task ', ...
     'for the implementation plan. Returning stub.']);

result            = struct();
result.method     = 'jPCA';
result.status     = 'not_implemented';
result.references = {
    'Churchland MM, Cunningham JP, Kaufman MT, Foster JD, Nuyujukian P, Ryu SI, Shenoy KV. (2012). Neural population dynamics during reaching. Nature 487:51–56.', ...
    'Reference codepack: Churchland Lab MATLAB jPCA distribution.', ...
    'Python port: https://github.com/bantin/jPCA' };
end
