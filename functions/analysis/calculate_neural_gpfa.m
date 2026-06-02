function result = calculate_neural_gpfa(neurons, fireRate, condition, opt) %#ok<INUSD>
% calculate_neural_gpfa  PLACEHOLDER for GPFA (Yu et al. 2009).
%
% PURPOSE (when implemented):
%   Gaussian-Process Factor Analysis: jointly smooth spike trains and
%   reduce dimensionality so single-trial neural trajectories can be
%   extracted (no need to average across trials). Useful when
%   trial-to-trial variability is the question itself, or when the
%   paradigm has too few trials per condition for trial-averaged PCA to
%   be meaningful.
%
% STATUS:
%   NOT YET IMPLEMENTED. This file is a placeholder so the
%   calculate_population_dynamics wrapper can route opt.popDyn.GPFA=true
%   to a stable call site that emits a clear warning. Activate by
%   completing roadmap task: see docs/audit_calculate_neural_dynamics.md §6.
%
% PLAN (for the eventual implementation):
%   1. Bring a MATLAB GPFA codepack into toolboxes/gpfa/. Candidates:
%      - Byron Yu Lab DataHigh-bundled GPFA
%        (https://users.ece.cmu.edu/~byronyu/software.shtml).
%      - https://github.com/wrongu/gpfa
%      - https://github.com/aecker/gpfa
%   2. Convert per-cluster spike timestamps into the GPFA input format
%      (struct array with .spikes binary matrices per trial).
%   3. Expose key GPFA knobs in opt.popDyn (binWidth, xDim latent
%      dimensionality, kernSD), validated in set_default.
%   4. Run GPFA; collect latent trajectories per single trial.
%   5. Plot single-trial trajectories, coloured by condition; optionally
%      also a condition-averaged trajectory.
%   6. Return result.method='GPFA' with the latent trajectories, model
%      parameters, and the bin width used.
%
% USAGE (today):
%   result = calculate_neural_gpfa(neurons, fireRate, condition, opt)
%     -> issues NGL:notImplemented warning, returns a stub struct.
%
% SEE ALSO:
%   calculate_population_dynamics, calculate_neural_pca, calculate_neural_jpca.
%
% Last modified 29.05.2026 (Jesus)

warning('NGL:notImplemented', ...
    ['calculate_neural_gpfa is a placeholder. GPFA is not yet integrated; ', ...
     'see docs/audit_calculate_neural_dynamics.md §6 and the roadmap task ', ...
     'for the implementation plan. Returning stub.']);

result            = struct();
result.method     = 'GPFA';
result.status     = 'not_implemented';
result.references = {
    'Yu BM, Cunningham JP, Santhanam G, Ryu SI, Shenoy KV, Sahani M. (2009). Gaussian-process factor analysis for low-dimensional single-trial analysis of neural population activity. J. Neurophysiol. 102:614–35.', ...
    'Byron Yu Lab software: https://users.ece.cmu.edu/~byronyu/software.shtml', ...
    'MATLAB implementations: https://github.com/wrongu/gpfa  |  https://github.com/aecker/gpfa' };
end
