function neuralDynamics = calculate_population_dynamics(neurons, fireRate, trialdef, condition, opt)
% calculate_population_dynamics  Wrapper that dispatches enabled
%                                population-dim-reduction methods over
%                                one area's spike data and aggregates the
%                                results into a single struct.
%
% PURPOSE:
%   NGL02 calls this once per area (multi-area mode) or once for the
%   session (single-area mode). It looks at the opt.popDyn.<method>
%   flags, runs each enabled method, and collects everything into
%       neuralDynamics.pca         (calculate_neural_pca)
%       neuralDynamics.jPCA        (calculate_neural_jpca; placeholder)
%       neuralDynamics.GPFA        (calculate_neural_gpfa; placeholder)
%       neuralDynamics.trialEmbed  (calculate_neural_trialEmbedding)
%   Methods that error out are caught, logged via warning, and reported
%   in their result struct so one broken method doesn't bring the rest
%   down with it.
%
% USAGE:
%   neuralDynamics = calculate_population_dynamics( ...
%                       neurons, fireRate, trialdef, condition, opt)
%
% INPUTS:
%   neurons   - struct from sort2trials (standard mode). Passed through
%               to the method functions; some only use the shape.
%   fireRate  - struct from calculate_fireRate_general; primary input
%               (already binned and rate-normalised by calcFireRate).
%   trialdef  - cell from trialdefGen; available to methods that need it
%               (kept in signature; not used by every method).
%   condition - per-trial condition struct; used by the method functions
%               for per-condition grouping when
%               opt.popDyn.conditionVar is set.
%   opt       - resolved options struct. opt.popDyn.* selects which
%               methods run; opt.area / opt.alignto / opt.analysis /
%               opt.SavFileName are forwarded to method functions for
%               plot titles and save paths.
%
% OUTPUT:
%   neuralDynamics - struct with one sub-field per method that ran. Each
%                    sub-field is the return value of the method
%                    function (see those headers for shape).
%                    .meta carries the opt.popDyn snapshot used.
%
% MULTI-AREA:
%   This wrapper is area-agnostic. NGL02's doSpikething block loops over
%   input.areaMap.uniqueAreas and calls this wrapper per area, building
%   neuralDynamics.<area> in the caller. Single-area runs call it once
%   with a flat fireRate.
%
% ADDING A NEW METHOD:
%   1. Drop opt.popDyn.<newName>=false into default_opt and validate in
%      set_default.
%   2. Add a calculate_neural_<newName>.m alongside the others.
%   3. Add a dispatch block below; copy one of the existing ones.
%
% SEE ALSO:
%   calculate_neural_pca, calculate_neural_jpca, calculate_neural_gpfa,
%   calculate_neural_trialEmbedding.
%
% Last modified 29.05.2026 (Jesus)

neuralDynamics      = struct();
neuralDynamics.meta = opt.popDyn;

ran = {};

%% PCA — trial-averaged smoothed-rate PCA with time-resolved trajectories.
if opt.popDyn.pca
    try
        neuralDynamics.pca = calculate_neural_pca(neurons, fireRate, condition, opt);
        ran{end+1} = 'pca';
    catch ME
        warning('NGL:popDyn:pcaFailed', ...
            'calculate_neural_pca failed: %s', ME.message);
        neuralDynamics.pca = struct('method','PCA','status','error', ...
                                    'message', ME.message);
    end
end

%% jPCA — rotational dynamics. Currently a placeholder.
if opt.popDyn.jPCA
    try
        neuralDynamics.jPCA = calculate_neural_jpca(neurons, fireRate, condition, opt);
        ran{end+1} = 'jPCA';
    catch ME
        warning('NGL:popDyn:jpcaFailed', ...
            'calculate_neural_jpca failed: %s', ME.message);
        neuralDynamics.jPCA = struct('method','jPCA','status','error', ...
                                     'message', ME.message);
    end
end

%% GPFA — single-trial smooth trajectories. Currently a placeholder.
if opt.popDyn.GPFA
    try
        neuralDynamics.GPFA = calculate_neural_gpfa(neurons, fireRate, condition, opt);
        ran{end+1} = 'GPFA';
    catch ME
        warning('NGL:popDyn:gpfaFailed', ...
            'calculate_neural_gpfa failed: %s', ME.message);
        neuralDynamics.GPFA = struct('method','GPFA','status','error', ...
                                     'message', ME.message);
    end
end

%% Trial-state embedding — legacy view (per-trial points in low-dim space).
if opt.popDyn.trialEmbed
    try
        neuralDynamics.trialEmbed = calculate_neural_trialEmbedding(neurons, fireRate, condition, opt);
        ran{end+1} = 'trialEmbed';
    catch ME
        warning('NGL:popDyn:trialEmbedFailed', ...
            'calculate_neural_trialEmbedding failed: %s', ME.message);
        neuralDynamics.trialEmbed = struct('method', opt.popDyn.trialEmbedMethod, ...
                                           'status','error','message', ME.message);
    end
end

%% Report.
if isempty(ran)
    warning('NGL:popDyn:nothingEnabled', ...
        ['opt.popDyn.do=true but no method flag was on. Enable at least one of ', ...
         'opt.popDyn.pca / jPCA / GPFA / trialEmbed.']);
else
    fprintf('calculate_population_dynamics: ran %s\n', strjoin(ran, ', '));
end
end
