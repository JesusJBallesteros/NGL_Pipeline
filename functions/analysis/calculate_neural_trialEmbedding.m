function result = calculate_neural_trialEmbedding(neurons, fireRate, condition, opt) %#ok<INUSL>
% calculate_neural_trialEmbedding  Per-trial population-state embedding via
%                                  PCA, t-SNE, or UMAP. The legacy
%                                  "neural dynamics" view, renamed to
%                                  reflect what it actually does.
%
% PURPOSE:
%   For every trial, build a feature vector by flattening the trial's
%   smoothed firing-rate matrix (Nclust x Nbins) into a single row, then
%   project all trials into a low-dim embedding so trial-to-trial
%   structure (drift, condition clusters, outliers) can be inspected
%   visually. This is NOT a neural-trajectory analysis — time is folded
%   into the feature axis. For trajectory-through-time analyses, use
%   calculate_neural_pca (and, eventually, jPCA/GPFA).
%
% USAGE:
%   result = calculate_neural_trialEmbedding(neurons, fireRate, opt)
%
% INPUTS:
%   neurons   - struct from sort2trials (used only to size things).
%   fireRate  - struct from calculate_fireRate_general. Required:
%                 .sps  {Nclust x 1} of [Ntrials x Nbins] matrices.
%   opt       - resolved options struct. Used subfields:
%                 .popDyn.trialEmbedMethod  'PCA' | 'tSNE' | 'UMAP'
%                 .popDyn.smoothSigma       Gaussian sigma (seconds)
%                 .popDyn.nComponents       output dimensionality (>=2)
%                 .analysis, .area, .alignto, .SavFileName (for title/save)
%
% OUTPUT:
%   result    - struct with:
%                 .method        'PCA' | 'tSNE' | 'UMAP'
%                 .reducedData   [Ntrials x K]
%                 .explained     [K x 1] (PCA only; [] for tSNE/UMAP)
%                 .area          area tag
%
% CHANGES FROM THE OLD calculate_neural_dynamics.m:
%   - The synth-Poisson surrogate has been removed. Real data only. If a
%     null-distribution comparison is ever needed, build a dedicated
%     calculate_neural_null_surrogate function.
%   - The dead "real-data" branch is now the only branch and is built on
%     fireRate.sps directly (no more rebinning from raw spike times with
%     inconsistent units).
%   - The misnamed "trajectories" plot has been removed (the function
%     never produced trajectories; trajectories live in
%     calculate_neural_pca).
%   - Figure saved to opt.analysis/plots/population_dynamics/, not pwd.
%   - Returns a populated struct instead of an unassigned output.
%
% SEE ALSO:
%   calculate_population_dynamics, calculate_neural_pca, smooth_spikes,
%   fireRate_to_tensor.
%
% Last modified 29.05.2026 (Jesus) - port from calculate_neural_dynamics (#18)

method = opt.popDyn.trialEmbedMethod;
K      = max(2, opt.popDyn.nComponents);

%% Build the [Nclust x Nbins x Ntrials] tensor and smooth along time.
if isfield(opt,'stepSz_ms'), stepSz_ms = opt.stepSz_ms; else, stepSz_ms = 20; end
binSize_s = stepSz_ms / 1000;

% NOTE (#26): fireRate.sps is {Nclust x Nalign}; pick alignment via
% opt.popDyn.alignIdx (default 1).
alignIdx   = opt.popDyn.alignIdx;
rateTensor = fireRate_to_tensor(fireRate, alignIdx);
rateTensor = smooth_spikes(rateTensor, opt.popDyn.smoothSigma, binSize_s);

%% Optional: drop aborted trials (#19, default true via opt.popDyn.dropAborted).
% Mirrors the historic calculate_fireRate_general filter for consistency
% across the popDyn family.
if isfield(opt.popDyn,'dropAborted') && opt.popDyn.dropAborted ...
        && isstruct(condition) && isfield(condition,'aborted')
    validMask  = applyTrialFilter(condition, 'allInitiated');
    rateTensor = rateTensor(:, :, validMask);
end
[Nclust, Nbins, Ntrials] = size(rateTensor);

% Flatten each trial into a single feature vector -> [Ntrials x (Nclust*Nbins)].
X = reshape(permute(rateTensor, [3 1 2]), Ntrials, Nclust * Nbins);

%% Dispatch.
explained = [];
switch lower(method)
    case 'pca'
        Kavail = min(K, min(size(X)));
        if Kavail < K
            warning('NGL:trialEmbedding:fewerComponents', ...
                'Requested %d components but only %d available.', K, Kavail);
        end
        [~, score, ~, ~, explained] = pca(X, 'NumComponents', Kavail);
        reducedData = score;

    case 'tsne'
        rng('default');
        perplexity = min(30, max(5, floor(Ntrials/4)));
        reducedData = tsne(X, 'NumDimensions', K, 'Perplexity', perplexity, ...
                           'Standardize', true);

    case 'umap'
        if exist('run_umap','file') ~= 2
            error('NGL:UMAPmissing', ...
                ['opt.popDyn.trialEmbedMethod=''UMAP'' but run_umap is not on the ', ...
                 'MATLAB path. Install UMAP for MATLAB from File Exchange #71902.']);
        end
        [reducedData, ~, ~] = run_umap(X, 'n_components', K);

    otherwise
        error('NGL:trialEmbedding:badMethod', ...
            'opt.popDyn.trialEmbedMethod=''%s'' is not recognised.', method);
end

%% Plot.
areaTag = ''; if isfield(opt,'area'),    areaTag = opt.area; end
align   = ''; if isfield(opt,'alignto') && numel(opt.alignto) >= alignIdx, align = opt.alignto{alignIdx}; end
titleStr = sprintf('Trial-state embedding (%s)  |  area %s  |  align %s  |  %d trials', ...
                   method, areaTag, align, Ntrials);

fig = figure('Visible','off','Position',[100 100 700 600]);
if size(reducedData,2) >= 3
    scatter3(reducedData(:,1), reducedData(:,2), reducedData(:,3), 30, ...
             1:Ntrials, 'filled');
    xlabel('Dim 1'); ylabel('Dim 2'); zlabel('Dim 3');
    view(3);
else
    scatter(reducedData(:,1), reducedData(:,2), 30, 1:Ntrials, 'filled');
    xlabel('Dim 1'); ylabel('Dim 2');
end
title(titleStr); grid on; box off;
cb = colorbar; cb.Label.String = 'trial index';

if isfield(opt,'analysis') && ~isempty(opt.analysis)
    outDir = fullfile(opt.analysis, 'plots', 'population_dynamics');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    stem = '';
    if isfield(opt,'SavFileName'), stem = opt.SavFileName; end
    exportgraphics(fig, fullfile(outDir, sprintf('%s_trialEmbed_%s.png', stem, method)));
end
close(fig);

%% Pack result.
result.method      = method;
result.reducedData = reducedData;
result.explained   = explained;
result.area        = areaTag;
end
