function result = calculate_neural_pca(neurons, fireRate, condition, opt) %#ok<INUSL>
% calculate_neural_pca  Population-PCA on trial-averaged, smoothed firing rates.
%                       Produces time-resolved neural-state trajectories,
%                       one per condition (or one overall if no condition
%                       grouping is requested).
%
% PURPOSE:
%   The canonical "first-look" neural-dynamics view used throughout the
%   systems-neuroscience literature: smooth each cluster's binned firing
%   rate, average across trials within each condition, project the
%   resulting [Nclust x Nbins x Ncond] tensor onto the leading principal
%   components, and plot the projection as one trajectory per condition
%   through PC space. Time runs along the trajectory.
%
% USAGE:
%   result = calculate_neural_pca(neurons, fireRate, condition, opt)
%
% INPUTS:
%   neurons   - struct from sort2trials (kept in signature; only used to
%               check that ROI labels are available for the title).
%   fireRate  - struct from calculate_fireRate_general. Required:
%                 .sps   {Nclust x 1} of [Ntrials x Nbins] firing-rate
%                        matrices, ALREADY BINNED and rate-normalised by
%                        calcFireRate.
%   condition - per-trial condition struct (from conditions_script). If
%               opt.popDyn.conditionVar names a field, trials with the
%               same value of condition.(conditionVar) are grouped into
%               one trajectory. If conditionVar is '' or the field is
%               absent, all valid trials become a single trajectory.
%   opt       - resolved options struct. Used subfields:
%                 .popDyn.smoothSigma   Gaussian sigma in SECONDS
%                 .popDyn.nComponents   number of PCs to keep (>=2)
%                 .popDyn.conditionVar  field name on `condition` for grouping
%                 .alignto              {1} used for the figure title
%                 .area                 used for the figure title
%                 .analysis             save directory
%                 .SavFileName          filename stem
%   param     - (optional) only param.binSize and param.stepSz are read
%               to recover the time axis; both fall back to calcFireRate
%               defaults (binSize=200 ms, stepSz=20 ms) if absent.
%
% OUTPUT:
%   result    - struct with:
%                 .method        'PCA'
%                 .components    [Nclust x K] principal axes
%                 .explained     [K x 1] percent variance per PC
%                 .scores        [Nbins x K x Ncond] projected trajectories
%                 .timeAxis      [1 x Nbins] time in seconds (or relative)
%                 .conditions    {Ncond x 1} cell of condition labels
%                 .meanFR        the trial-averaged [Nclust x Nbins x Ncond]
%                                tensor before projection (handy for QC)
%
% PLOTS:
%   - 2-D projection on first two PCs, one line per condition (time as
%     line progression; markers at start/end).
%   - 3-D projection on first three PCs (if nComponents>=3).
%   - Saved to opt.analysis/plots/population_dynamics/.
%
% NOTES:
%   - This function is the modern replacement for the trial-state view
%     in the legacy calculate_neural_dynamics.m. The trial-state view
%     lives on as calculate_neural_trialEmbedding.m (also called from
%     calculate_population_dynamics when opt.popDyn.trialEmbed=true).
%   - Multi-area awareness lives in the calling wrapper, not here. This
%     function operates on a single area's fireRate at a time.
%
% SEE ALSO:
%   calculate_population_dynamics, smooth_spikes, fireRate_to_tensor,
%   calculate_neural_jpca, calculate_neural_gpfa, calculate_neural_trialEmbedding.
%
% Last modified 29.05.2026 (Jesus)

%% Recover bin width to size the smoothing kernel.
% calcFireRate uses param.stepSz (ms) as the bin step in fireRate.sps. We
% read it back from the saved param if present on opt; otherwise default.
if isfield(opt,'stepSz_ms'),  stepSz_ms = opt.stepSz_ms;
else,                          stepSz_ms = 20;
end
binSize_s = stepSz_ms / 1000;

K = opt.popDyn.nComponents;

%% Build [Nclust x Nbins x Ntrials] tensor and smooth along time.
rateTensor = fireRate_to_tensor(fireRate);
rateTensor = smooth_spikes(rateTensor, opt.popDyn.smoothSigma, binSize_s);
[Nclust, Nbins, Ntrials] = size(rateTensor);

%% Group trials by condition.
% conditionVar names a field on `condition` that gives a per-trial label.
% Empty / missing field -> one "all" group with every trial.
condVar = opt.popDyn.conditionVar;
if isempty(condVar) || ~isstruct(condition) || ~isfield(condition, condVar)
    groupLabels = repmat({'all'}, Ntrials, 1);
else
    raw = condition.(condVar);
    if iscell(raw)
        groupLabels = raw(:);
    elseif isnumeric(raw) || islogical(raw)
        groupLabels = arrayfun(@(v) sprintf('%g', v), raw(:), 'uni', 0);
    elseif iscategorical(raw)
        groupLabels = cellstr(raw(:));
    elseif isstring(raw)
        groupLabels = cellstr(raw(:));
    else
        warning('NGL:calculate_neural_pca:badCondType', ...
            ['condition.%s is of unsupported type %s; falling back to a ', ...
             'single group.'], condVar, class(raw));
        groupLabels = repmat({'all'}, Ntrials, 1);
    end
end

% Drop trials labelled NaN/'' (treated as "exclude").
isUsable = ~cellfun(@(v) isempty(v) || (ischar(v) && any(strcmp(v, {'NaN','<missing>'}))), groupLabels);
groupLabels = groupLabels(isUsable);
rateTensor  = rateTensor(:, :, isUsable);

[uniqueGroups, ~, gIdx] = unique(groupLabels, 'stable');
Ncond = numel(uniqueGroups);

%% Trial-average per group -> [Nclust x Nbins x Ncond].
meanFR = zeros(Nclust, Nbins, Ncond);
for g = 1:Ncond
    sel = (gIdx == g);
    if any(sel)
        meanFR(:, :, g) = mean(rateTensor(:, :, sel), 3, 'omitnan');
    end
end

%% Reshape for PCA: neurons are observations (cols), (time x cond) are rows.
% pca(X) treats X rows as observations, cols as variables. We want PCs
% over the neuron axis, so X is [(Nbins*Ncond) x Nclust].
X = reshape(permute(meanFR, [2 3 1]), Nbins * Ncond, Nclust);

Kavail = min(K, min(size(X)));
if Kavail < K
    warning('NGL:calculate_neural_pca:fewerComponents', ...
        'Requested %d components but only %d available; using %d.', K, Kavail, Kavail);
end
[coeff, scoreMat, ~, ~, explained] = pca(X, 'NumComponents', Kavail);

% Reshape scores back to [Nbins x Kavail x Ncond] for per-condition trajectories.
scores = reshape(scoreMat, Nbins, Ncond, Kavail);
scores = permute(scores, [1 3 2]);   % [Nbins x Kavail x Ncond]

%% Recover a time axis. calcFireRate centres bins around alignment via
%   param.interval; without that info, return a relative index in seconds.
timeAxis = (0:Nbins-1) * binSize_s;

%% Plot trajectories.
if isfield(opt,'analysis') && ~isempty(opt.analysis)
    outDir = fullfile(opt.analysis, 'plots', 'population_dynamics');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
else
    outDir = '';
end

areaTag = ''; if isfield(opt,'area'), areaTag = opt.area; end
align   = ''; if isfield(opt,'alignto') && ~isempty(opt.alignto), align = opt.alignto{1}; end
titleStr = sprintf('Population PCA  |  area %s  |  align %s  |  %d conditions, %d/%d PCs', ...
                   areaTag, align, Ncond, Kavail, K);

cmap = lines(Ncond);

% 2-D figure (PC1 vs PC2). Always available.
fig2 = figure('Visible','off','Position',[100 100 700 600]);
hold on;
for g = 1:Ncond
    x = scores(:, 1, g);
    y = scores(:, 2, g);
    plot(x, y, '-', 'Color', cmap(g,:), 'LineWidth', 1.5);
    plot(x(1),   y(1),   'o', 'MarkerFaceColor', cmap(g,:), 'MarkerEdgeColor', 'k');
    plot(x(end), y(end), 's', 'MarkerFaceColor', cmap(g,:), 'MarkerEdgeColor', 'k');
end
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(min(2,end))));
title(titleStr); legend(uniqueGroups,'Location','bestoutside'); grid on; box off;
if ~isempty(outDir)
    exportgraphics(fig2, fullfile(outDir, sprintf('%s_pca_PC1PC2.png', getfield_default(opt,'SavFileName','session'))));
end
close(fig2);

% 3-D figure when available.
if Kavail >= 3
    fig3 = figure('Visible','off','Position',[100 100 700 600]);
    hold on;
    for g = 1:Ncond
        plot3(scores(:,1,g), scores(:,2,g), scores(:,3,g), '-', ...
              'Color', cmap(g,:), 'LineWidth', 1.5);
    end
    xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
    title(titleStr); legend(uniqueGroups,'Location','bestoutside'); grid on; view(3);
    if ~isempty(outDir)
        exportgraphics(fig3, fullfile(outDir, sprintf('%s_pca_3D.png', getfield_default(opt,'SavFileName','session'))));
    end
    close(fig3);
end

%% Pack result.
result.method     = 'PCA';
result.components = coeff;
result.explained  = explained(1:Kavail);
result.scores     = scores;
result.timeAxis   = timeAxis;
result.conditions = uniqueGroups;
result.meanFR     = meanFR;
result.area       = areaTag;
end

function v = getfield_default(s, fld, default)
    if isfield(s, fld) && ~isempty(s.(fld)), v = s.(fld); else, v = default; end
end
