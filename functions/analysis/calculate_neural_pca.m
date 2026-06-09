function result = calculate_neural_pca(neurons, fireRate, condition, opt) %#ok<INUSL>
% calculate_neural_pca  Per-session population PCA via the shared NGL04
%                       infrastructure: real single-trial trajectories
%                       overlaid on the condition mean, plus a trial-
%                       bootstrap CI tube. Iterates over every alignment
%                       in opt.alignto, every distinct cluster-label value
%                       observed in the current session, and every entry
%                       in opt.popDyn.pcaConditions.
%
% PURPOSE:
%   Replaces the old trial-averaged PCA with the same engine NGL04_PCA
%   uses cross-subject (buildSingleSessionPool + calculate_pca_from_pool
%   + plot_pca_state_space). At single-session level the per-trial
%   projection becomes meaningful (each grey trace is a real recorded
%   trial), so the 'singleTrials' figure variant automatically uses
%   traj_trial instead of the cross-session marginal fallback.
%
% USAGE (called from calculate_population_dynamics):
%   result = calculate_neural_pca(neurons, fireRate, condition, opt);
%
% INPUTS:
%   neurons    - struct from sort2trials. Required:
%                  .<alignName>{c,1}   {Ntrials x 1} cell of spike vectors (ms)
%                  .<labelField>{c}    per-cluster label value, for one of
%                                       HumanLabel / KSLabel / bc_unitType /
%                                       phyLabel (resolved via
%                                       opt.fireRatePlot.labelPriority).
%                  .ROI{c}             area tag (forwarded to titles).
%   fireRate   - struct from calculate_fireRate_general. KEPT IN
%                SIGNATURE for back-compat but no longer consumed; the
%                new path re-bins from neurons.<alignName> directly.
%   condition  - per-trial condition struct. Must contain .aborted for
%                'allInitiated'; any other token requires the named field.
%   opt        - resolved options struct. Used fields:
%                  .alignto               cell of alignment names
%                  .binSize_ms, .stepSz_ms
%                  .popDyn.smoothSigma
%                  .popDyn.nComponents
%                  .popDyn.pcaConditions  cell of condition tokens (single
%                                          or 'X vs Y'); default
%                                          {'allInitiated'}.
%                  .fireRatePlot.interval window around alignment, ms
%                  .fireRatePlot.labelPriority resolution order for which
%                                          label field carries the
%                                          per-cluster identity
%                  .pcaPlot.{nBootstrap, rngSeed, sessionAlpha, ciAlpha,
%                            ciStride, variants, outDir}
%                  .analysis              base dir for output plots
%                  .area                  area tag (multi-area runs)
%                  .SavFileName           filename stem for the session
%
% OUTPUT:
%   result    - struct:
%                 .method       'PCA-pool'
%                 .labelField   which label field was iterated
%                 .labelValues  cell of label values iterated
%                 .alignments   echo of opt.alignto
%                 .pcaConditions echo of opt.popDyn.pcaConditions
%                 .scans        {Nalign x Nlabel x Nconds} cell of PCA
%                               results from calculate_pca_from_pool.
%                               Empty cells indicate "no matching
%                               clusters/trials" combinations.
%                 .files        cell of saved PNG paths
%                 .area         opt.area
%
% PLOTS:
%   For each non-empty (alignment, label, conditionEntry) scan, every
%   variant in opt.pcaPlot.variants ('singleTrials' and/or 'ciTube')
%   produces one PNG under
%       <opt.analysis>/plots/population_dynamics/
%       <SavFileName>_pca_<align>_<label>_<condEntry>_<variant>.png
%
% NOTES:
%   - Multi-area awareness lives in the calling wrapper: NGL02 sets
%     opt.area per area before calling this. We just forward it into the
%     title prefix.
%   - opt.popDyn.alignIdx is IGNORED here on purpose: this implementation
%     iterates all alignments so the user sees the full picture in one
%     pass. The old single-alignment behaviour is gone.
%
% SEE ALSO:
%   buildSingleSessionPool, calculate_pca_from_pool, plot_pca_state_space,
%   calculate_population_dynamics.
%
% Last modified 09.06.2026 (Jesus)

result = struct();
result.method        = 'PCA-pool';

%% Resolve the label field via priority order.
% Pick the first label field that's populated on the current neurons
% struct. ROI is never used as a label here (it's the area tag).
labelPriority = opt.fireRatePlot.labelPriority;
labelField    = '';
for f = labelPriority
    if isfield(neurons, f{1}) && ~isempty(neurons.(f{1}))
        labelField = f{1};
        break
    end
end
if isempty(labelField)
    warning('NGL:calculate_neural_pca:noLabels', ...
        ['No populated label field on neurons (checked %s). Falling back to ', ...
         'treating every cluster as labelValue=''all''.'], ...
        strjoin(labelPriority, ' > '));
    labelField   = 'pseudoLabel';
    pseudoVals   = repmat({'all'}, numel(neurons.ROI), 1);
    neurons.pseudoLabel = pseudoVals;
end

% Distinct label values present in this session.
rawLabels = neurons.(labelField);
rawLabels = rawLabels(~cellfun(@isempty, rawLabels));
rawLabels = cellfun(@(v) char(string(v)), rawLabels, 'uni', false);
labelValues = unique(rawLabels);
if isempty(labelValues)
    warning('NGL:calculate_neural_pca:noLabelValues', ...
        'No usable label values found on neurons.%s; nothing to plot.', labelField);
    result.labelField   = labelField;
    result.labelValues  = {};
    result.alignments   = opt.alignto;
    result.pcaConditions = opt.popDyn.pcaConditions;
    result.scans        = {};
    result.files        = {};
    return
end

condEntries = opt.popDyn.pcaConditions;
nA = numel(opt.alignto);
nL = numel(labelValues);
nC = numel(condEntries);

%% Output directory.
if isfield(opt,'pcaPlot') && isfield(opt.pcaPlot,'outDir') && ~isempty(opt.pcaPlot.outDir)
    outDir = opt.pcaPlot.outDir;
else
    outDir = fullfile(opt.analysis, 'plots', 'population_dynamics');
end
if ~isfolder(outDir), mkdir(outDir); end

stem = getfield_default(opt, 'SavFileName', 'session');
areaTag = ''; if isfield(opt,'area'), areaTag = opt.area; end

%% Shared PCA params (same for every scan in this session).
pcaParams = struct( ...
    'intervalMs',    opt.fireRatePlot.interval, ...
    'binSize_ms',    opt.binSize_ms,            ...
    'stepSz_ms',     opt.stepSz_ms,             ...
    'smoothSigma_s', opt.popDyn.smoothSigma,    ...
    'nComponents',   opt.popDyn.nComponents,    ...
    'nBootstrap',    opt.pcaPlot.nBootstrap,    ...
    'smpRate',       1000,                      ...
    'rngSeed',       opt.pcaPlot.rngSeed,       ...
    'includeTrials', true);

variants = opt.pcaPlot.variants;
nV       = numel(variants);

%% Iterate alignments × labels × condition entries.
result.labelField    = labelField;
result.labelValues   = labelValues(:);
result.alignments    = opt.alignto(:);
result.pcaConditions = condEntries(:);
result.scans         = cell(nA, nL, nC);
result.files         = cell(nA, nL, nC, nV);
result.area          = areaTag;

for aIdx = 1:nA
    alignName = opt.alignto{aIdx};
    if ~isfield(neurons, alignName)
        warning('NGL:calculate_neural_pca:missingAlign', ...
            'neurons has no alignment ''%s''; skipping.', alignName);
        continue
    end
    for lIdx = 1:nL
        labelValue = labelValues{lIdx};
        for cIdx = 1:nC
            condEntry = strtrim(condEntries{cIdx});
            % 'X vs Y' -> two parts; single -> one part.
            condParts = regexp(condEntry, '\s+vs\s+', 'split', 'ignorecase');
            condParts = cellfun(@strtrim, condParts, 'uni', false);

            % Build one pool per condition part on this session.
            condPools = cell(numel(condParts), 1);
            for p = 1:numel(condParts)
                condPools{p} = buildSingleSessionPool( ...
                    neurons, condition, alignName, condParts{p}, ...
                    labelValue, labelField);
            end

            % Skip if no cluster matched in any part.
            if all(cellfun(@(pp) pp.nClust == 0, condPools))
                continue
            end

            try
                pcaResult = calculate_pca_from_pool(condPools, condParts, pcaParams);
            catch ME
                warning('NGL:calculate_neural_pca:scanFailed', ...
                    '(%s | %s | %s): %s', alignName, labelValue, condEntry, ME.message);
                continue
            end
            result.scans{aIdx, lIdx, cIdx} = pcaResult;

            titlePref = sprintf('PCA | %s | %s | %s', alignName, labelValue, condEntry);
            if ~isempty(areaTag)
                titlePref = sprintf('%s | area %s', titlePref, areaTag);
            end
            for v = 1:nV
                variant = variants{v};
                areaSlug = sanitize(areaTag);
                if isempty(areaSlug), areaSlug = 'area'; end
                tag      = sprintf('%s_pca_%s_%s_%s_%s_%s.png', stem, ...
                             areaSlug, sanitize(alignName), sanitize(labelValue), ...
                             sanitize(condEntry), variant);
                outFile  = fullfile(outDir, tag);
                plot_pca_state_space(pcaResult, struct( ...
                    'variant',      variant,                  ...
                    'titlePrefix',  titlePref,                ...
                    'outFile',      outFile,                  ...
                    'sessionAlpha', opt.pcaPlot.sessionAlpha, ...
                    'ciAlpha',      opt.pcaPlot.ciAlpha,      ...
                    'ciStride',     opt.pcaPlot.ciStride));
                result.files{aIdx, lIdx, cIdx, v} = outFile;
                fprintf('calculate_neural_pca: %s\n', outFile);
            end
        end
    end
end

end

function v = getfield_default(s, fld, default)
    if isfield(s, fld) && ~isempty(s.(fld)), v = s.(fld); else, v = default; end
end

function s = sanitize(s)
    s = regexprep(char(string(s)), '[^\w\-]', '_');
    if isempty(s), s = '_'; end
end
