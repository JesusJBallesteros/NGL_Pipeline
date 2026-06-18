%% NGL04_PCA. Cross-subject, cross-session population PCA state-space plotter.
%
% PURPOSE:
%   Sibling of NGL04_fireRate. Takes the same 3-cell `request` and
%   produces neural state-space (PCA) trajectories instead of PSTH lines.
%   For every (alignment, label) cell of the request the script:
%
%     * pools per-trial spike vectors per cluster across (subject, session)
%       reusing the firepools cache written by NGL04_fireRate when the
%       same (alignment, condition, label) combination is available,
%     * fits PCA on the pooled, smoothed condition-mean rates,
%     * projects per-session marginal trajectories into the same basis,
%     * computes a trial-bootstrap CI envelope around each condition mean,
%     * renders TWO figure variants per (alignment, label):
%         singleTrials -> condition mean + per-session grey traces
%         ciTube       -> condition mean + bootstrap CI ribbon (2D) /
%                         per-axis CI crosshairs (3D)
%
%   Multi-area aware: if aggregated cells nest per-area sub-structs
%   matching input.Areas, the script iterates over those areas and
%   produces a full set of PNGs + a result entry per area, with the
%   area name appended to every filename. Single-area runs behave
%   exactly as before.
%
% USAGE (from NGL_SetAndRunMe section 4.2):
%   request = {'correct vs incorrect', 'good', 'stimOn2'};
%   NGL04_PCA
%
%   request = {'correct', 'itiOn vs stimOn2', 'good'};
%   NGL04_PCA
%
% INPUT (workspace variable):
%   request - 1x3 cell with the same semantics as NGL04_fireRate.
%             Conditions become trajectory colours within a PC space.
%             Alignment and label each produce a separate PC space
%             (separate output figure).
%
% PIPELINE:
%   00.  NGL00_Prep + Areas recovery + set_default + findSessions
%   01.  loadAggregatedSpikes
%   01a. Multi-area discovery (input.Areas)
%   for each area:
%     02.  buildRequestCatSets + parseFireRateRequest (per-area view)
%     03.  Build / load firepools cache (area-tagged key, shared with NGL04_fireRate)
%     04.  For each (alignment, label): calculate_pca_from_pool over conditions
%     05.  Plot two variants per result via plot_pca_state_space
%   06.  Save PCA results + figure paths as <fname>[_<area>]_pca.mat
%
% OUTPUT (workspace + on disk):
%   result - struct with:
%       .areas         cell of area tags processed (= {''} in single-area)
%       .byArea.<area> per-area sub-struct with .pools / .levels / .pca / .files / .matFile
%     In single-area mode the per-area fields are ALSO mirrored at the
%     top level for legacy callers.
%   PNGs   - <input.analysis>/plots/PCA/<encoded-request>[_<area>]_aA_lL_<variant>.png
%   MAT    - <input.analysis>/plots/PCA/<encoded-request>[_<area>]_pca.mat
%
% DEPENDENCIES:
%   functions/analysis/{loadAggregatedSpikes, buildRequestCatSets,
%   parseFireRateRequest, buildFireRatePool, fireRatePoolCacheKey,
%   loadFireRatePoolCache, saveFireRatePoolCache, encodeFireRateRequest,
%   calculate_pca_from_pool, requestSubplotTitle,
%   detectMultiAreaFields, flattenAggregatedForArea};
%   functions/plotting/plot_pca_state_space;
%   toolboxes/BDPAT_NGL/calcFireRate.
%
% Last modified 18.06.2026 (Jesus) - multi-area support: per-area
%                                     iteration + area-tagged outputs +
%                                     area in firepools cache key.

%% 00. Standard scaffolding.
NGL00_Prep

% Areas recovery (same fallback pattern as NGL04_fireRate)
if ~isfield(input,'Areas') || isempty(input.Areas)
    analysisCodePath = fullfile(input.datadrive, input.studyName, 'analysisCode');
    if ~contains(analysisCodePath, ':\') && ~isempty(input.datadrive)
        analysisCodePath = fullfile([input.datadrive(1) ':\'], input.studyName, 'analysisCode');
    end
    try
        masterInfo = loadPreprocInfo(analysisCodePath, 'master');
        if isfield(masterInfo,'Areas') && ~isempty(masterInfo.Areas)
            input.Areas = masterInfo.Areas;
        end
    catch ME
        if ~strcmp(ME.identifier, 'NGL:loadPreprocInfo:notFound')
            warning('NGL04_PCA:preflight', 'Could not read master preprocInfo: %s', ME.message);
        end
    end
end

[input, opt] = set_default(input, opt);
input.sessions = findSessions(input);

assert(exist('request','var') == 1 && iscell(request) && numel(request) == 3, ...
    'NGL04_PCA:badRequest', ...
    ['NGL04_PCA requires a workspace cell `request` of length 3, ', ...
     'e.g. request = {''correct vs incorrect'', ''good'', ''stim2''}.']);

%% 01. Load aggregated cells.
aggregated = loadAggregatedSpikes(input);

%% 01a. Multi-area discovery.
areasToRun = {};
if isfield(input,'Areas') && ~isempty(input.Areas)
    areasToRun = detectMultiAreaFields(aggregated, input.Areas);
end
if isempty(areasToRun)
    areasToRun = {''};
end

% Shared paths once.
cacheDir = opt.fireRatePlot.cacheDir;
if isempty(cacheDir)
    cacheDir = fullfile(input.analysis, 'cache', 'firepools');
end
sourceFile = fullfile(input.analysis, 'aggregated.mat');
if ~isfile(sourceFile), sourceFile = ''; end

outDir = opt.pcaPlot.outDir;
if isempty(outDir)
    outDir = fullfile(input.analysis, 'plots', 'PCA');
end
if ~isfolder(outDir), mkdir(outDir); end
fname_root = encodeFireRateRequest(request);

pcaParams = struct( ...
    'intervalMs',    opt.pcaPlot.interval,      ...
    'binSize_ms',    opt.pcaPlot.binSize_ms,    ...
    'stepSz_ms',     opt.pcaPlot.stepSz_ms,     ...
    'smoothSigma_s', opt.pcaPlot.smoothSigma_s, ...
    'nComponents',   opt.pcaPlot.nComponents,   ...
    'nBootstrap',    opt.pcaPlot.nBootstrap,    ...
    'smpRate',       1000,                      ...
    'rngSeed',       opt.pcaPlot.rngSeed);

result        = struct();
result.areas  = areasToRun;
result.byArea = struct();

isMultiArea = numel(areasToRun) > 1 || (numel(areasToRun) == 1 && ~isempty(areasToRun{1}));
if isMultiArea
    fprintf('NGL04_PCA: multi-area mode, iterating over %s\n', strjoin(areasToRun, ', '));
end

for areaIdx = 1:numel(areasToRun)
    areaTag = areasToRun{areaIdx};
    if isempty(areaTag)
        aggView  = aggregated;
        areaSlug = '';
        areaLbl  = 'all';
    else
        aggView  = flattenAggregatedForArea(aggregated, areaTag);
        areaSlug = ['_' areaTag];
        areaLbl  = areaTag;
        fprintf('\nNGL04_PCA: ===== area %s =====\n', areaTag);
    end

    %% 02. Build category sets + parse request.
    catSets = buildRequestCatSets(aggView, opt);
    parsed  = parseFireRateRequest(request, catSets);

    nA = numel(parsed.alignment);
    nC = numel(parsed.condition);
    nL = numel(parsed.label);
    fprintf('NGL04_PCA: %d alignment(s) x %d condition(s) x %d label(s); varying = {%s}.\n', ...
            nA, nC, nL, strjoin(parsed.varying, ', '));

    %% 03. Build / load pools. SHARED cache with NGL04_fireRate (area-tagged).
    pools  = cell(nA, nC, nL);
    levels = repmat(struct('alignment','','condition','','label',''), nA, nC, nL);
    for aIdx = 1:nA
        for cIdx = 1:nC
            for lIdx = 1:nL
                aL = parsed.alignment{aIdx};
                cL = parsed.condition{cIdx};
                lL = parsed.label{lIdx};
                lF = parsed.labelField{lIdx};
                levels(aIdx, cIdx, lIdx).alignment = aL;
                levels(aIdx, cIdx, lIdx).condition = cL;
                levels(aIdx, cIdx, lIdx).label     = lL;

                cKey         = fireRatePoolCacheKey(aL, cL, lL, lF, areaTag);
                [pool, cHit] = loadFireRatePoolCache(cacheDir, cKey, sourceFile);
                if ~cHit
                    pool = buildFireRatePool(aggView, aL, cL, lL, lF);
                    saveFireRatePoolCache(cacheDir, cKey, pool);
                end
                pools{aIdx, cIdx, lIdx} = pool;
                hitTag = ternaryChar(cHit, '[cache]', '[built]');
                fprintf('  %s [%s] (%s | %s | %s [%s]): %d trials | %d clust | %d sess | %d subj\n', ...
                        hitTag, areaLbl, aL, cL, lL, lF, ...
                        pool.nTrials, pool.nClust, pool.nSess, pool.nSubj);
            end
        end
    end

    %% 04. PCA per (alignment, label) over the Ncond conditions.
    pcaCells  = cell(nA, nL);
    fileCells = cell(nA, nL, numel(opt.pcaPlot.variants));
    for aIdx = 1:nA
        for lIdx = 1:nL
            condPools  = cell(nC, 1);
            condLabels = cell(nC, 1);
            for cIdx = 1:nC
                condPools{cIdx}  = pools{aIdx, cIdx, lIdx};
                condLabels{cIdx} = parsed.condition{cIdx};
            end

            pca = calculate_pca_from_pool(condPools, condLabels, pcaParams);
            pcaCells{aIdx, lIdx} = pca;

            fprintf('NGL04_PCA: (area %s, a=%d, l=%d) PCA done: %d clusters / %d sessions, %d PCs (top 3 var: %.1f%% %.1f%% %.1f%%).\n', ...
                areaLbl, aIdx, lIdx, ...
                numel(pca.clusters), numel(pca.sessionKeys), numel(pca.explained), ...
                pca.explained(1), pca.explained(min(2,end)), pca.explained(min(3,end)));

            %% 05. Plot the configured variants.
            titlePref = sprintf('PCA | %s', requestSubplotTitle(parsed, aIdx));
            if nL > 1
                titlePref = sprintf('%s | %s', titlePref, parsed.label{lIdx});
            end
            if isMultiArea
                titlePref = sprintf('%s | area %s', titlePref, areaTag);
            end
            for v = 1:numel(opt.pcaPlot.variants)
                variant = opt.pcaPlot.variants{v};
                fSuffix = sprintf('%s_a%d_l%d_%s.png', areaSlug, aIdx, lIdx, variant);
                outFile = fullfile(outDir, [fname_root, fSuffix]);
                plot_pca_state_space(pca, struct( ...
                    'variant',      variant,                  ...
                    'titlePrefix',  titlePref,                ...
                    'outFile',      outFile,                  ...
                    'sessionAlpha', opt.pcaPlot.sessionAlpha, ...
                    'ciAlpha',      opt.pcaPlot.ciAlpha,      ...
                    'ciStride',     opt.pcaPlot.ciStride));
                fileCells{aIdx, lIdx, v} = outFile;
                fprintf('NGL04_PCA: %s\n', outFile);
            end
        end
    end

    %% 06. Per-area result + .mat dump.
    matFile = fullfile(outDir, [fname_root, areaSlug, '_pca.mat']);
    areaResult = struct( ...
        'pools',   {pools},   ...
        'levels',  levels,    ...
        'pca',     {pcaCells},...
        'files',   {fileCells},...
        'parsed',  parsed,    ...
        'matFile', matFile);
    save(matFile, 'areaResult', '-v7.3');
    fprintf('NGL04_PCA: %s\n', matFile);

    areaKey = areaTag;
    if isempty(areaKey), areaKey = 'all'; end
    result.byArea.(areaKey) = areaResult;

    % Back-compat mirror for single-area / legacy callers.
    if ~isMultiArea
        result.pools   = pools;
        result.levels  = levels;
        result.pca     = pcaCells;
        result.files   = fileCells;
        result.matFile = matFile;
    end
end

% Local helpers (UI-only).
function s = ternaryChar(cond, a, b)
    if cond, s = a; else, s = b; end
end
