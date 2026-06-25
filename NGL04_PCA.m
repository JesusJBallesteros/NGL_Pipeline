%% NGL04_PCA. Per-subject (× per-area) population PCA state-space plotter.
%
% PURPOSE:
%   Sibling of NGL04_fireRate. Takes the same 3-cell `request` and
%   produces neural state-space (PCA) trajectories instead of PSTH lines.
%   Per (alignment, label) the script:
%
%     * pools per-trial spike vectors per cluster, reusing the firepools
%       cache written by NGL04_fireRate when the same
%       (subject, area, alignment, condition, label) combination is
%       already on disk,
%     * fits PCA on the pooled, smoothed condition-mean rates,
%     * projects per-session marginal trajectories into the same basis,
%     * computes a trial-bootstrap CI envelope around each condition mean,
%     * renders TWO figure variants per (alignment, label):
%         singleTrials -> condition mean + per-session grey traces
%         ciTube       -> condition mean + bootstrap CI ribbon (2D) /
%                         per-axis CI crosshairs (3D)
%
%   Per-subject pooling: iterates each subject in input.subjects and
%   produces one set of outputs per subject. The cache lives in
%   <cacheDir>/<subject>/ so subject A's pools never contaminate
%   subject B's request. Aggregated.mat is loaded LAZILY — when every
%   pool and the parsed-request the run needs are already cached, the
%   large file is never touched.
%
%   Multi-area aware: when input.Areas declares more than one area,
%   the script also iterates per area, with the area name appended to
%   every filename.
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
%   01.  Set up shared paths (no aggregated load yet)
%   for each subject:
%     02.  Try cached parsed-request; load aggregated only if miss
%     03.  Parse `request` into per-factor level lists
%     for each area:
%       04.  Build / load pools (cache HIT or restrict-aggregated build)
%       05.  PCA per (alignment, label) via calculate_pca_from_pool
%       06.  Plot each variant via plot_pca_state_space
%       07.  Save .mat with PCA results for this (subject, area)
%
% OUTPUT (workspace + on disk):
%   result - struct with:
%       .subjects     cell of subject names processed
%       .areas        cell of area tags processed (= {''} for single-area)
%       .bySubject.(subj).byArea.(area)  per-(subj,area) sub-struct with
%             .pools (cell of pool structs) / .levels / .pca / .files /
%             .parsed / .matFile
%     Single-subject runs mirror to result.byArea; single-subject single-
%     area runs additionally mirror to result.{pools, pca, files, matFile}
%     at the top level for legacy callers.
%   PNGs - <input.analysis>/plots/PCA/
%            <encoded-request>_<subject>[_<area>]_aA_lL_<variant>.png
%   MAT  - <input.analysis>/plots/PCA/
%            <encoded-request>_<subject>[_<area>]_pca.mat
%
% CACHE LAYOUT (matches NGL04_fireRate):
%   <cacheDir>/<subject>/
%       <encoded-request>__parsed.mat       parsed request struct
%       <align>__<cond>__<label>_<labelField>[__area_<NCL>].mat   pool struct
%   Pools written by NGL04_fireRate are re-used here and vice-versa.
%
% DEPENDENCIES:
%   functions/analysis/{loadAggregatedSpikes, buildRequestCatSets,
%   parseFireRateRequest, buildFireRatePool, fireRatePoolCacheKey,
%   loadFireRatePoolCache, saveFireRatePoolCache, encodeFireRateRequest,
%   calculate_pca_from_pool, requestSubplotTitle,
%   restrictAggregatedToSubject, loadParsedRequestCache,
%   saveParsedRequestCache, flattenAggregatedForArea};
%   functions/plotting/plot_pca_state_space;
%   toolboxes/BDPAT_NGL/calcFireRate.
%
% Last modified 25.06.2026 (Jesus) - per-subject iteration; per-subject
%                                     cache folder (subject in path);
%                                     lazy aggregated load skipped when
%                                     parsed + all pools are cached.
%                                     Mirrors the NGL04_fireRate pattern.

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

%% 01. Shared paths + lazy aggregated handle.
baseCacheDir = opt.fireRatePlot.cacheDir;
if isempty(baseCacheDir)
    baseCacheDir = fullfile(input.analysis, 'cache', 'firepools');
end
sourceFile = fullfile(input.analysis, 'aggregated.mat');
if ~isfile(sourceFile), sourceFile = ''; end

outDir = opt.pcaPlot.outDir;
if isempty(outDir)
    outDir = fullfile(input.analysis, 'plots', 'PCA');
end
if ~isfolder(outDir), mkdir(outDir); end
fname_root = encodeFireRateRequest(request);

% Area discovery from input.Areas (no aggregated probe needed).
if isfield(input,'Areas') && ~isempty(input.Areas)
    areasToRun = unique(input.Areas, 'stable');
else
    areasToRun = {''};
end
isMultiArea = numel(areasToRun) > 1 || (numel(areasToRun) == 1 && ~isempty(areasToRun{1}));

subjects    = input.subjects;
nSubj       = numel(subjects);
isMultiSubj = nSubj > 1;

% PCA params (constant across subjects / areas).
pcaParams = struct( ...
    'intervalMs',    opt.pcaPlot.interval,      ...
    'binSize_ms',    opt.pcaPlot.binSize_ms,    ...
    'stepSz_ms',     opt.pcaPlot.stepSz_ms,     ...
    'smoothSigma_s', opt.pcaPlot.smoothSigma_s, ...
    'nComponents',   opt.pcaPlot.nComponents,   ...
    'nBootstrap',    opt.pcaPlot.nBootstrap,    ...
    'smpRate',       1000,                      ...
    'rngSeed',       opt.pcaPlot.rngSeed);

result          = struct();
result.subjects = {subjects.name};
result.areas    = areasToRun;
result.bySubject = struct();

aggregated = [];   % lazy

if isMultiSubj
    fprintf('NGL04_PCA: %d subjects in input.subjects; will iterate per subject.\n', nSubj);
end
if isMultiArea
    fprintf('NGL04_PCA: multi-area mode, iterating over %s\n', strjoin(areasToRun, ', '));
end

%% Outer subject loop
for sIdx = 1:nSubj
    subj         = subjects(sIdx).name;
    subjCacheDir = fullfile(baseCacheDir, subj);
    if ~isfolder(subjCacheDir), mkdir(subjCacheDir); end
    fprintf('\nNGL04_PCA: ===== subject %s =====\n', subj);

    %% 02. Parsed-request cache (skip aggregated load when possible).
    [parsed, pHit] = loadParsedRequestCache(subjCacheDir, request, sourceFile);
    if ~pHit
        aggregated = localEnsureAggregated(aggregated, input);
        aggView    = restrictAggregatedToSubject(aggregated, sIdx);
        catSets    = buildRequestCatSets(aggView, opt);
        parsed     = parseFireRateRequest(request, catSets);
        saveParsedRequestCache(subjCacheDir, request, parsed);
        fprintf('  parsed request built from aggregated (cache miss).\n');
    else
        fprintf('  parsed request loaded from cache.\n');
    end

    nA = numel(parsed.alignment);
    nC = numel(parsed.condition);
    nL = numel(parsed.label);
    fprintf('NGL04_PCA: %d alignment(s) x %d condition(s) x %d label(s); varying = {%s}.\n', ...
            nA, nC, nL, strjoin(parsed.varying, ', '));

    result.bySubject.(subj).byArea = struct();

    %% Per-area loop 
    for areaIdx = 1:numel(areasToRun)
        areaTag = areasToRun{areaIdx};
        if isempty(areaTag)
            areaSlug = '';
            areaLbl  = 'all';
        else
            areaSlug = ['_' areaTag];
            areaLbl  = areaTag;
            fprintf('\nNGL04_PCA: ----- area %s -----\n', areaTag);
        end

        %% 03. Build / load pools. SHARED cache layout with NGL04_fireRate.
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
                    [pool, cHit] = loadFireRatePoolCache(subjCacheDir, cKey, sourceFile);
                    if ~cHit
                        aggregated = localEnsureAggregated(aggregated, input);
                        aggView    = restrictAggregatedToSubject(aggregated, sIdx);
                        if ~isempty(areaTag)
                            aggView = flattenAggregatedForArea(aggView, areaTag);
                        end
                        pool = buildFireRatePool(aggView, aL, cL, lL, lF);
                        saveFireRatePoolCache(subjCacheDir, cKey, pool);
                    end
                    pools{aIdx, cIdx, lIdx} = pool;
                    hitTag = ternaryChar(cHit, '[cache]', '[built]');
                    fprintf('  %s [%s/%s] (%s | %s | %s [%s]): %d trials | %d clust | %d sess\n', ...
                            hitTag, subj, areaLbl, aL, cL, lL, lF, ...
                            pool.nTrials, pool.nClust, pool.nSess);
                end
            end
        end

        %% 04-05. PCA per (alignment, label) + plotting.
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

                fprintf('NGL04_PCA: (%s | %s | a=%d, l=%d) PCA done: %d clusters / %d sessions, %d PCs (top 3 var: %.1f%% %.1f%% %.1f%%).\n', ...
                    subj, areaLbl, aIdx, lIdx, ...
                    numel(pca.clusters), numel(pca.sessionKeys), numel(pca.explained), ...
                    pca.explained(1), pca.explained(min(2,end)), pca.explained(min(3,end)));

                %% Plot the configured variants.
                titlePref = sprintf('PCA | %s | %s', subj, requestSubplotTitle(parsed, aIdx));
                if nL > 1
                    titlePref = sprintf('%s | %s', titlePref, parsed.label{lIdx});
                end
                if isMultiArea
                    titlePref = sprintf('%s | area %s', titlePref, areaTag);
                end
                for v = 1:numel(opt.pcaPlot.variants)
                    variant = opt.pcaPlot.variants{v};
                    fSuffix = sprintf('_%s%s_a%d_l%d_%s.png', subj, areaSlug, aIdx, lIdx, variant);
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

        %% 07. Per-(subject, area) result + .mat dump.
        matFile = fullfile(outDir, [fname_root, '_', subj, areaSlug, '_pca.mat']);
        areaResult = struct( ...
            'pools',   {pools},   ...
            'levels',  levels,    ...
            'pca',     {pcaCells},...
            'files',   {fileCells},...
            'parsed',  parsed,    ...
            'matFile', matFile);
        save(matFile, 'areaResult', '-v7.3');
        fprintf('NGL04_PCA: %s\n', matFile);

        areaKey = areaTag; if isempty(areaKey), areaKey = 'all'; end
        result.bySubject.(subj).byArea.(areaKey) = areaResult;

        % Single-subject single-area legacy mirror.
        if ~isMultiSubj && ~isMultiArea
            result.pools   = pools;
            result.levels  = levels;
            result.pca     = pcaCells;
            result.files   = fileCells;
            result.matFile = matFile;
        end
    end

    % Single-subject convenience mirror.
    if ~isMultiSubj
        result.byArea = result.bySubject.(subj).byArea;
    end
end

%% Local helpers (UI-only utilities). Anything related to request
% parsing / pool building lives under functions/analysis/ and is
% shared with NGL04_fireRate.

function aggregated = localEnsureAggregated(aggregated, input)
% Lazy-load aggregated.mat on the first cache miss. Once loaded, the
% handle is reused across all subjects + areas in the same run.
    if isempty(aggregated)
        fprintf('NGL04_PCA: loading aggregated.mat (first cache miss this run)...\n');
        aggregated = loadAggregatedSpikes(input);
    end
end

function s = ternaryChar(cond, a, b)
    if cond, s = a; else, s = b; end
end
