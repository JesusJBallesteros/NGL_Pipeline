%% NGL04_fireRate. Per-subject (× per-area) FR PSTH plotter.
%
% PURPOSE:
%   Consume the aggregated cells produced by NGL03_acrossSession and
%   render averaged firing-rate line plots (mean + error shade) via
%   plotPSTH from toolboxes/BDPAT_NGL. A user-set `request` triplet
%   selects which condition, cluster-label, and alignment to pool over.
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
%   every filename. Single-area runs behave exactly as before.
%
% USAGE (from NGL_SetAndRunMe section 3.x):
%   request = {'correct vs incorrect', 'good', 'stim2'};
%   NGL04_fireRate
%
%   request = {'correct', 'itiOn vs stim2', 'good'};
%   NGL04_fireRate
%
%   request = {'correct', 'good vs mua', 'stim2'};
%   NGL04_fireRate
%
% INPUT (workspace variable):
%   request - 1x3 cell. Each entry is one of:
%               * single value (e.g. 'correct', 'good', 'stim2')
%               * comparison   (e.g. 'correct vs incorrect')
%             At most ONE entry may contain ' vs '. Categories are
%             inferred from content so the order is flexible.
%
% PIPELINE:
%   00.  NGL00_Prep + Areas recovery + set_default + findSessions
%   01.  Set up paths / output dir / cache dir (no aggregated load yet)
%   for each subject:
%     02.  Try cached parsed-request; load aggregated only if miss
%     03.  Parse `request` into per-factor level lists
%     for each area:
%       04.  Pool per-trial spike cells (cache HIT or build from aggregated)
%       05.  Plot (alignment -> subplots, condition x label -> overlays)
%       06.  Save PNG (subject + area suffixed)
%       07.  Example waveforms PNG (same naming)
%
% OUTPUT (workspace + on disk):
%   result - struct with:
%       .subjects     cell of subject names processed
%       .areas        cell of area tags processed (={''} for single-area)
%       .bySubject.(subj).byArea.(area)  per-(subj,area) sub-struct with
%             .pooled / .levels / .upperY / .figFile / .figFile_waveforms
%     For single-subject runs the per-area fields are mirrored at
%     result.byArea.<area>; for single-subject + single-area runs they
%     are also mirrored at result.{pooled,figFile,...} for legacy callers.
%   PNGs - <input.analysis>/plots/fireRate/
%            <encoded-request>_<subject>[_<area>].png
%            <encoded-request>_<subject>[_<area>]_waveforms.png
%
% CACHE LAYOUT (changed 23.06.2026):
%   <cacheDir>/<subject>/
%       <encoded-request>__parsed.mat       parsed request struct
%       <align>__<cond>__<label>_<labelField>[__area_<NCL>].mat   pool struct
%   Subject A and B no longer share cache files. If you upgraded from
%   the pre-23-June layout (<cacheDir>/<key>.mat), the old files will
%   simply never be hit and you can delete them.
%
% DEPENDENCIES:
%   plotPSTH, calcFireRate, nanMeanSterrHistogram (toolboxes/BDPAT_NGL/);
%   aggregated.mat (or per-subject files) from NGL03_acrossSession;
%   functions/analysis/{loadAggregatedSpikes, buildRequestCatSets,
%   parseFireRateRequest, buildFireRatePool, encodeFireRateRequest,
%   requestSubplotTitle, requestTraceLabel, restrictAggregatedToSubject,
%   loadParsedRequestCache, saveParsedRequestCache,
%   detectMultiAreaFields, flattenAggregatedForArea}.
%
% Last modified 23.06.2026 (Jesus) - per-subject iteration; per-subject
%                                     cache folder (subject in path);
%                                     lazy aggregated load skipped when
%                                     parsed + all pools are cached.

%% 00. Standard scaffolding.
NGL00_Prep

% Areas recovery
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
            warning('NGL04:preflight', 'Could not read master preprocInfo: %s', ME.message);
        end
    end
end

[input, opt] = set_default(input, opt);
input.sessions = findSessions(input);

assert(exist('request','var') == 1 && iscell(request) && numel(request) == 3, ...
    'NGL04:badRequest', ...
    ['NGL04_fireRate requires a workspace cell `request` of length 3, ', ...
     'e.g. request = {''correct vs incorrect'', ''good'', ''stim2''}.']);

%% 01. Shared paths + lazy aggregated handle.
% aggregated stays [] until the first cache miss; if every parsed-
% request + pool for every (subject, area) is on disk, the big
% aggregated.mat is never loaded.
baseCacheDir = localFireRatePlotField(opt, 'cacheDir', '');
if isempty(baseCacheDir)
    baseCacheDir = fullfile(input.analysis, 'cache', 'firepools');
end
sourceFile = fullfile(input.analysis, 'aggregated.mat');
if ~isfile(sourceFile), sourceFile = ''; end

outDir = fullfile(input.analysis, 'plots', 'fireRate');
if ~exist(outDir, 'dir'), mkdir(outDir); end
fname = encodeFireRateRequest(request);

% Area discovery from input.Areas (no aggregated probe needed). The
% smart detector in detectMultiAreaFields was useful when input.Areas
% wasn't carried by preprocInfo; the recovery block above now ensures
% input.Areas is set whenever the original NGL01 ran multi-area.
if isfield(input,'Areas') && ~isempty(input.Areas)
    areasToRun = unique(input.Areas, 'stable');
else
    areasToRun = {''};
end
isMultiArea = numel(areasToRun) > 1 || (numel(areasToRun) == 1 && ~isempty(areasToRun{1}));

subjects    = input.subjects;
nSubj       = numel(subjects);
isMultiSubj = nSubj > 1;

result          = struct();
result.subjects = {subjects.name};
result.areas    = areasToRun;
result.bySubject = struct();

% Workspace-aware lazy aggregated handle. If a previous interactive
% run already loaded aggregated.mat into the base workspace, reuse it
% when its fingerprint (studyName + analysis path + subject count +
% source mtime) matches the current input. Otherwise reset to [] so
% the first cache miss triggers a fresh load. This makes
% "tweak request → re-run NGL04" iterations free.
if exist('aggregated','var') && isstruct(aggregated) ...
        && isfield(aggregated, 'srcFingerprint') ...
        && isequaln(aggregated.srcFingerprint, localAggFingerprint(input))
    fprintf('NGL04_fireRate: reusing aggregated already in workspace (fingerprint match).\n');
else
    aggregated = [];
end

if isMultiSubj
    fprintf('NGL04_fireRate: %d subjects in input.subjects; will iterate per subject.\n', nSubj);
end
if isMultiArea
    fprintf('NGL04_fireRate: multi-area mode, iterating over %s\n', strjoin(areasToRun, ', '));
end

%% Outer subject loop
for sIdx = 1:nSubj
    subj         = subjects(sIdx).name;
    subjCacheDir = fullfile(baseCacheDir, subj);
    if ~isfolder(subjCacheDir), mkdir(subjCacheDir); end
    fprintf('\nNGL04_fireRate: ===== subject %s =====\n', subj);

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
    fprintf('NGL04_fireRate: %d alignment(s) x %d condition(s) x %d label(s); varying = {%s}.\n', ...
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
            fprintf('\nNGL04_fireRate: ----- area %s -----\n', areaTag);
        end

        %% 04. Pool cells (cache HIT or build).
        pooled = cell(nA, nC, nL);
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
                    pooled{aIdx, cIdx, lIdx} = pool;
                    hitTag = ternaryChar(cHit, '[cache]', '[built]');
                    fprintf('  %s [%s/%s] (%s | %s | %s [%s]): %d trials | %d clust | %d sess\n', ...
                            hitTag, subj, areaLbl, aL, cL, lL, lF, ...
                            pool.nTrials, pool.nClust, pool.nSess);
                end
            end
        end

        %% 05. Plot via plotPSTH.
        plotCfg = struct();
        plotCfg.interval  = localFireRatePlotField(opt, 'interval',   [-500 4000]);
        plotCfg.binSize   = localFireRatePlotField(opt, 'binSize_ms', opt.binSize_ms);
        plotCfg.stepSz    = localFireRatePlotField(opt, 'stepSz_ms',  opt.stepSz_ms);
        plotCfg.smooth    = localFireRatePlotField(opt, 'smoothPlot', true);
        plotCfg.errAlpha  = localFireRatePlotField(opt, 'errAlpha',   0.4);
        plotCfg.smpRate   = 1000;
        plotCfg.busyWarn  = localFireRatePlotField(opt, 'busyWarnTraces', 4);

        nTracesPerSubplot = nC * nL;
        plotCfg.palette   = lines(max(2, nTracesPerSubplot));
        if nTracesPerSubplot > plotCfg.busyWarn
            warning('NGL04:busyPlot', ...
                'Request produces %d overlaid traces per subplot (threshold %d). Consider narrowing the request.', ...
                nTracesPerSubplot, plotCfg.busyWarn);
        end

        fig    = figure('Visible','off','Position',[100 100 max(900, 450*nA) 500]);
        upperY = nan(nA, nC, nL);
        axHandles = gobjects(nA, 1);

        intervalMs = plotCfg.interval;
        stepMs     = plotCfg.stepSz;
        sLo        = floor(intervalMs(1)/1000);
        sHi        = ceil(intervalMs(2)/1000);
        tickSec    = sLo:0.5:sHi;
        tickBins   = (tickSec*1000 - intervalMs(1)) / stepMs + 1;
        tickLab    = arrayfun(@(t) ternaryChar(abs(mod(t,1))<1e-9, sprintf('%g', t), ''), ...
                              tickSec, 'uni', false);
        alignBin   = (0 - intervalMs(1)) / stepMs + 1;

        for aIdx = 1:nA
            if nA > 1
                axHandles(aIdx) = subplot(1, nA, aIdx);
            else
                axHandles(aIdx) = gca;
            end
            hold on
            traceLabels = cell(0,1);
            for cIdx = 1:nC
                for lIdx = 1:nL
                    traceIdx = (cIdx-1) * nL + lIdx;
                    pool     = pooled{aIdx, cIdx, lIdx};
                    if isempty(pool.cells)
                        warning('NGL04:emptyCell', ...
                            'No matching trials for (%s | %s | %s | %s | %s); skipping trace.', ...
                            subj, areaLbl, parsed.alignment{aIdx}, parsed.condition{cIdx}, parsed.label{lIdx});
                        continue
                    end
                    upY = plotPSTH(pool.cells, ...
                                  plotCfg.stepSz, plotCfg.binSize, ...
                                  plotCfg.interval, plotCfg.smpRate, ...
                                  'plotcol',    plotCfg.palette(traceIdx, :), ...
                                  'meanline',   '-', ...
                                  'smoothplot', plotCfg.smooth, ...
                                  'erralpha',   plotCfg.errAlpha);
                    upperY(aIdx, cIdx, lIdx) = upY;
                end
                baseLabel = requestTraceLabel(parsed.varying, parsed, cIdx, lIdx);
                traceLabels{end+2, 1} = sprintf('%s; n: %d, c: %d, tr: %d', ...
                    baseLabel, pool.nSess, pool.nClust, pool.nTrials);
            end

            ax = axHandles(aIdx);
            ax.XTick      = tickBins;
            ax.XTickLabel = tickLab;
            title(requestSubplotTitle(parsed, aIdx));
            xlabel('t since event (s)');
            if aIdx == 1
                ylabel('spikes/s');
            end
            xl = xline(alignBin, '--k', parsed.alignment{aIdx});
            xl.LabelHorizontalAlignment = 'right';
            xl.LabelVerticalAlignment   = 'top';
            xl.LabelOrientation         = 'horizontal';
            xl.FontSize                 = 10;
            if numel(traceLabels) >= 1 && aIdx == nA
                traceLabels(1:2:3) = {''};
                legend(traceLabels, 'Location', 'south');
                legend Box off
            end
            box off
            hold off
        end

        % Suptitle: subject (and area in multi-area).
        if isMultiArea
            sgtitle(fig, sprintf('%s | area %s', subj, areaTag));
        else
            sgtitle(fig, subj);
        end

        % Harmonise y-axis across all subplots.
        maxY = max(upperY(:), [], 'omitnan');
        if isfinite(maxY) && maxY > 0
            for aIdx = 1:nA
                ylim(axHandles(aIdx), [0, maxY * 1.05]);
                if aIdx > 1
                    axHandles(aIdx).YAxis.Visible = 'off';
                    axHandles(aIdx).YTickLabel = {};
                    axHandles(aIdx).YLabel.String = '';
                end
            end
        end

        %% 06. Save main PSTH PNG.
        figFile = fullfile(outDir, [fname '_' subj areaSlug '.png']);
        exportgraphics(fig, figFile, 'Resolution', 300);
        close(fig);
        fprintf('NGL04_fireRate: %s\n', figFile);

        %% 07. Example waveforms diagnostic (same naming).
        figFile_waveforms = '';
        anyWF = false;
        for k = 1:numel(pooled)
            if ~isempty(pooled{k}.waveforms), anyWF = true; break, end
        end
        if anyWF
            figWF = figure('Visible','off','Position',[100 100 max(900, 450*nA) 320]);
            for aIdx = 1:nA
                if nA > 1, subplot(1, nA, aIdx); end
                hold on
                legHandles = gobjects(0); legNames = {};
                for cIdx = 1:nC
                    for lIdx = 1:nL
                        traceIdx = (cIdx-1) * nL + lIdx;
                        pool     = pooled{aIdx, cIdx, lIdx};
                        if isempty(pool.waveforms), continue, end
                        picks    = localPickRandom(numel(pool.waveforms), 4);
                        col      = plotCfg.palette(traceIdx, :);
                        hLast    = [];
                        for k = picks
                            wf = pool.waveforms{k};
                            hLast = plot(wf, 'Color', [col, 0.5], 'LineWidth', 1);
                        end
                        if ~isempty(hLast)
                            legHandles(end+1) = hLast;
                            legNames{end+1}   = sprintf('%s (%d out of %d)', ...
                                requestTraceLabel(parsed.varying, parsed, cIdx, lIdx), ...
                                numel(picks), numel(pool.waveforms));
                        end
                    end
                end
                title(['Waveforms | ' parsed.alignment{aIdx}]);
                xlabel('sample'); ylabel('amplitude');
                if ~isempty(legHandles)
                    legend(legHandles, legNames, 'Location', 'southeast');
                    legend Box off
                end
                box off
                hold off
            end
            if isMultiArea
                sgtitle(figWF, sprintf('Waveforms | %s | area %s', subj, areaTag));
            else
                sgtitle(figWF, sprintf('Waveforms | %s', subj));
            end
            figFile_waveforms = fullfile(outDir, [fname '_' subj areaSlug '_waveforms.png']);
            exportgraphics(figWF, figFile_waveforms, 'Resolution', 300);
            close(figWF);
            fprintf('NGL04_fireRate: %s\n', figFile_waveforms);
        end

        %% Stash per-(subject, area) outputs onto result.
        areaKey = areaTag; if isempty(areaKey), areaKey = 'all'; end
        areaResult = struct( ...
            'pooled',            {pooled}, ...
            'levels',            levels,   ...
            'upperY',            upperY,   ...
            'figFile',           figFile,  ...
            'figFile_waveforms', figFile_waveforms, ...
            'parsed',            parsed);
        result.bySubject.(subj).byArea.(areaKey) = areaResult;

        % Single-subject single-area legacy mirror.
        if ~isMultiSubj && ~isMultiArea
            result.pooled            = pooled;
            result.levels            = levels;
            result.upperY            = upperY;
            result.figFile           = figFile;
            result.figFile_waveforms = figFile_waveforms;
        end
    end

    % Single-subject convenience mirror.
    if ~isMultiSubj
        result.byArea = result.bySubject.(subj).byArea;
    end
end

%% Local helpers (UI-only utilities). Anything related to request
% parsing / pool building lives under functions/analysis/ and is
% shared with NGL04_PCA.

function aggregated = localEnsureAggregated(aggregated, input)
% Lazy-load aggregated.mat on the first cache miss. Once loaded, the
% handle is reused across all subjects + areas in the same run AND
% (because the script preserves it in the base workspace) across
% subsequent runs whose fingerprint matches.
    if isempty(aggregated)
        fprintf('NGL04_fireRate: loading aggregated.mat (first cache miss this run)...\n');
        aggregated = loadAggregatedSpikes(input);
        aggregated.srcFingerprint = localAggFingerprint(input);
    end
end

function fp = localAggFingerprint(input)
% Cheap identity tag stamped on aggregated. Compared on subsequent
% NGL04 runs to decide whether the in-workspace copy is still valid.
% sourceMtime catches the "user re-ran NGL03 since last load" case.
    fp = struct( ...
        'studyName',   '', ...
        'analysis',    '', ...
        'nSubj',       0, ...
        'sourceMtime', NaN);
    if isfield(input,'studyName'), fp.studyName = input.studyName; end
    if isfield(input,'analysis'),  fp.analysis  = input.analysis;  end
    if isfield(input,'subjects'),  fp.nSubj     = numel(input.subjects); end
    src = fullfile(fp.analysis, 'aggregated.mat');
    if isfile(src)
        d = dir(src);
        if ~isempty(d), fp.sourceMtime = d.datenum; end
    end
end

function v = localFireRatePlotField(opt, fname, dflt)
% Read opt.fireRatePlot.(fname) with a safe default fallback so we
% don't crash when the user partially populated opt.fireRatePlot.
    if isfield(opt,'fireRatePlot') && isfield(opt.fireRatePlot, fname)
        v = opt.fireRatePlot.(fname);
    else
        v = dflt;
    end
end

function s = ternaryChar(cond, a, b)
    if cond, s = a; else, s = b; end
end

function picks = localPickRandom(N, k)
    if N <= k
        picks = 1:N;
    else
        picks = sort(randperm(N, k));
    end
end
