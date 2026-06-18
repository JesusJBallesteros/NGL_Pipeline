%% NGL04_fireRate. Cross-subject, cross-session FR PSTH plotter.
%
% PURPOSE:
%   Consume the aggregated cells produced by NGL03_acrossSession and
%   render averaged firing-rate line plots (mean + error shade) via
%   plotPSTH from toolboxes/BDPAT_NGL. A user-set `request` triplet
%   selects which condition, cluster-label, and alignment to pool over.
%   Exactly one of the three may be a 'X vs Y' comparison; the plot
%   layout is chosen accordingly.
%
%   Multi-area aware: if aggregated.allspike{x,y} carries per-area
%   nested sub-structs (matching entries in input.Areas), the script
%   iterates over those areas and produces one set of outputs per area,
%   with the area name appended to every filename. Single-area runs
%   behave exactly as before.
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
%   01.  Load aggregated cells
%   01a. Multi-area discovery (input.Areas)
%   for each area:
%     02.  Build category sets for inference on the per-area view
%     03.  Parse `request` into per-factor level lists
%     04.  Pool per-trial spike cells across matching (subj, sess, c)
%     05.  Plot:
%            varying ALIGNMENT -> two subplots side-by-side
%            else              -> overlaid traces, different colors
%     06.  Save PNG (area-suffixed in multi-area mode)
%     07.  Example waveforms PNG (same naming)
%
% OUTPUT (workspace + on disk):
%   result - struct with:
%       .areas         cell of area tags processed (= {''} in single-area)
%       .byArea.<area> per-area sub-struct with .pooled / .levels /
%                      .upperY / .figFile / .figFile_waveforms
%     In single-area mode the per-area fields are ALSO mirrored at the
%     top level (result.pooled, result.figFile, ...) for legacy callers.
%   PNGs - <input.analysis>/plots/fireRate/
%            <encoded-request>.png                 (single area)
%            <encoded-request>_<area>.png          (multi-area)
%            <encoded-request>[_<area>]_waveforms.png
%
% DEPENDENCIES:
%   plotPSTH, calcFireRate, nanMeanSterrHistogram (toolboxes/BDPAT_NGL/);
%   aggregated.mat (or per-subject files) from NGL03_acrossSession;
%   functions/analysis/{loadAggregatedSpikes, buildRequestCatSets,
%   parseFireRateRequest, buildFireRatePool, encodeFireRateRequest,
%   requestSubplotTitle, requestTraceLabel,
%   detectMultiAreaFields, flattenAggregatedForArea}.
%
% Last modified 18.06.2026 (Jesus) - multi-area support: per-area
%                                     iteration + area-tagged outputs +
%                                     area in firepools cache key.

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

%% 01. Load aggregated cells.
% Prefer the study-level aggregated.mat (subject-aggregation). Fall back
% to per-subject files if only session-aggregation ran.
aggregated = loadAggregatedSpikes(input);
% Required fields in `aggregated`: allspike, allneurons, allcondition.
% Optional: allevents, alltrialdef.

%% 01a. Multi-area discovery.
% If aggregated.allspike{*}'s cells carry per-area nested sub-structs
% matching input.Areas, iterate over those areas. Otherwise run once
% with areaTag = '' (single-area / legacy path).
areasToRun = {};
if isfield(input,'Areas') && ~isempty(input.Areas)
    areasToRun = detectMultiAreaFields(aggregated, input.Areas);
end
if isempty(areasToRun)
    areasToRun = {''};           % single-area / flat aggregated
end

% Resolve shared paths once.
cacheDir = localFireRatePlotField(opt, 'cacheDir', '');
if isempty(cacheDir)
    cacheDir = fullfile(input.analysis, 'cache', 'firepools');
end
sourceFile = fullfile(input.analysis, 'aggregated.mat');
if ~isfile(sourceFile), sourceFile = ''; end

outDir = fullfile(input.analysis, 'plots', 'fireRate');
if ~exist(outDir, 'dir'), mkdir(outDir); end
fname = encodeFireRateRequest(request);

result        = struct();
result.areas  = areasToRun;
result.byArea = struct();

isMultiArea = numel(areasToRun) > 1 || (numel(areasToRun) == 1 && ~isempty(areasToRun{1}));
if isMultiArea
    fprintf('NGL04_fireRate: multi-area mode, iterating over %s\n', ...
            strjoin(areasToRun, ', '));
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
        fprintf('\nNGL04_fireRate: ===== area %s =====\n', areaTag);
    end

    %% 02. Build category sets used by content-based category inference.
    catSets = buildRequestCatSets(aggView, opt);

    %% 03. Parse `request` into per-factor level lists.
    parsed = parseFireRateRequest(request, catSets);

    nA = numel(parsed.alignment);
    nC = numel(parsed.condition);
    nL = numel(parsed.label);
    fprintf('NGL04_fireRate: %d alignment(s) x %d condition(s) x %d label(s); varying = {%s}.\n', ...
            nA, nC, nL, strjoin(parsed.varying, ', '));

    %% 04. Pool per-trial spike cells across matching (subj, sess, c).
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
                [pool, cHit] = loadFireRatePoolCache(cacheDir, cKey, sourceFile);
                if ~cHit
                    pool = buildFireRatePool(aggView, aL, cL, lL, lF);
                    saveFireRatePoolCache(cacheDir, cKey, pool);
                end
                pooled{aIdx, cIdx, lIdx} = pool;
                hitTag = ternaryChar(cHit, '[cache]', '[built]');
                fprintf('  %s [%s] (%s | %s | %s [%s]): %d trials | %d clust | %d sess | %d subj\n', ...
                        hitTag, areaLbl, aL, cL, lL, lF, ...
                        pool.nTrials/pool.nClust, pool.nClust, pool.nSess, pool.nSubj);
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
                        'No matching trials for (%s | %s | %s | %s); skipping this trace.', ...
                        areaLbl, parsed.alignment{aIdx}, parsed.condition{cIdx}, parsed.label{lIdx});
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
            traceLabels{end+2, 1} = sprintf('%s; N: %d, n: %d, c: %d, tr: %d', ...
                baseLabel, pool.nSubj, pool.nSess, pool.nClust, pool.nTrials/pool.nClust);
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

    % Suptitle includes the area in multi-area mode.
    if isMultiArea
        sgtitle(fig, sprintf('area %s', areaTag));
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
    figFile = fullfile(outDir, [fname areaSlug '.png']);
    exportgraphics(fig, figFile, 'Resolution', 300);
    close(fig);
    fprintf('NGL04_fireRate: %s\n', figFile);

    %% 07. Example waveforms — separate diagnostic figure.
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
                        legHandles(end+1) = hLast; %#ok<SAGROW>
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
        if isMultiArea, sgtitle(figWF, sprintf('Waveforms | area %s', areaTag)); end
        figFile_waveforms = fullfile(outDir, [fname areaSlug '_waveforms.png']);
        exportgraphics(figWF, figFile_waveforms, 'Resolution', 300);
        close(figWF);
        fprintf('NGL04_fireRate: %s\n', figFile_waveforms);
    else
        warning('NGL04:noWaveforms', ...
            ['No contributing cluster in area %s carries a .waveform field; skipping ', ...
             'example-waveforms diagnostic figure. (Set opt.getwF=true in NGL01.)'], areaLbl);
    end

    %% Stash per-area outputs onto result.
    areaKey = areaTag;
    if isempty(areaKey), areaKey = 'all'; end
    result.byArea.(areaKey) = struct( ...
        'pooled',            {pooled}, ...
        'levels',            levels,   ...
        'upperY',            upperY,   ...
        'figFile',           figFile,  ...
        'figFile_waveforms', figFile_waveforms, ...
        'parsed',            parsed);

    % Back-compat mirror for single-area / legacy callers.
    if ~isMultiArea
        result.pooled            = pooled;
        result.levels            = levels;
        result.upperY            = upperY;
        result.figFile           = figFile;
        result.figFile_waveforms = figFile_waveforms;
    end
end

% Local helpers (UI-only utilities). Anything related to the request
% parsing / pool building lives under functions/analysis/ and is shared
% with NGL04_PCA.

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
% Tiny inline ternary for tick-label assembly.
    if cond, s = a; else, s = b; end
end

function picks = localPickRandom(N, k)
% Return a random size-min(N,k) index vector into 1:N (no replacement).
    if N <= k
        picks = 1:N;
    else
        picks = sort(randperm(N, k));
    end
end
