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
%   02.  Build category sets for inference
%   03.  Parse `request` into per-factor level lists
%   04.  Pool per-trial spike cells across all matching (subj, sess, c)
%   05.  Plot:
%          varying ALIGNMENT -> two subplots side-by-side
%          else              -> overlaid traces, different colors
%        Each trace drawn by plotPSTH(pooled, ..., 'plotcol', ...).
%   06.  Save PNG; return `result` struct on workspace.
%
% OUTPUT (workspace + on disk):
%   result - struct with .pooled (cell of pooled spike-time cells per
%            level), .levels (struct describing each level), .meanLines
%            and .errLines (returned by plotPSTH), and .figFile (PNG path).
%   PNG    - <input.analysis>/plots/NGL04_fireRate/<encoded-request>.png
%
% DEPENDENCIES:
%   plotPSTH, calcFireRate, nanMeanSterrHistogram (toolboxes/BDPAT_NGL/);
%   aggregated.mat (or per-subject files) from NGL03_acrossSession;
%   functions/analysis/{loadAggregatedSpikes, buildRequestCatSets,
%   parseFireRateRequest, buildFireRatePool, encodeFireRateRequest,
%   requestSubplotTitle, requestTraceLabel}.
%
% MULTI-AREA NOTE:
%   v1 supports flat (single-area) aggregated cells. For multi-area
%   aggregation (allspike{x,y}.<Area>.KSLabel), add opt.fireRatePlot.area
%   to pick the area and update Section 04.
%
% Last modified 09.06.2026 (Jesus) - extracted request/pool helpers to
%                                     functions/analysis/ for reuse by
%                                     NGL04_PCA (task #30).

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

%% 02. Build category sets used by content-based category inference.
catSets = buildRequestCatSets(aggregated, opt);

%% 03. Parse `request` into per-factor level lists.
%
% Each of the three factors gets a cell of 1 OR 2 entries depending on
% whether that slot was a single value or a 'X vs Y' comparison. ZERO,
% ONE, TWO, or THREE factors can vary simultaneously (revised
% 02.06.2026).
%
% parsed.alignment {1x1 or 1x2} cell of alignment names (members of opt.alignto)
% parsed.condition {1x1 or 1x2} cell of condition field names
% parsed.label     {1x1 or 1x2} cell of label values (with .field tag)
% parsed.varying   cell of factor names that have >1 levels (e.g. {'alignment','condition'})
parsed = parseFireRateRequest(request, catSets);

nA = numel(parsed.alignment);
nC = numel(parsed.condition);
nL = numel(parsed.label);
fprintf('NGL04_fireRate: %d alignment(s) x %d condition(s) x %d label(s); varying = {%s}.\n', ...
        nA, nC, nL, strjoin(parsed.varying, ', '));

%% 04. Pool per-trial spike cells across all matching (subj, sess, c).
%
% Cartesian product over the three factors. Each (aIdx, cIdx, lIdx)
% cell receives a pool struct. Per-pool cache lives in
% opt.fireRatePlot.cacheDir (default <input.analysis>/cache/firepools)
% and is shared with NGL04_PCA; staleness checked against the
% aggregated.mat mtime.
cacheDir = localFireRatePlotField(opt, 'cacheDir', '');
if isempty(cacheDir)
    cacheDir = fullfile(input.analysis, 'cache', 'firepools');
end
sourceFile = fullfile(input.analysis, 'aggregated.mat');
if ~isfile(sourceFile), sourceFile = ''; end

result        = struct();
result.pooled = cell(nA, nC, nL);
result.levels = repmat(struct('alignment','','condition','','label',''), nA, nC, nL);
for aIdx = 1:nA
    for cIdx = 1:nC
        for lIdx = 1:nL
            aL  = parsed.alignment{aIdx};
            cL  = parsed.condition{cIdx};
            lL  = parsed.label{lIdx};
            lF  = parsed.labelField{lIdx};
            result.levels(aIdx, cIdx, lIdx).alignment = aL;
            result.levels(aIdx, cIdx, lIdx).condition = cL;
            result.levels(aIdx, cIdx, lIdx).label     = lL;

            cKey         = fireRatePoolCacheKey(aL, cL, lL, lF);
            [pool, cHit] = loadFireRatePoolCache(cacheDir, cKey, sourceFile);
            if ~cHit
                pool = buildFireRatePool(aggregated, aL, cL, lL, lF);
                saveFireRatePoolCache(cacheDir, cKey, pool);
            end
            result.pooled{aIdx, cIdx, lIdx} = pool;
            hitTag = ternaryChar(cHit, '[cache]', '[built]');
            fprintf('  %s (%s | %s | %s [%s]): %d trials | %d clust | %d sess | %d subj\n', ...
                    hitTag, aL, cL, lL, lF, pool.nTrials, pool.nClust, pool.nSess, pool.nSubj);
        end
    end
end

%% 05. Plot via plotPSTH.
%
% LAYOUT RULE:
%   ALIGNMENT controls SUBPLOT layout (side-by-side, one per alignment
%   level, because time references differ across alignments).
%   CONDITION x LABEL controls OVERLAY within each subplot (one trace
%   per (condition, label) combination).
%
% Units: spike times in aggregated.allneurons are MILLISECONDS, so pass
% smpRate = 1000 to plotPSTH; binSize/stepSz/interval are also in ms.
plotCfg = struct();
plotCfg.interval  = localFireRatePlotField(opt, 'interval',   [-500 4000]);
plotCfg.binSize   = localFireRatePlotField(opt, 'binSize_ms', opt.binSize_ms);
plotCfg.stepSz    = localFireRatePlotField(opt, 'stepSz_ms',  opt.stepSz_ms);
plotCfg.smooth    = localFireRatePlotField(opt, 'smoothPlot', true);
plotCfg.errAlpha  = localFireRatePlotField(opt, 'errAlpha',   0.4);
plotCfg.smpRate   = 1000;       % aggregated spike times are ms; treat ms as "samples"
plotCfg.busyWarn  = localFireRatePlotField(opt, 'busyWarnTraces', 4);

nTracesPerSubplot = nC * nL;
plotCfg.palette   = lines(max(2, nTracesPerSubplot));
if nTracesPerSubplot > plotCfg.busyWarn
    warning('NGL04:busyPlot', ...
        'Request produces %d overlaid traces per subplot (threshold %d). Consider narrowing the request.', ...
        nTracesPerSubplot, plotCfg.busyWarn);
end

fig = figure('Visible','off','Position',[100 100 max(900, 450*nA) 500]);
result.upperY = nan(nA, nC, nL);
axHandles     = gobjects(nA, 1);

% Pre-compute seconds-based tick positions (in bin units, which is what
% plotPSTH's axes use).
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
            traceIdx = (cIdx-1) * nL + lIdx;    % stable per-(c,l) colour index
            pool     = result.pooled{aIdx, cIdx, lIdx};
            if isempty(pool.cells)
                warning('NGL04:emptyCell', ...
                    'No matching trials for (%s | %s | %s); skipping this trace.', ...
                    parsed.alignment{aIdx}, parsed.condition{cIdx}, parsed.label{lIdx});
                continue
            end
            upperY = plotPSTH(pool.cells, ...
                              plotCfg.stepSz, plotCfg.binSize, ...
                              plotCfg.interval, plotCfg.smpRate, ...
                              'plotcol',    plotCfg.palette(traceIdx, :), ...
                              'meanline',   '-', ...
                              'smoothplot', plotCfg.smooth, ...
                              'erralpha',   plotCfg.errAlpha);
            result.upperY(aIdx, cIdx, lIdx) = upperY;
        end
        % Enriched legend entry: trace label + pool composition.
        baseLabel = requestTraceLabel(parsed.varying, parsed, cIdx, lIdx);
        traceLabels{end+2, 1} = sprintf('%s; N: %d, n: %d, c: %d, tr: %d', ...
            baseLabel, pool.nSubj, pool.nSess, pool.nClust, pool.nTrials);
    end

    % Axes formatting: seconds, integer-only numerals, alignment xline.
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

% Harmonise y-axis across all subplots so they're visually comparable.
% First subplot keeps full y-axis ticks/label; later subplots hide them.
maxY = max(result.upperY(:), [], 'omitnan');
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

%% 06. Save and report (main PSTH).
outDir = fullfile(input.analysis, 'plots', 'fireRate');
if ~exist(outDir, 'dir'), mkdir(outDir); end
fname = encodeFireRateRequest(request);
result.figFile = fullfile(outDir, [fname '.png']);
exportgraphics(fig, result.figFile, 'Resolution', 300);
close(fig);
fprintf('NGL04_fireRate: %s\n', result.figFile);

%% 07. Example waveforms — separate diagnostic figure.
% One panel per alignment, all traces overlaid in their PSTH colour.
% Each trace contributes up to 4 randomly picked cluster mean-waveforms.
% Skipped entirely if no contributing cluster carries a waveform (i.e.,
% opt.getwF was false during NGL01).
anyWF = false;
for k = 1:numel(result.pooled)
    if ~isempty(result.pooled{k}.waveforms), anyWF = true; break, end
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
                pool     = result.pooled{aIdx, cIdx, lIdx};
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
    result.figFile_waveforms = fullfile(outDir, [fname '_waveforms.png']);
    exportgraphics(figWF, result.figFile_waveforms, 'Resolution', 300);
    close(figWF);
    fprintf('NGL04_fireRate: %s\n', result.figFile_waveforms);
else
    warning('NGL04:noWaveforms', ...
        ['No contributing cluster carries a .waveform field; skipping ', ...
         'example-waveforms diagnostic figure. (Set opt.getwF=true in NGL01.)']);
end

% Local helpers (UI-only utilities). Anything related to the request
% parsing / pool building lives under functions/analysis/ and is shared
% with NGL04_PCA. Anything related to common option lookup or array
% utilities stays local because it's plotting-specific.

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