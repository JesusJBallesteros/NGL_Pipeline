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
%   aggregated.mat (or per-subject files) from NGL03_acrossSession.
%
% MULTI-AREA NOTE:
%   v1 supports flat (single-area) aggregated cells. For multi-area
%   aggregation (allspike{x,y}.<Area>.KSLabel), add opt.fireRatePlot.area
%   to pick the area and update Section 04. See task #27.
%
% Last modified <date> (Jesus) - skeleton; bodies TODO (task #27)

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
aggregated = localLoadAggregated(input);
% Required fields in `aggregated`: allspike, allneurons, allcondition.
% Optional: allevents, alltrialdef.

%% 02. Build category sets used by content-based category inference.
catSets = struct();
catSets.alignments  = opt.alignto;                              % from set_default
catSets.conditions  = localCollectConditionFields(aggregated);  % union of condition fieldnames
catSets.labelPool   = localCollectClusterLabels(aggregated);    % union of label-string values
catSets.labelPriority = localFireRatePlotField(opt, 'labelPriority', ...
                          {'HumanLabel','KSLabel','bc_unitType'});

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
parsed = localParseRequest(request, catSets);

nA = numel(parsed.alignment);
nC = numel(parsed.condition);
nL = numel(parsed.label);
fprintf('NGL04_fireRate: %d alignment(s) x %d condition(s) x %d label(s); varying = {%s}.\n', ...
        nA, nC, nL, strjoin(parsed.varying, ', '));

%% 04. Pool per-trial spike cells across all matching (subj, sess, c).
%
% Cartesian product over the three factors. Each (aIdx, cIdx, lIdx)
% cell receives a pooled {1 x N} cell of per-trial spike vectors (ms).
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
            pool = localPoolCells(aggregated, aL, cL, lL, lF, catSets);
            result.pooled{aIdx, cIdx, lIdx} = pool;
            fprintf('  (%s | %s | %s [%s]): %d trials | %d clust | %d sess | %d subj\n', ...
                    aL, cL, lL, lF, pool.nTrials, pool.nClust, pool.nSess, pool.nSubj);
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

fig = figure('Visible','on','Position',[100 100 max(900, 450*nA) 500]);
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
        baseLabel = localTraceLabel(parsed.varying, parsed, cIdx, lIdx);
        traceLabels{end+2, 1} = sprintf('%s; N: %d, n: %d, c: %d, tr: %d', ...
            baseLabel, pool.nSubj, pool.nSess, pool.nClust, pool.nTrials);
    end

    % Axes formatting: seconds, integer-only numerals, alignment xline.
    ax = axHandles(aIdx);
    ax.XTick      = tickBins;
    ax.XTickLabel = tickLab;
    title(localSubplotTitle(parsed, aIdx));
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
outDir = localFireRatePlotField(opt, 'outDir', ...
            fullfile(input.analysis, 'plots', 'NGL04_fireRate'));
if ~exist(outDir, 'dir'), mkdir(outDir); end
fname = localEncodeRequest(request);
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
                    legHandles(end+1) = hLast; %#ok<SAGROW>
                    legNames{end+1}   = sprintf('%s (n=%d shown of %d)', ...
                        localTraceLabel(parsed.varying, parsed, cIdx, lIdx), ...
                        numel(picks), numel(pool.waveforms)); %#ok<SAGROW>
                end
            end
        end
        title(['Waveforms | ' parsed.alignment{aIdx}]);
        xlabel('sample'); ylabel('amplitude');
        if ~isempty(legHandles)
            legend(legHandles, legNames, 'Location', 'best');
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

% ======================================================================
% Local helpers (TODO bodies — see task #27 for the implementation plan).
% ======================================================================

function aggregated = localLoadAggregated(input)
% Load NGL03_acrossSession output. Prefer study-level aggregated.mat;
% fall back to stitching per-subject files into a (Nsubj x maxSess)
% struct of cells.
    studyFile = fullfile(input.analysis, 'aggregated.mat');
    if isfile(studyFile)
        fprintf('NGL04_fireRate: loading study-level %s\n', studyFile);
        aggregated = load(studyFile);
        return
    end
    fprintf('NGL04_fireRate: study-level aggregated.mat not found; stitching per-subject files.\n');
    aggregated = struct();
    for x = 1:input.nsubjects
        sname    = input.subjects(x).name;
        subjFile = fullfile(input.analysis, sname, [sname '_aggregated.mat']);
        if ~isfile(subjFile)
            warning('NGL04:missingPerSubject', ...
                'Per-subject aggregated file missing for %s: %s', sname, subjFile);
            continue
        end
        S  = load(subjFile);
        fn = fieldnames(S);
        for k = 1:numel(fn)
            field = fn{k};
            row   = S.(field);     % {1 x nSess_for_this_subj}
            if ~iscell(row), continue, end
            if ~isfield(aggregated, field), aggregated.(field) = {}; end
            for y = 1:numel(row)
                aggregated.(field){x, y} = row{y};
            end
        end
    end
end

function fields = localCollectConditionFields(aggregated)
% Union of fieldnames seen across aggregated.allcondition cells.
    fields = {};
    if ~isfield(aggregated, 'allcondition'), return; end
    cells = aggregated.allcondition;
    for k = 1:numel(cells)
        c = cells{k};
        if isempty(c) || ~isstruct(c), continue, end
        fields = union(fields, fieldnames(c));
    end
    fields = fields(:);
end

function labelPool = localCollectClusterLabels(aggregated)
% Build a struct keyed by label-field name, each value the cell of
% unique non-empty string values observed across aggregated.allspike.
    labelPool = struct( ...
        'HumanLabel',  {{}}, ...
        'KSLabel',     {{}}, ...
        'bc_unitType', {{}}, ...
        'phyLabel',    {{}});
    if ~isfield(aggregated, 'allspike'), return; end
    cells = aggregated.allspike;
    poolFields = fieldnames(labelPool)';
    for k = 1:numel(cells)
        spk = cells{k};
        if isempty(spk) || ~isstruct(spk), continue, end
        for f = poolFields
            field = f{1};
            if ~isfield(spk, field) || ~iscell(spk.(field)), continue, end
            vals = spk.(field);
            vals = vals(~cellfun(@isempty, vals));
            valsChar = cellfun(@(v) char(string(v)), vals, 'uni', false);
            labelPool.(field) = union(labelPool.(field), valsChar);
        end
    end
end

function v = localFireRatePlotField(opt, fname, dflt)
    if isfield(opt,'fireRatePlot') && isfield(opt.fireRatePlot, fname)
        v = opt.fireRatePlot.(fname);
    else
        v = dflt;
    end
end

function parsed = localParseRequest(request, catSets)
% Categorise each entry and split on 'vs'. Categories: alignment,
% condition, label. Priority for ambiguous tokens: alignment > condition
% > label. For label tokens, the label field is resolved via
% catSets.labelPriority.
    assert(iscell(request) && numel(request) == 3, ...
        'NGL04:badRequest', 'request must be a 3-element cell.');

    entries     = cell(3, 1);
    entryCats   = cell(3, 1);
    entryFields = cell(3, 1);
    for k = 1:3
        s     = strtrim(request{k});
        parts = regexp(s, '\s+vs\s+', 'split', 'ignorecase');
        parts = cellfun(@strtrim, parts, 'uni', false);
        entries{k}     = parts;
        entryCats{k}   = cell(1, numel(parts));
        entryFields{k} = cell(1, numel(parts));
        for p = 1:numel(parts)
            [cat, field] = localCategoriseToken(parts{p}, catSets);
            entryCats{k}{p}   = cat;
            entryFields{k}{p} = field;
        end
        if ~all(strcmp(entryCats{k}, entryCats{k}{1}))
            error('NGL04:badRequest', ...
                'Entry ''%s'' mixes categories (%s); both sides of ''vs'' must be the same kind.', ...
                request{k}, strjoin(entryCats{k}, ', '));
        end
    end

    byCat       = struct('alignment', {{}}, 'condition', {{}}, 'label', {{}});
    labelFields = {};
    for k = 1:3
        cat = entryCats{k}{1};
        if ~isempty(byCat.(cat))
            error('NGL04:badRequest', ...
                ['Two entries both resolve to category ''%s'': ''%s'' and ''%s''. ', ...
                 'Each of the three slots must be a different category.'], ...
                cat, strjoin(byCat.(cat), ' vs '), strjoin(entries{k}, ' vs '));
        end
        byCat.(cat) = entries{k};
        if strcmp(cat, 'label'), labelFields = entryFields{k}; end
    end

    for c = {'alignment','condition','label'}
        if isempty(byCat.(c{1}))
            error('NGL04:badRequest', ...
                'No entry in request resolves to category ''%s''.', c{1});
        end
    end

    parsed = struct();
    parsed.alignment  = byCat.alignment;
    parsed.condition  = byCat.condition;
    parsed.label      = byCat.label;
    parsed.labelField = labelFields;
    parsed.varying    = {};
    catNames = {'alignment','condition','label'};
    for c = catNames
        if numel(byCat.(c{1})) > 1
            parsed.varying{end+1} = c{1};
        end
    end
end

function [cat, field] = localCategoriseToken(token, catSets)
    field = '';
    if any(strcmp(token, catSets.alignments))
        cat = 'alignment'; return
    end
    if any(strcmp(token, catSets.conditions))
        cat = 'condition'; return
    end
    for f = catSets.labelPriority
        fld = f{1};
        if isfield(catSets.labelPool, fld) && ...
                any(strcmp(token, catSets.labelPool.(fld)))
            cat   = 'label';
            field = fld;
            return
        end
    end
    error('NGL04:unknownToken', ...
        ['Cannot categorise ''%s''. It is not in opt.alignto, not a ', ...
         'condition fieldname, and not a known cluster label value ', ...
         '(checked in priority %s).'], token, ...
        strjoin(catSets.labelPriority, ' > '));
end

function pool = localPoolCells(aggregated, alignName, condField, labelValue, labelField, catSets) %#ok<INUSL>
% Walk every aggregated.allspike{x,y} cell. Match clusters whose
% spike.(labelField){c} == labelValue. For each match, pull
% allneurons{x,y}.(alignName){c}, filter by allcondition{x,y}.(condField),
% append to the pooled list and capture the cluster's mean waveform.
% Returns a struct with cells, counts, and per-cluster mean waveforms.
    pool = struct( ...
        'cells',     {{}},  ...
        'nSubj',     0,     ...
        'nSess',     0,     ...
        'nClust',    0,     ...
        'nTrials',   0,     ...
        'waveforms', {{}});

    if ~isfield(aggregated,'allspike') || ~isfield(aggregated,'allneurons') ...
            || ~isfield(aggregated,'allcondition')
        warning('NGL04:missingAggField', ...
            'Aggregated file lacks one of allspike/allneurons/allcondition; pool empty.');
        return
    end

    [nSubj, nSess] = size(aggregated.allspike);
    warnedAlign    = false(nSubj, nSess);
    warnedCond     = false(nSubj, nSess);

    cells       = {};
    waveforms   = {};
    subjFlag    = false(nSubj, 1);
    sessFlag    = false(nSubj, nSess);
    clustCount  = 0;

    for x = 1:nSubj
        for y = 1:nSess
            spk = localSafeIdx(aggregated.allspike,     x, y);
            neu = localSafeIdx(aggregated.allneurons,   x, y);
            cnd = localSafeIdx(aggregated.allcondition, x, y);
            if isempty(spk) || isempty(neu) || isempty(cnd), continue, end
            if ~isstruct(spk) || ~isfield(spk, labelField), continue, end
            if ~isstruct(neu) || ~isfield(neu, alignName)
                if ~warnedAlign(x, y)
                    warning('NGL04:missingAlign', ...
                        'allneurons{%d,%d} missing alignment ''%s''; skipping.', ...
                        x, y, alignName);
                    warnedAlign(x, y) = true;
                end
                continue
            end
            if ~isstruct(cnd) || ~isfield(cnd, condField)
                if ~warnedCond(x, y)
                    warning('NGL04:missingCond', ...
                        'allcondition{%d,%d} missing field ''%s''; skipping.', ...
                        x, y, condField);
                    warnedCond(x, y) = true;
                end
                continue
            end

            mask = logical(cnd.(condField)(:))';
            if ~isfield(spk,'label') || isempty(spk.label), continue, end
            Nclust = numel(spk.label);

            for c = 1:Nclust
                if c > numel(spk.(labelField)), continue, end
                lab = spk.(labelField){c};
                if isempty(lab) || ~strcmp(char(string(lab)), labelValue), continue, end
                if c > numel(neu.(alignName)), continue, end
                trialCells = neu.(alignName){c};
                if isempty(trialCells), continue, end
                if numel(trialCells) ~= numel(mask)
                    warning('NGL04:shapeMismatch', ...
                        'Trial count mismatch at (%d,%d) cluster %d (neurons %d vs condition %d); skipping.', ...
                        x, y, c, numel(trialCells), numel(mask));
                    continue
                end
                filtered = trialCells(mask);
                cells    = [cells; filtered(:)]; %#ok<AGROW>

                clustCount     = clustCount + 1;
                subjFlag(x)    = true;
                sessFlag(x, y) = true;

                if isfield(spk,'waveform') && c <= numel(spk.waveform) && ~isempty(spk.waveform{c})
                    wf = spk.waveform{c};
                    if iscell(wf), wf = wf{1}; end
                    if isnumeric(wf) && ~isempty(wf)
                        if isvector(wf)
                            waveforms{end+1, 1} = wf(:); %#ok<AGROW>
                        else
                            waveforms{end+1, 1} = mean(wf, 2, 'omitnan'); %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end

    pool.cells     = cells;
    pool.nSubj     = sum(subjFlag);
    pool.nSess     = sum(sessFlag(:));
    pool.nClust    = clustCount;
    pool.nTrials   = numel(cells);
    pool.waveforms = waveforms;
end

function v = localSafeIdx(arr, x, y)
    [nx, ny] = size(arr);
    if x > nx || y > ny, v = []; else, v = arr{x, y}; end
end

function s = localSubplotTitle(parsed, aIdx)
% Subplot title shows the alignment for this subplot plus any factor
% that is NOT varying (fixed across the whole figure). Varying non-
% alignment factors (condition, label) are described in the legend.
    parts = {};
    parts{end+1} = parsed.alignment{aIdx};
    if numel(parsed.condition) == 1, parts{end+1} = parsed.condition{1}; end
    if numel(parsed.label)     == 1, parts{end+1} = parsed.label{1};     end
    s = strjoin(parts, ' | ');
end

function s = localTraceLabel(varying, parsed, cIdx, lIdx)
% Legend label for one overlaid trace: only the factors that vary
% within a subplot (condition and/or label, never alignment).
    parts = {};
    if any(strcmp(varying, 'condition'))
        parts{end+1} = parsed.condition{cIdx};
    end
    if any(strcmp(varying, 'label'))
        parts{end+1} = parsed.label{lIdx};
    end
    if isempty(parts)
        if numel(parsed.condition) == 1, parts{end+1} = parsed.condition{1}; end
        if numel(parsed.label)     == 1, parts{end+1} = parsed.label{1};     end
    end
    s = strjoin(parts, ' | ');
end

function s = localEncodeRequest(request)
    raw = cellfun(@(c) strrep(c, ' vs ', '_vs_'), request, 'uni', 0);
    s   = strrep(strjoin(raw, '__'), ' ', '_');
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