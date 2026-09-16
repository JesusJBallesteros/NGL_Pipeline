function makeFixtureExamples(varargin)
%MAKEFIXTUREEXAMPLES  Regenerate the example run log and figures kept in the repo.
%
% PURPOSE:
%   `examples/` holds what a correct run looks like: the console output of
%   building the fixture and checking it, and one figure per analysis. They are
%   committed so anyone can see the expected result without running anything -
%   and so a change that alters it shows up as a diff, in words, in review.
%
%   The examples are produced by this script and never edited by hand. The data
%   they came from is built in a temporary folder and thrown away.
%
% USAGE:
%   makeFixtureExamples                      % all three sessions
%   makeFixtureExamples('sessions', {'19850214'})
%   makeFixtureExamples('root', 'D:\scratch\MRX')
%
% INPUTS (name/value):
%   'root'     where to build the throwaway study (default: a temp folder)
%   'sessions' which sessions to build (default: all three)
%   'keep'     true to leave the generated data behind (default false)
%
% WRITES (under tests/fixture/examples/):
%   expected_output.txt   the run log: generation, then every check
%   example_contrast.png  correct vs incorrect, NCL, aligned to stimOn2
%   example_csd.png       current source density down shank 1
%   example_tagging.png   tagging spectrum of the 1.3 Hz stream
%   example_psth.png      one PSTH per task, driven unit against a quiet one
%   example_pecks.png     peck detection: dms response, arena screen and feeder
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).

    p = inputParser;
    p.addParameter('root',     '');
    p.addParameter('sessions', {'19850214', '19860704', '19870321'});
    p.addParameter('keep',     false);
    p.parse(varargin{:});
    a = p.Results;

    here = fileparts(mfilename('fullpath'));
    outDir = fullfile(here, 'examples');
    if ~isfolder(outDir), mkdir(outDir); end
    root = a.root;
    if isempty(root), root = fullfile(tempdir, 'MRXFIX_examples'); end
    if isfolder(root), rmdir(root, 's'); end

    %% 1. Build and check, capturing everything that is printed.
    log = {};
    log{end+1} = sprintf(['Example run of the MRX fixture.\n' ...
        'Produced by tests/fixture/makeFixtureExamples.m - do not edit by hand.\n' ...
        'MATLAB %s\n'], version('-release'));
    log{end+1} = localBanner('GENERATE');
    log{end+1} = evalc(sprintf('makeFakeStudy(root, ''sessions'', %s);', ...
                               localCellLiteral(a.sessions)));
    for k = 1:numel(a.sessions)
        log{end+1} = localBanner(['CHECK ' a.sessions{k}]); %#ok<AGROW>
        log{end+1} = evalc(sprintf( ...
            'checkFixture(root, ''session'', ''%s'', ''strict'', false);', ...
            a.sessions{k})); %#ok<AGROW>
    end

    %% 2. One figure per analysis, from the first session.
    log{end+1} = localBanner('FIGURES');
    files = {};
    try
        log{end+1} = evalc('files = localFigures(root, a.sessions{1}, outDir);');
        for f = files(:)'
            log{end+1} = sprintf('wrote %s\n', f{1}); %#ok<AGROW>
        end
    catch ME
        log{end+1} = sprintf('FIGURES FAILED: %s\n', ME.message);
    end

    txt = fullfile(outDir, 'expected_output.txt');
    fid = fopen(txt, 'w');
    fprintf(fid, '%s', localStripFormatting(strjoin(log, newline)));
    fclose(fid);
    fprintf('makeFixtureExamples: wrote %s\n', txt);

    if ~a.keep && isfolder(root), rmdir(root, 's'); end
end

% ======================= figures =======================
function files = localFigures(root, session, outDir)
% The three analyses that have a picture worth keeping. Deliberately the same
% calls a user would make - if the API changes, this stops working, which is
% the point of keeping it runnable rather than pasting screenshots.
    S = localSession(root, session);
    opt = S.opt;
    opt.lfp.plot.Resolution = 150;      % keep the committed PNGs small
    files = {};

    %% (a) Event-centered power contrast: correct vs incorrect, NCL, stimOn2.
    FT = localLoad(fullfile(S.trial, [session '_stimOn2.mat']));
    cond = localSubset(S.condition, FT.trialinfo);
    isDMS = ismember(FT.trialinfo(:), localTrialsFor(S, 'contrast'));
    TFR = computeTrialparsedTFR(FT, cond, struct(), opt, 'stimOn2');
    labels = FT.label(strcmp(FT.chanArea, 'NCL'));
    band = ft_selectdata(struct('channel', {labels}), TFR{1});
    spec = parseTrialContrast(struct('A', cond.correct(:) > 0 & isDMS, ...
                                     'B', cond.incorrect(:) > 0 & isDMS, ...
                                     'labelA', 'correct', 'labelB', 'incorrect'), ...
                              [], size(band.powspctrm, 1));
    res = computeTFRcontrast(band, spec, opt, 'area', 'NCL', 'align', 'stimOn2');
    fig = plotTFRcontrast(res, opt, 'save', false);
    files{end+1} = localExport(fig, fullfile(outDir, 'example_contrast.png'), opt);

    %% (b) Current source density down shank 1, aligned to stimOn1.
    FT1 = localLoad(fullfile(S.trial, [session '_stimOn1.mat']));
    FT1 = localWindow(FT1, [-0.2 0.4]);
    shanks = lfpChannelGeometry(FT1, opt.fixtureInput, opt, ...
                                'chanMap', S.truth.chanMapFile);
    sh = shanks([shanks.shank] == S.truth.lfp.csdShank);
    use = ismember(FT1.trialinfo, localTrialsFor(S, 'csd'));
    csd = computeCSD(FT1, sh, opt, 'trials', use(:));
    fig = plotCSD(csd, opt, 'align', 'stimOn1', 'save', false);
    files{end+1} = localExport(fig, fullfile(outDir, 'example_csd.png'), opt);

    %% (c) Tagging spectrum of the 1.3 Hz stream.
    FTc = localLoad(fullfile(S.preproc, [session '_FTcont.mat']));
    blk = S.blocks(strcmp({S.blocks.label}, 'nft'));
    trials = [S.truth.blocks(strcmp({S.truth.blocks.task}, 'nft')).trials];
    in = [trials.rate] == 1.3;
    win = [max(blk(1).tStart, min([trials(in).tITI])), max([trials(in).tEnd])];
    epochs = nftEpochs(FTc, win, 1.3);
    spc = computeTaggingSpectrum(epochs, opt);
    resp = taggingResponse(spc, 1.3, opt);
    fig = plotTaggingSpectrum(resp, opt, 'save', false);
    files{end+1} = localExport(fig, fullfile(outDir, 'example_tagging.png'), opt);

    %% (d) Spikes: one PSTH per task, through sort2trials + calcFireRate.
    files{end+1} = localPSTHfigure(S, opt, outDir);

    %% (e) Pecks: what the accelerometer shows at a screen interaction.
    files{end+1} = localPeckFigure(S, outDir);
end

% ======================= spikes =======================
function f = localPSTHfigure(S, opt, outDir)
% One column per task, aligned to stimOn1: the driven unit against a unit with
% no stimulus job, over the same trials. Three tasks, one unit, one answer -
% and the flat line underneath is what says the answer means anything.
    spike = localLoad(fullfile(S.spike, 'spike.mat'));
    neurons = sort2trials(spike, S.trialdef, opt);
    ids = cellfun(@str2double, spike.label);
    uD = find(ids == S.truth.spikes.drivenUnits(1), 1);
    uQ = find(~ismember(ids, S.truth.spikes.drivenUnits) & ...
              ids ~= S.truth.spikes.rewardUnit & ...
              ids ~= S.truth.spikes.taggedUnit, 1);
    param = struct('binSize', 50, 'stepSz', 10, 'interval', [-500 500], ...
                   'smpRate', 1000, 'baseline', 500, 'plot', false);

    tasks = {'dms', 'arena', 'nft'};
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 640]);
    tl = tiledlayout(fig, 2, numel(tasks), 'TileSpacing', 'compact', ...
                     'Padding', 'compact');
    title(tl, sprintf('%s / %s   spikes aligned to stimOn1', S.truth.subject, ...
          S.truth.session), 'FontWeight', 'bold', 'Interpreter', 'none');

    raster = cell(1, numel(tasks)); rates = cell(1, numel(tasks));
    for k = 1:numel(tasks)
        idx = localTrialsInBlock(S, tasks{k});
        raster{k} = neurons.stimOn1{uD}(idx);
        [tAxis, rD] = localRate(neurons.stimOn1{uD}, idx, param);
        [~,     rQ] = localRate(neurons.stimOn1{uQ}, idx, param);
        rates{k} = [rD; rQ];

        % Raster: every trial of the block, one row, cropped to the window.
        ax = nexttile(tl, k); hold(ax, 'on');
        for i = 1:numel(raster{k})
            v = raster{k}{i};
            v = v(v >= param.interval(1) & v <= param.interval(2));
            if isempty(v), continue; end
            plot(ax, v, i * ones(size(v)), '.', 'MarkerSize', 4, 'Color', [.15 .15 .15]);
        end
        xline(ax, 0, 'k-');
        xlim(ax, param.interval); ylim(ax, [0 numel(raster{k}) + 1]);
        set(ax, 'TickDir', 'out', 'Box', 'off');
        title(ax, sprintf('%s  (%d trials)', tasks{k}, numel(idx)), 'Interpreter', 'none');
        if k == 1, ylabel(ax, 'trial'); end
    end

    yl = [0, max(cellfun(@(r) max(r(:)), rates)) * 1.15];
    for k = 1:numel(tasks)
        ax = nexttile(tl, numel(tasks) + k); hold(ax, 'on');
        plot(ax, tAxis, rates{k}(1, :), 'LineWidth', 1.6, 'Color', [0.16 0.35 0.60]);
        plot(ax, tAxis, rates{k}(2, :), 'LineWidth', 1.2, 'Color', [0.65 0.65 0.65]);
        xline(ax, 0, 'k-');
        xlim(ax, param.interval); ylim(ax, yl);
        set(ax, 'TickDir', 'out', 'Box', 'off');
        xlabel(ax, 'time from stimOn1 (ms)');
        if k == 1
            ylabel(ax, 'spikes/s');
            legend(ax, {sprintf('unit %d (driven)', ids(uD)), ...
                        sprintf('unit %d (not driven)', ids(uQ))}, ...
                   'Location', 'northwest', 'Box', 'off');
        end
    end
    f = fullfile(outDir, 'example_psth.png');
    exportgraphics(fig, f, 'Resolution', 150);
    close(fig);
end

function [t, rate] = localRate(trialSpikes, idx, param)
    sel = trialSpikes(idx);
    sel(cellfun(@isempty, sel)) = {NaN};
    fr = calcFireRate(sel, param, struct());
    if iscell(fr), fr = cell2mat(fr); end
    rate = mean(fr, 1, 'omitnan');
    t = param.interval(1) + param.binSize/2 + (0:size(fr, 2)-1) * param.stepSz;
end

function idx = localTrialsInBlock(S, label)
    k = find(strcmp({S.blocks.label}, label), 1);
    if isempty(k), idx = []; else, idx = S.blocks(k).trialIdx; end
end

% ======================= pecks =======================
function f = localPeckFigure(S, outDir)
% What a peck looks like in the accelerometer, in both tasks that have one.
% Top row: the jerk magnitude around a single interaction, with the detected
% time marked. Bottom row: every peck of that task, averaged - the perched
% bird's peck against the same movement made while walking.
    L = load(fullfile(S.preproc, 'MotionData_raw.mat'), '-mat');
    raw = L.raw; fs = raw.fs;
    acc = [raw.acc.X(:), raw.acc.Y(:), raw.acc.Z(:)];
    [b, a] = butter(3, 20 / (fs/2), 'high');
    aHP = filtfilt(b, a, acc);
    jerk = zeros(size(aHP));
    jerk(2:end-1, :) = (aHP(3:end, :) - aHP(1:end-2, :)) / (2 / fs);
    mag = sqrt(sum(jerk.^2, 2));
    thr = median(mag) + 15 * mad(mag, 1);
    detected = localThreshold(mag, thr, round(0.1 * fs)) / fs;

    % Pecks by task, from the block table: in dms the single response peck,
    % in the arena the screen interaction and then the feeder on the way back.
    trials = [S.truth.blocks.trials];
    tasks = {'dms', 'arena'};
    what = {'response peck', 'screen peck'};
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 620]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf(['%s / %s   peck detection on the jerk magnitude ' ...
        '(threshold = median + 15 MAD)'], S.truth.subject, S.truth.session), ...
        'FontWeight', 'bold', 'Interpreter', 'none');

    for k = 1:numel(tasks)
        isT = strcmp({trials.task}, tasks{k});
        pk = arrayfun(@(t) localFirstPeck(t), trials(isT), 'UniformOutput', false);
        pk = [pk{:}];
        pk = pk(pk > 1 & pk < numel(mag)/fs - 1);
        assert(~isempty(pk), 'no pecks in the %s block', tasks{k});

        % (top) one interaction, in full
        ax = nexttile(tl, k); hold(ax, 'on');
        c = pk(min(3, numel(pk)));
        w = round(0.4 * fs);
        i0 = round(c * fs) - w; i1 = round(c * fs) + w;
        tt = ((i0:i1) - round(c*fs)) / fs * 1000;
        plot(ax, tt, mag(i0:i1), 'Color', [0.16 0.35 0.60], 'LineWidth', 1.1);
        yline(ax, thr, '--', 'threshold', 'Color', [0.55 0.55 0.55]);
        d = detected(detected > c - 0.4 & detected < c + 0.4);
        plot(ax, (d - c) * 1000, repmat(max(mag(i0:i1)) * 1.05, size(d)), ...
             'v', 'MarkerFaceColor', [0.75 0.22 0.17], 'MarkerEdgeColor', 'none');
        xlim(ax, [-400 400]); set(ax, 'TickDir', 'out', 'Box', 'off');
        title(ax, sprintf('%s: one %s', tasks{k}, what{k}), 'Interpreter', 'none');
        xlabel(ax, 'time from peck (ms)');
        if k == 1, ylabel(ax, 'jerk magnitude'); end

        % (bottom) every peck of this task, averaged
        ax = nexttile(tl, 2 + k); hold(ax, 'on');
        w2 = round(0.25 * fs);
        M = zeros(numel(pk), 2*w2 + 1);
        for i = 1:numel(pk)
            j0 = round(pk(i) * fs) - w2;
            M(i, :) = mag(j0:j0 + 2*w2);
        end
        tt2 = (-w2:w2) / fs * 1000;
        plot(ax, tt2, mean(M, 1), 'Color', [0.16 0.35 0.60], 'LineWidth', 1.6);
        plot(ax, tt2, prctile(M, [10 90]), ':', 'Color', [0.16 0.35 0.60]);
        yline(ax, thr, '--', 'Color', [0.55 0.55 0.55]);
        xlim(ax, [-250 250]); set(ax, 'TickDir', 'out', 'Box', 'off');
        hit = mean(arrayfun(@(x) any(abs(detected - x) <= 0.05), pk));
        title(ax, sprintf('%d pecks averaged - %.0f%% detected', numel(pk), 100*hit), ...
              'Interpreter', 'none');
        xlabel(ax, 'time from peck (ms)');
        if k == 1, ylabel(ax, 'jerk magnitude'); end
    end
    f = fullfile(outDir, 'example_pecks.png');
    exportgraphics(fig, f, 'Resolution', 150);
    close(fig);
end

function p = localFirstPeck(tr)
    if isempty(tr.pecks), p = []; else, p = tr.pecks(1); end
end

function idx = localThreshold(x, thr, refr)
% One event per excursion above threshold, then a refractory sweep.
    above = x(:) > thr;
    d = diff([false; above; false]);
    st = find(d == 1); sp = find(d == -1) - 1;
    idx = zeros(numel(st), 1);
    for k = 1:numel(st)
        [~, rel] = max(x(st(k):sp(k)));
        idx(k) = st(k) + rel - 1;
    end
    keep = true(size(idx)); last = -Inf;
    for k = 1:numel(idx)
        if idx(k) - last < refr, keep(k) = false; else, last = idx(k); end
    end
    idx = idx(keep);
end

function f = localExport(fig, f, opt)
    exportgraphics(fig, f, 'Resolution', opt.lfp.plot.Resolution);
    close(fig);
end

% ======================= small helpers =======================
function S = localSession(root, session)
% The same loading checkFixture does, kept separate so the two cannot drift
% into disagreeing about where a file lives.
    S.preproc = fullfile(root, 'data', 'preprocessing', 'MRX', session);
    S.trial   = fullfile(root, 'data', 'trialSorted',   'MRX', session);
    S.spike   = fullfile(root, 'data', 'spikeSorted',   'MRX', session);
    S.truth     = load(fullfile(S.preproc, 'fixture_truth.mat'));
    S.condition = localLoad(fullfile(S.trial, 'condition.mat'));
    S.blocks    = localLoad(fullfile(S.trial, 'blocks.mat'));
    S.trialdef  = localLoad(fullfile(S.trial, 'trialdef.mat'));
    opt = generateDefaultsFromSchema(optSchema());
    opt.SavFileName = session;
    opt.trialSorted = S.trial;
    opt.analysis    = fullfile(root, 'data', 'analysis', 'MRX', session);
    opt.KSchanMapFile = S.truth.chanMapFile;
    opt.TFRmethod   = 'wavelet';
    opt.freqInterest = {8:1:40};
    % A second wider on each side than the effect needs; the fixed trial-end
    % event is what makes every trial reach that far.
    opt.toi = [-2 1.8]; opt.timeResol = 0.05;
    opt.alignto = {'itiOn', 'stimOn1', 'stimOn2'};   % as trialdef was built
    opt.fixtureInput = struct('analysisCode', '', 'areaMap', ...
                              struct('chanMapPath', S.truth.chanMapFile));
    S.opt = opt;
end

function v = localLoad(file)
    L = load(file, '-mat');
    f = fieldnames(L);
    v = L.(f{1});
    if isstruct(v) && isfield(v, 'FT_data'), v = v.FT_data; end
end

function idx = localTrialsFor(S, analysis)
    idx = [];
    for k = 1:numel(S.blocks)
        label = S.blocks(k).label;
        if isfield(S.truth.taskAnalyses, label) && ...
                any(strcmp(S.truth.taskAnalyses.(label), analysis))
            idx = [idx, S.blocks(k).trialIdx]; %#ok<AGROW>
        end
    end
end

function cond = localSubset(cond, idx)
    f = fieldnames(cond);
    for k = 1:numel(f)
        v = cond.(f{k});
        if isnumeric(v) && numel(v) >= max(idx), cond.(f{k}) = v(idx); end
    end
end

function FT = localWindow(FT, win)
% Cut to a common window and keep only the trials that span it - a CSD averages
% across trials in time, so they have to be squared off first.
    FT = ft_redefinetrial(struct('toilim', win), FT);
    want = floor(diff(win) * FT.fsample);
    n = cellfun(@numel, FT.time);
    keep = n >= want;
    cut = min(n(keep));
    for k = find(keep(:)')
        FT.trial{k} = FT.trial{k}(:, 1:cut);
        FT.time{k}  = FT.time{k}(1:cut);
    end
    FT.trial = FT.trial(keep);
    FT.time  = FT.time(keep);
    FT.trialinfo = FT.trialinfo(keep);
end

function s = localBanner(what)
    s = sprintf('\n%s\n== %s\n%s\n', repmat('=', 1, 72), what, repmat('=', 1, 72));
end

function s = localCellLiteral(c)
    s = ['{''' strjoin(c, ''', ''') '''}'];
end

function s = localStripFormatting(s)
% Make the captured output readable as a text file. Three things get in the
% way: carriage returns (FieldTrip overwrites its progress lines), backspaces
% (MATLAB marks a warning with one), and FieldTrip's own running commentary,
% which is per-call, changes with its version, and would bury the lines this
% file exists to show. What the PIPELINE prints stays, its warnings included -
% a warning that belongs to the pipeline is exactly what a reference run
% should show.
    s = regexprep(s, '[^\n]*\r', '');
    s = regexprep(s, '\r', '');
    s = strrep(s, char(8), '');
    noise = {'^the call to', '^the input is', '^processing', '^computing', ...
             '^found \d', '^using a ', '^using "', '^selecting ', '^averaging ', ...
             '^total number of', '^number of ', '^estimated time', ...
             '^evaluating ', '^the returned probabilities', '^\s+In ', ...
             '^repairing ', '^reading ', '^converting ', '^removing ', ...
             '^resampling', '^original sampling rate', '^new sampling rate', ...
             '^constructing randomized', '^trial \d+, frequency', ...
             '^Starting parallel pool', '^Connected to parallel pool', ...
             '^time axes\]'};   % the tail of a wrapped FieldTrip warning
    % FieldTrip's warnings, which say nothing about this data.
    ftWarn = ['^\[Warning: (adding|not all trials|the data does not contain|' ...
              'reconstructing sampleinfo|correcting numerical|' ...
              'could not determine dimord)'];

    lines = strsplit(s, newline);
    keep = true(size(lines));
    k = 1;
    while k <= numel(lines)
        % A warning is a block: it runs until the line that closes the bracket,
        % and the decision to keep it belongs to its FIRST line.
        if ~isempty(regexp(lines{k}, ftWarn, 'once'))
            % Drop the warning and the struct dump it prints, line by line:
            % only the dump's own shapes ("  field: value", a bare closing
            % bracket, blanks) go. evalc interleaves the two output streams,
            % so a line of real output can land in the middle of a warning -
            % deleting a whole span between brackets would take it with it.
            keep(k) = false;
            k = k + 1;
            while k <= numel(lines)
                isDump = ~isempty(regexp(lines{k}, '^\s+\S+:\s', 'once')) || ...
                         ~isempty(regexp(lines{k}, '^\s*$', 'once'));
                isEnd  = ~isempty(regexp(lines{k}, '^\]\s*$', 'once'));
                if ~(isDump || isEnd), break; end
                keep(k) = false;
                k = k + 1;
                if isEnd, break; end
            end
            continue
        end
        for n = 1:numel(noise)
            if ~isempty(regexp(lines{k}, noise{n}, 'once'))
                keep(k) = false; break
            end
        end
        k = k + 1;
    end
    lines = lines(keep);
    % Collapse the blank runs the removals leave behind.
    blank = cellfun(@(l) isempty(strtrim(l)), lines);
    lines = lines(~(blank & [blank(2:end), false]));
    s = strjoin(lines, newline);
end
