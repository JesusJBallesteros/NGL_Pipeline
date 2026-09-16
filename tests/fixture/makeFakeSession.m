function truth = makeFakeSession(root, subject, session, varargin)
%MAKEFAKESESSION  Build one synthetic session with known answers, in the real layout.
%
% PURPOSE:
%   A test dataset whose results are known in advance. Every analysis in the
%   pipeline can then be checked against what was planted, instead of against
%   what looks plausible - and everyone checks against the same thing, because
%   the generator is seeded and its output is byte-identical on any machine.
%
%   The repository holds this generator, not the data it makes: a session is
%   tens of megabytes, which does not belong in git, and a fixed seed gives
%   the same bytes anyway.
%
% USAGE:
%   truth = makeFakeSession('D:\TESTSTUDY', 'MRX', '19850214')
%   truth = makeFakeSession(root, 'MRX', '19850214', 'tasks', {'dms','arena'})
%   truth = makeFakeSession(root, 'MRX', '19850214', 'length', 'full')
%
% INPUTS:
%   root    - study folder; the standard tree is created under <root>\data\.
%   subject - subject ID. Use MRX for the shipped fixture.
%   session - session name, normally a date. The fixture uses PRE-1990 dates
%             so a test session can never be mistaken for a real recording.
%   Name/value pairs:
%     'tasks'  which blocks, in order (default all three):
%                'dms'   delay-match-to-sample, perched, pecking a screen
%                'arena' freely moving, 3 of 6 screens in a hexagonal arena
%                'nft'   passive rapid stimulus stream (frequency tagging)
%     'length' 'short' (default, ~8 min) or 'full' (~20 min)
%     'seed'   random seed (default derived from the session name, so each
%              session differs but is reproducible)
%     'fs'     LFP sample rate (default 937.5 Hz, the pipeline's auto value)
%
% WRITES (the layout set_default builds):
%   data\preprocessing\<subj>\<sess>\  <sess>_FTcont.mat, EventRecord.mat,
%                                      MotionData_raw.mat
%   data\trialSorted\<subj>\<sess>\    events.mat, trialdef.mat, condition.mat,
%                                      <sess>_<align>.mat (one per alignment)
%   data\spikeSorted\<subj>\<sess>\    spike.mat
%   data\analysis\<subj>\<sess>\       (empty; analyses write here)
%   ...\<sess>\fixture_truth.mat       what was planted, for the tests
%
% OUTPUT:
%   truth - the same struct saved as fixture_truth.mat. Its fields say what
%           every analysis should find; see makeFakeStudy for the summary.
%
% NOTES:
%   * Timing is deliberately irregular. Response latencies are drawn per
%     trial with a floor (nobody reacts instantly) and a ceiling (the response
%     window), so trial lengths vary the way real ones do and code that
%     assumes a fixed trial length fails here - which is the point.
%   * The recording carries realistic trouble: one dead channel, one noisy
%     channel, occasional movement artifacts, and a unit whose amplitude
%     drifts across the session. A fixture without them would pass tests that
%     real data fails.
%   * Spikes are written as spike.mat directly, in the shape NGL02_postPhy
%     produces after Kilosort and Phy. Running a sorter on synthetic data
%     would test the sorter, not the pipeline.
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).

    p = inputParser;
    p.addParameter('tasks',  {'dms', 'arena', 'nft'});
    p.addParameter('length', 'short');
    p.addParameter('seed',   []);
    p.addParameter('fs',     30000 / 32);      % 937.5 Hz, the pipeline default
    p.parse(varargin{:});
    a = p.Results;
    if ischar(a.tasks), a.tasks = {a.tasks}; end
    seed = a.seed;
    if isempty(seed), seed = mod(str2double(session), 2^31 - 1); end
    rng(seed, 'twister');

    cfg = fixtureConfig(a.length, a.fs);
    cfg.subject = subject; cfg.session = session; cfg.seed = seed;

    %% Geometry: the real 32-channel, 2-shank map the lab uses.
    [geom, cfg.chanMapFile] = localGeometry();
    cfg.nChan = numel(geom.label);

    %% 1. Trial structure per block, then the event record they imply.
    blkCell = cell(1, numel(a.tasks));
    tCursor = cfg.leadIn;
    for b = 1:numel(a.tasks)
        blkCell{b} = localBuildBlock(a.tasks{b}, tCursor, cfg);
        tCursor = blkCell{b}.tEnd + cfg.blockGap;
    end
    blocks = [blkCell{:}];
    duration = tCursor + cfg.leadIn;

    [EventRecord, eventList] = localEventRecord(blocks, cfg);

    %% 2-4. The signals, and the record of what was planted in them.
    [lfp, spike, motion, sigTruth] = makeFixtureSignals(blocks, geom, duration, cfg);

    %% 5. Assemble and write.
    paths = localPaths(root, subject, session);
    FT_data = makeFixtureFT(lfp, geom, cfg);
    save(fullfile(paths.preproc, [session '_FTcont.mat']), 'FT_data', '-v7.3');
    save(fullfile(paths.preproc, 'EventRecord.mat'), 'EventRecord', '-v7.3');
    save(fullfile(paths.preproc, 'MotionData_raw.mat'), '-struct', 'motion', '-v7.3');
    % (-struct writes the single variable `raw`, the name GetMotionSensors uses)

    [events, trialdef, condition] = localTrialTables(blocks, eventList, cfg);
    save(fullfile(paths.trial, 'events.mat'), 'events', '-v7.3');
    save(fullfile(paths.trial, 'trialdef.mat'), 'trialdef', '-v7.3');
    save(fullfile(paths.trial, 'condition.mat'), 'condition', '-v7.3');
    aligns = fieldnames(events);
    for k = 1:numel(aligns)
        S = struct('FT_data', makeFixtureFT(lfp, geom, cfg, trialdef, aligns{k}));
        save(fullfile(paths.trial, [session '_' aligns{k} '.mat']), '-struct', 'S', '-v7.3');
        clear S
    end
    save(fullfile(paths.spike, 'spike.mat'), 'spike', '-v7.3');

    truth = struct('subject', subject, 'session', session, 'seed', seed, ...
                   'fs', cfg.fs, 'duration', duration, 'blocks', blocks, ...
                   'chanMapFile', cfg.chanMapFile, 'lfp', sigTruth.lfp, ...
                   'spikes', sigTruth.spikes, 'imu', sigTruth.imu, 'nTrials', ...
                   sum(arrayfun(@(b) numel(b.trials), blocks)));
    save(fullfile(paths.preproc, 'fixture_truth.mat'), '-struct', 'truth');

    fprintf(['makeFakeSession: %s / %s - %.0f s, %d channels, %d trials in %d ' ...
             'block(s), %d units\n'], subject, session, duration, cfg.nChan, ...
            truth.nTrials, numel(blocks), numel(spike.label));
end

% ======================= configuration =======================
function cfg = fixtureConfig(len, fs)
    cfg = struct();
    cfg.fs = fs;
    cfg.leadIn = 5;                 % quiet seconds at each end
    cfg.blockGap = 8;
    switch lower(len)
        case 'short', cfg.nDMS = 40; cfg.nArena = 24; cfg.nftSeconds = 120;
        case 'full',  cfg.nDMS = 90; cfg.nArena = 60; cfg.nftSeconds = 300;
        otherwise, error('makeFakeSession:length', 'length must be ''short'' or ''full''.');
    end
    % Reserved event codes (eventDefinitions.m; do not change).
    cfg.code = struct('itiOn', 0, 'stimOn1', 1, 'stimOn2', 2, 'bhv', 3, ...
                      'end1', 4, 'rwd', 7, 'preIni', 8, 'end2', 10, ...
                      'pun', 11, 'end3', 15);
    % Project codes for the tagging stimuli (>= 16, as the convention requires).
    cfg.nftCodes = 8001:8006;
    cfg.bands = struct('theta', [4 8], 'beta', [15 30]);
end

function [geom, mapFile] = localGeometry()
    here = fileparts(mfilename('fullpath'));
    mapFile = fullfile(fileparts(fileparts(here)), 'channelmaps', ...
                       'chanMap_ATLAS-E32+R-50-S2-L10-200NT_INTAN.mat');
    assert(isfile(mapFile), 'makeFakeSession:noMap', ...
        'channel map not found: %s', mapFile);
    m = load(mapFile);
    geom.ycoords = m.ycoords(:);
    geom.kcoords = m.kcoords(:);
    n = numel(geom.ycoords);
    geom.label = arrayfun(@(k) sprintf('LFP_%02d', k), (1:n)', 'UniformOutput', false);
    geom.area = repmat({'NCL'}, n, 1);
    geom.area(geom.kcoords == 2) = {'STR'};
end

function paths = localPaths(root, subject, session)
    paths.preproc = fullfile(root, 'data', 'preprocessing', subject, session);
    paths.trial   = fullfile(root, 'data', 'trialSorted',   subject, session);
    paths.spike   = fullfile(root, 'data', 'spikeSorted',   subject, session);
    paths.analysis= fullfile(root, 'data', 'analysis',      subject, session);
    paths.raw     = fullfile(root, 'data', 'raw', subject, session);
    f = fieldnames(paths);
    for k = 1:numel(f)
        if ~isfolder(paths.(f{k})), mkdir(paths.(f{k})); end
    end
end

% ======================= trial structure =======================
function blk = localBuildBlock(task, t0, cfg)
% Each task has its own trial anatomy and its own way of running late, but
% every trial carries the SAME fields (localNewTrial): blocks are concatenated
% downstream, and a struct array cannot hold trials of differing shape.
    trials = {};
    t = t0;
    switch lower(task)
        case 'dms'
            for k = 1:cfg.nDMS
                tr = localNewTrial('dms');
                tr.iti    = 0.8 + 0.4 * rand;              % variable baseline
                tr.sample = 0.35;
                tr.delay  = 1.0 + 0.5 * rand;              % memory delay
                tr.test   = 0.35;
                % The outcome comes first: it decides whether there is a
                % response at all, and so how long the trial runs.
                u = rand;
                if     u < 0.62, tr.outcome = 'correct';
                elseif u < 0.85, tr.outcome = 'incorrect';
                elseif u < 0.95, tr.outcome = 'omission';
                else,            tr.outcome = 'aborted';
                end
                switch tr.outcome
                    case 'omission'
                        tr.rt = NaN; tr.resolve = 2.0;      % the response window
                    case 'aborted'
                        tr.rt = 0.05 + 0.1 * rand;          % pecked too early
                        tr.resolve = tr.rt;
                    otherwise
                        % Floor: nobody reacts instantly. Ceiling: the window.
                        tr.rt = min(0.18 + 0.35 * abs(randn), 2.0);
                        tr.resolve = tr.rt;
                end
                tr.tITI   = t;
                tr.tStim1 = t + tr.iti;
                tr.tStim2 = tr.tStim1 + tr.sample + tr.delay;
                tr.tResp  = tr.tStim2 + tr.test + tr.resolve;
                tr.tEnd   = tr.tResp + 0.4 + 0.3 * rand;
                if ~strcmp(tr.outcome, 'omission'), tr.pecks = tr.tResp; end
                t = tr.tEnd + 0.2;
                trials{end+1} = tr; %#ok<AGROW>
            end

        case 'arena'
            screens = [1 3 5];                             % 3 of 6, alternating
            for k = 1:cfg.nArena
                tr = localNewTrial('arena');
                tr.iti    = 1.0 + 0.6 * rand;
                tr.screen = screens(mod(k - 1, numel(screens)) + 1);
                tr.travel = 1.6 + 1.4 * rand;              % walking takes time
                tr.sample = 0.3;
                u = rand;
                if     u < 0.70, tr.outcome = 'correct';
                elseif u < 0.88, tr.outcome = 'incorrect';
                else,            tr.outcome = 'omission';
                end
                if strcmp(tr.outcome, 'omission')
                    tr.rt = NaN; tr.resolve = 4.0;
                else
                    tr.rt = min(0.2 + 0.5 * abs(randn), 3.0);
                    tr.resolve = tr.rt;
                end
                tr.tITI    = t;
                tr.tStim1  = t + tr.iti;
                tr.tStim2  = tr.tStim1 + tr.sample + tr.travel;   % reaches screen
                tr.tResp   = tr.tStim2 + tr.resolve;
                tr.tReturn = tr.tResp + 1.4 + 1.0 * rand;         % walks back
                tr.tEnd    = tr.tReturn;
                if ~strcmp(tr.outcome, 'omission')
                    tr.pecks = [tr.tResp, tr.tReturn - 0.2];      % screen, feeder
                end
                t = tr.tEnd + 0.3;
                trials{end+1} = tr; %#ok<AGROW>
            end

        case 'nft'
            % Passive: no response, a fast stream of short stimuli. Two streams,
            % as the real design has - 1.3 Hz, then 2.6 Hz paired into syllables
            % so structure sits at half the component rate.
            rates = [1.3 2.6];
            per = cfg.nftSeconds / numel(rates);
            for r = 1:numel(rates)
                period = 1 / rates(r);
                for k = 1:floor(per / period)
                    tr = localNewTrial('nft');
                    tr.outcome = 'passive';
                    tr.rate    = rates(r);
                    tr.paired  = rates(r) > 2;             % the 2.6 Hz block pairs
                    tr.stimIdx = mod(k - 1, 6) + 1;
                    tr.code    = cfg.nftCodes(tr.stimIdx);
                    tr.tITI    = t;
                    tr.tStim1  = t + period - 0.23;        % 230 ms stimulus
                    tr.tEnd    = t + period;
                    t = tr.tEnd;
                    trials{end+1} = tr; %#ok<AGROW>
                end
            end
        otherwise
            error('makeFakeSession:task', ...
                'unknown task ''%s''; use dms, arena or nft.', task);
    end
    blk = struct('task', task, 'tStart', t0, 'tEnd', t, 'trials', [trials{:}]);
end

function tr = localNewTrial(task)
% The canonical trial. Every field exists for every task; the ones a task does
% not use stay NaN or empty, so trials concatenate and a consumer can test a
% field rather than the task name.
    tr = struct('task', task, 'outcome', '', 'iti', NaN, 'sample', NaN, ...
                'delay', NaN, 'test', NaN, 'travel', NaN, 'screen', NaN, ...
                'rate', NaN, 'paired', false, 'stimIdx', NaN, 'code', [], ...
                'rt', NaN, 'resolve', NaN, 'tITI', NaN, 'tStim1', NaN, ...
                'tStim2', NaN, 'tResp', NaN, 'tReturn', NaN, 'tEnd', NaN, ...
                'pecks', []);
end

% ======================= events =======================
function [EventRecord, eventList] = localEventRecord(blocks, cfg)
% The digital-input record NGL01 would extract: one row per code transition.
    times = []; codes = [];
    for b = 1:numel(blocks)
        for k = 1:numel(blocks(b).trials)
            tr = blocks(b).trials(k);
            times(end+1) = tr.tITI - 0.03;   codes(end+1) = cfg.code.preIni; %#ok<AGROW>
            times(end+1) = tr.tITI;          codes(end+1) = cfg.code.itiOn;  %#ok<AGROW>
            times(end+1) = tr.tStim1;        codes(end+1) = cfg.code.stimOn1;%#ok<AGROW>
            if isfield(tr, 'code') && ~isempty(tr.code)    % tagging stimulus id
                times(end+1) = tr.tStim1 + 0.004; codes(end+1) = tr.code;    %#ok<AGROW>
            end
            if ~isnan(tr.tStim2)
                times(end+1) = tr.tStim2;    codes(end+1) = cfg.code.stimOn2;%#ok<AGROW>
            end
            if ~isnan(tr.tResp) && ~strcmp(tr.outcome, 'omission')
                times(end+1) = tr.tResp;     codes(end+1) = cfg.code.bhv;    %#ok<AGROW>
            end
            % Outcome marker, then the trial-end code. EVERY trial ends with
            % one of end1/end2/end3: trialdefGen pairs itiOn with those to
            % count trials, so a trial without an end code is not a trial.
            switch tr.outcome
                case 'correct'
                    times(end+1) = tr.tEnd - 0.2; codes(end+1) = cfg.code.rwd;  %#ok<AGROW>
                    endCode = cfg.code.end1;
                case 'incorrect'
                    times(end+1) = tr.tEnd - 0.2; codes(end+1) = cfg.code.pun;  %#ok<AGROW>
                    endCode = cfg.code.end2;
                case 'omission',  endCode = cfg.code.end2;   % no response given
                case 'aborted',   endCode = cfg.code.end3;   % broken off early
                otherwise,        endCode = cfg.code.end1;   % passive stream
            end
            times(end+1) = tr.tEnd; codes(end+1) = endCode; %#ok<AGROW>
        end
    end
    [times, order] = sort(times(:));
    codes = codes(order);
    codes = codes(:);

    EventRecord = struct();
    EventRecord.EventType          = codes;
    EventRecord.EventNumber        = (1:numel(codes))';
    EventRecord.TimeStamp          = round(times * 30000);      % raw samples
    EventRecord.TimeMsFromMidnight = times * 1000 + 9 * 3600 * 1000;  % 09:00 start
    EventRecord.TimeSecFromMidnight= times + 9 * 3600;
    EventRecord.TimeSource         = nan(size(codes));
    EventRecord.Details            = nan(size(codes));
    EventRecord.TimeBreak          = {[], []};
    eventList = struct('time', times, 'code', codes);
end

function [events, trialdef, condition] = localTrialTables(blocks, eventList, ~)
% events / trialdef / condition, in the shapes trialdefGen and NGL02 produce.
    trials = [blocks.trials];
    nTr = numel(trials);
    aligns = {'itiOn', 'stimOn1', 'stimOn2'};
    t0 = struct('itiOn', [trials.tITI], 'stimOn1', [trials.tStim1], ...
                'stimOn2', [trials.tStim2]);
    starts = [trials.tITI] - 0.5;
    ends   = [trials.tEnd];

    trialdef = cell(2, numel(aligns));
    events = struct();
    for i = 1:numel(aligns)
        zero = t0.(aligns{i});
        trialdef{1, i} = aligns{i};
        trialdef{2, i} = [starts(:) * 1000, ends(:) * 1000, zero(:) * 1000];
        ev = struct('code', {cell(nTr, 1)}, 'time', {cell(nTr, 1)});
        for k = 1:nTr
            in = eventList.time >= starts(k) & eventList.time <= ends(k);
            ev.code{k} = eventList.code(in);
            ev.time{k} = eventList.time(in) - zero(k);     % seconds, 0 = align
        end
        events.(aligns{i}) = ev;
    end
    % A trial with no stimOn2 (the tagging stream) keeps its row with a NaN
    % offset, which is what trialdefGen leaves when an alignment event is
    % missing from a trial. Row k is therefore trial k in every alignment,
    % and it is the trial parser that drops what it cannot cut.

    outcome = {trials.outcome};
    condition = struct();
    condition.correct   = double(strcmp(outcome, 'correct'))';
    condition.incorrect = double(strcmp(outcome, 'incorrect'))';
    condition.omission  = double(strcmp(outcome, 'omission'))';
    condition.aborted   = double(strcmp(outcome, 'aborted'))';
    condition.passive   = double(strcmp(outcome, 'passive'))';
    condition.responded = double(~isnan([trials.rt]))';
    condition.rt        = [trials.rt]';
    condition.block     = localBlockIndex(blocks)';
    condition.task      = localTaskIndex(blocks)';
end

function idx = localBlockIndex(blocks)
    idx = [];
    for b = 1:numel(blocks)
        idx = [idx, b * ones(1, numel(blocks(b).trials))]; %#ok<AGROW>
    end
end

function idx = localTaskIndex(blocks)
% 1 = dms, 2 = arena, 3 = nft, so a contrast can select a task.
    names = {'dms', 'arena', 'nft'};
    idx = [];
    for b = 1:numel(blocks)
        v = find(strcmp(names, blocks(b).task));
        idx = [idx, v * ones(1, numel(blocks(b).trials))]; %#ok<AGROW>
    end
end
