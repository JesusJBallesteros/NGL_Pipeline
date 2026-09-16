function report = checkFixture(root, varargin)
%CHECKFIXTURE  Run the pipeline's analyses on a fixture session and check the answers.
%
% PURPOSE:
%   The fixture is only worth having if the analyses recover what was planted.
%   This runs them - contrast, CSD, frequency tagging, jerk detection - and
%   compares each result against fixture_truth.mat, which was written by the
%   generator and not by any previous run. A check fails when the analysis
%   misses the planted effect OR when it finds one where nothing was planted;
%   both matter, and the second is the one a test usually forgets.
%
% USAGE:
%   checkFixture('D:\TESTSTUDY')
%   report = checkFixture(root, 'session', '19850214')
%   report = checkFixture(root, 'checks', {'trials','imu'}, 'strict', false)
%
% INPUTS:
%   root - the study folder given to makeFakeStudy.
%   Name/value pairs:
%     'subject'  default 'MRX'
%     'session'  default the first session found under the subject
%     'checks'   which to run; default all of
%                {'layout','isolation','trials','events','blocks','channels',
%                 'spikes','psth','imu','csd','contrast','nft'}
%     'strict'   true (default) errors if any check fails; false only reports
%
% OUTPUT:
%   report - struct array, one per check: .name .pass .detail
%
% NOTES:
%   * Runs the real functions, not copies of them. When the pipeline changes
%     and the planted answer stops coming back, this is what says so.
%   * The contrast and CSD checks take the slowest path (a wavelet TFR over
%     all trials), so the whole run is a couple of minutes on a short session.
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).

    p = inputParser;
    p.addParameter('subject', 'MRX');
    p.addParameter('session', '');
    p.addParameter('checks',  {'layout','isolation','trials','events','blocks', ...
                               'channels','spikes','psth','imu','csd', ...
                               'contrast','nft'});
    p.addParameter('strict',  true);
    p.parse(varargin{:});
    a = p.Results;
    if ischar(a.checks), a.checks = {a.checks}; end

    localAddPaths();
    S = localLoadSession(root, a.subject, a.session);
    fprintf('checkFixture: %s / %s (%d trials, %.0f s)\n', ...
            S.subject, S.session, S.truth.nTrials, S.truth.duration);

    report = struct('name', {}, 'pass', {}, 'detail', {});
    for k = 1:numel(a.checks)
        name = a.checks{k};
        try
            switch lower(name)
                case 'layout',   r = localCheckLayout(S);
                case 'isolation',r = localCheckIsolation(S);
                case 'trials',   r = localCheckTrials(S);
                case 'events',   r = localCheckEvents(S);
                case 'blocks',   r = localCheckBlocks(S);
                case 'channels', r = localCheckChannels(S);
                case 'spikes',   r = localCheckSpikes(S);
                case 'psth',     r = localCheckPSTH(S);
                case 'imu',      r = localCheckIMU(S);
                case 'csd',      r = localCheckCSD(S);
                case 'contrast', r = localCheckContrast(S);
                case 'nft',      r = localCheckNFT(S);
                otherwise, error('checkFixture:check', 'unknown check ''%s''.', name);
            end
        catch ME
            % A check that throws is a failure, not a crash of the run: the
            % remaining checks still carry information.
            r = struct('pass', false, 'detail', sprintf('threw: %s', ME.message));
        end
        r.name = name;
        report(end+1) = orderfields(r, {'name','pass','detail'}); %#ok<AGROW>
        fprintf('  [%s] %-9s %s\n', localTick(r.pass), name, r.detail);
    end

    nFail = sum(~[report.pass]);
    fprintf('checkFixture: %d/%d checks passed\n', numel(report) - nFail, numel(report));
    if nFail > 0 && a.strict
        error('checkFixture:failed', '%d fixture check(s) failed.', nFail);
    end
end

% ======================= setup =======================
function localAddPaths()
% The fixture lives inside the toolbox, so the toolbox is what it tests.
    here = fileparts(mfilename('fullpath'));
    reporoot = fileparts(fileparts(here));
    addpath(here, fullfile(reporoot, 'functions'));
    addpath(genpath(fullfile(reporoot, 'functions')));
    % calcFireRate and friends ship in a bundled toolbox, not under functions/.
    bd = fullfile(reporoot, 'toolboxes', 'BDPAT_NGL');
    if isfolder(bd), addpath(bd); end
    ftDir = fullfile(reporoot, 'toolboxes', 'fieldtrip_light');
    if isfolder(ftDir) && isempty(which('ft_freqanalysis'))
        addpath(ftDir); ft_defaults;
    end
end

function S = localLoadSession(root, subject, session)
    subjDir = fullfile(root, 'data', 'preprocessing', subject);
    assert(isfolder(subjDir), 'checkFixture:noSubject', ...
        'no fixture for subject %s under %s; run makeFakeStudy first.', subject, root);
    if isempty(session)
        d = dir(subjDir); d = d([d.isdir] & ~startsWith({d.name}, '.'));
        assert(~isempty(d), 'checkFixture:noSession', 'no sessions under %s.', subjDir);
        session = d(1).name;
    end
    S.root = root; S.subject = subject; S.session = session;
    S.preproc  = fullfile(root, 'data', 'preprocessing', subject, session);
    S.trial    = fullfile(root, 'data', 'trialSorted',   subject, session);
    S.spike    = fullfile(root, 'data', 'spikeSorted',   subject, session);
    S.analysis = fullfile(root, 'data', 'analysis',      subject, session);
    S.truth     = load(fullfile(S.preproc, 'fixture_truth.mat'));
    S.condition = localField(fullfile(S.trial, 'condition.mat'), 'condition');
    S.trialdef  = localField(fullfile(S.trial, 'trialdef.mat'),  'trialdef');
    S.events    = localField(fullfile(S.trial, 'events.mat'),    'events');
    S.blocks    = localField(fullfile(S.trial, 'blocks.mat'), 'sessionBlockTable');
    S.opt = localFixtureOpt(S);
end

function v = localField(file, name)
    L = load(file, '-mat');
    if isfield(L, name), v = L.(name); else, f = fieldnames(L); v = L.(f{1}); end
end

function [opt, input] = localFixtureOpt(S)
% A minimal options struct: schema defaults, then the few paths and TFR
% settings the analyses read. Going through the schema means the check uses
% the same defaults a user gets.
    opt = generateDefaultsFromSchema(optSchema());
    opt.SavFileName = S.session;
    opt.trialSorted = S.trial;
    opt.analysis    = S.analysis;
    opt.FolderProcDataMat = S.preproc;
    opt.KSchanMapFile = S.truth.chanMapFile;
    opt.TFRmethod  = 'wavelet';
    % One band from 8 Hz up: a 2 Hz wavelet needs more seconds than these
    % trials have, and the planted effects to test here are 20 Hz.
    opt.freqInterest = {8:1:40};
    % A second further out on each side than the effect needs. The trial-end
    % event sits at a fixed maximum, so every trial reaches +1.8 s whatever the
    % animal did - which is what makes a window this wide possible at all.
    opt.toi        = [-2 1.8];
    opt.timeResol  = 0.05;
    opt.alignto    = {'itiOn', 'stimOn1', 'stimOn2'};
    input = struct('analysisCode', '', 'areaMap', ...
                   struct('chanMapPath', S.truth.chanMapFile));
    opt.fixtureInput = input;
end

function s = localTick(pass)
    if pass, s = 'PASS'; else, s = 'FAIL'; end
end

function r = localResult(pass, fmt, varargin)
    r = struct('pass', logical(pass), 'detail', sprintf(fmt, varargin{:}));
end

% ======================= checks =======================
function r = localCheckLayout(S)
% The files a session must have, in the places the pipeline looks for them.
    want = {fullfile(S.preproc, [S.session '_FTcont.mat']), ...
            fullfile(S.preproc, 'EventRecord.mat'), ...
            fullfile(S.preproc, 'MotionData_raw.mat'), ...
            fullfile(S.trial, 'events.mat'), ...
            fullfile(S.trial, 'trialdef.mat'), ...
            fullfile(S.trial, 'condition.mat'), ...
            fullfile(S.trial, 'blocks.mat'), ...
            fullfile(S.trial, [S.session '_stimOn1.mat']), ...
            fullfile(S.spike, 'spike.mat')};
    missing = want(~cellfun(@isfile, want));
    if isempty(missing)
        r = localResult(true, '%d expected files present', numel(want));
    else
        [~, n, e] = cellfun(@fileparts, missing, 'UniformOutput', false);
        r = localResult(false, 'missing: %s', strjoin(strcat(n, e), ', '));
    end
end

function r = localCheckIsolation(S)
% The fixture must be impossible to mistake for data. Two halves: the marker
% is really on disk in every tree, and the guard that reads it behaves - 'all'
% skips it, naming it runs it, mixing it with a real subject is refused. The
% mixing rule is the one that matters, so it is tested rather than trusted.
    msg = {};
    trees = {S.preproc, S.trial, S.spike, S.analysis};
    for k = 1:numel(trees)
        subjFolder = fileparts(trees{k});
        if ~isSyntheticSubject(subjFolder)
            msg{end+1} = sprintf('no marker in %s', subjFolder); %#ok<AGROW>
        end
    end

    % A scratch study: this fixture subject beside an invented real one.
    tmp = fullfile(tempdir, ['NGLguard_' char(matlab.lang.internal.uuid)]);
    if isfolder(tmp), rmdir(tmp, 's'); end
    mkdir(fullfile(tmp, S.subject)); mkdir(fullfile(tmp, 'R01'));
    copyfile(fullfile(S.preproc, '..', 'SYNTHETIC_DATA.txt'), ...
             fullfile(tmp, S.subject, 'SYNTHETIC_DATA.txt'));
    mk = @(names) struct('name', names, 'isdir', {true});
    available = [mk({S.subject}), mk({'R01'})];

    % (1) 'all' drops it.
    got = guardSyntheticSubjects(struct('subjects', available), tmp, true, available);
    if ~isequal({got.subjects.name}, {'R01'})
        msg{end+1} = '''all'' did not drop the synthetic subject';
    end
    if ~isfield(got, 'fixtureRun') || got.fixtureRun
        msg{end+1} = '''all'' was marked as a fixture run';
    end

    % (2) named on its own, in its own study, it runs - and says so.
    solo = fullfile(tmp, 'solo');
    mkdir(fullfile(solo, S.subject));
    copyfile(fullfile(tmp, S.subject, 'SYNTHETIC_DATA.txt'), ...
             fullfile(solo, S.subject, 'SYNTHETIC_DATA.txt'));
    onlyFixture = mk({S.subject});
    got = evalc(['g = guardSyntheticSubjects(struct(''subjects'', onlyFixture), ' ...
                 'solo, false, onlyFixture);']); %#ok<NASGU>
    if ~g.fixtureRun, msg{end+1} = 'a fixture-only run was not flagged'; end

    % (3) mixed with a real subject: refused.
    mixed = [mk({S.subject}), mk({'R01'})];
    try
        guardSyntheticSubjects(struct('subjects', mixed), tmp, false, mixed);
        msg{end+1} = 'synthetic + real was ACCEPTED';
    catch ME
        if ~strcmp(ME.identifier, 'NGL:mixedSyntheticSubjects')
            msg{end+1} = sprintf('mixing failed with the wrong error (%s)', ME.identifier);
        end
    end

    % (4) fixture alone, but inside a study that holds a real subject: refused,
    %     because its study-level outputs would land on top of theirs.
    try
        guardSyntheticSubjects(struct('subjects', onlyFixture), tmp, false, available);
        msg{end+1} = 'a fixture inside a real study was ACCEPTED';
    catch ME
        if ~strcmp(ME.identifier, 'NGL:fixtureInRealStudy')
            msg{end+1} = sprintf('wrong error for fixture-in-real-study (%s)', ME.identifier);
        end
    end
    rmdir(tmp, 's');

    if isempty(msg)
        r = localResult(true, ['marked in %d trees; ''all'' skips it, alone it ' ...
            'runs flagged, mixed with real data it is refused'], numel(trees));
    else
        r = localResult(false, '%s', strjoin(msg, '; '));
    end
end

function r = localCheckTrials(S)
% Trial bookkeeping: counts agree, outcomes are exclusive, and latencies stay
% inside the response windows the design specifies.
    n = S.truth.nTrials;
    msg = {};
    if numel(S.condition.correct) ~= n
        msg{end+1} = sprintf('condition has %d rows, truth says %d', ...
                             numel(S.condition.correct), n);
    end
    if size(S.trialdef{2,1}, 1) ~= n
        msg{end+1} = sprintf('trialdef has %d rows', size(S.trialdef{2,1}, 1));
    end
    % Outcomes partition the trials: every trial has exactly one.
    oc = [S.condition.correct, S.condition.incorrect, S.condition.omission, ...
          S.condition.aborted, S.condition.passive];
    if any(sum(oc, 2) ~= 1)
        msg{end+1} = sprintf('%d trial(s) with no or several outcomes', ...
                             sum(sum(oc, 2) ~= 1));
    end
    % Latency: present exactly when a response was made, and within the window.
    resp = S.condition.responded == 1;
    if any(isnan(S.condition.rt(resp))) || any(~isnan(S.condition.rt(~resp)))
        msg{end+1} = 'rt and responded disagree';
    end
    trials = [S.truth.blocks.trials];
    for task = {'dms', 'arena'}
        isT = strcmp({trials.task}, task{1});
        cap = 2.0 * strcmp(task{1}, 'dms') + 3.0 * strcmp(task{1}, 'arena');
        rt = [trials(isT).rt]; rt = rt(~isnan(rt));
        % An aborted trial's peck lands BEFORE the test appears, so its latency
        % is negative by construction; everything else falls inside the window.
        if any(rt > cap + 1e-9)
            msg{end+1} = sprintf('%s rt above its %g s window', task{1}, cap); %#ok<AGROW>
        end
    end
    aborted = strcmp({trials.outcome}, 'aborted');
    if any(aborted) && any([trials(aborted).rt] >= 0)
        msg{end+1} = 'an aborted trial responded after the test appeared';
    end
    % The end event sits at a fixed maximum: two trials of the same task with
    % the same delay must end the same distance from their last stimulus, no
    % matter when the animal responded.
    isD = strcmp({trials.task}, 'dms');
    tail = [trials(isD).tEnd] - [trials(isD).tStim2];
    if std(tail) > 1e-9
        msg{end+1} = sprintf('dms trial end is not fixed (spread %.3f s)', std(tail));
    end
    % Trial lengths must still vary across trials: a fixture whose trials are
    % all identical would pass code that assumes they are.
    len = [trials.tEnd] - [trials.tITI];
    if std(len(isD)) < 0.05
        msg{end+1} = 'dms trial lengths are too uniform';
    end
    if isempty(msg)
        r = localResult(true, ['%d trials, %d correct / %d incorrect / %d omitted / ' ...
            '%d aborted / %d passive; dms length %.2f-%.2f s'], n, ...
            sum(S.condition.correct), sum(S.condition.incorrect), ...
            sum(S.condition.omission), sum(S.condition.aborted), ...
            sum(S.condition.passive), min(len), max(len));
        r.detail = sprintf('%s; dms ends %.2f s after stimOn2, always', ...
                           r.detail, mean(tail));
    else
        r = localResult(false, '%s', strjoin(msg, '; '));
    end
end

function r = localCheckEvents(S)
% The digital-input record: every trial opens with itiOn and closes with an
% end code, which is the pairing trialdefGen relies on to find trials at all.
    ER = localField(fullfile(S.preproc, 'EventRecord.mat'), 'EventRecord');
    codes = ER.EventType;
    t = ER.TimeSecFromMidnight;
    msg = {};
    if any(diff(t) < 0), msg{end+1} = 'event times are not monotonic'; end
    nStart = sum(codes == 0);
    nEnd   = sum(codes == 4 | codes == 10 | codes == 15);
    if nStart ~= S.truth.nTrials || nEnd ~= S.truth.nTrials
        msg{end+1} = sprintf('%d itiOn and %d end codes for %d trials', ...
                             nStart, nEnd, S.truth.nTrials);
    end
    % Every reward or punishment must follow a peck, closely. An outcome that
    % precedes its response, or floats seconds away from one, is a broken
    % record - and it is the kind of thing that still yields plausible plots.
    tPeck = t(codes == 3);
    tOut  = t(codes == 7 | codes == 11);
    if isempty(tPeck) || isempty(tOut)
        msg{end+1} = 'no pecks or no outcome markers';
    else
        lag = arrayfun(@(x) min([x - tPeck(tPeck <= x); Inf]), tOut);
        if any(lag > 0.5) || any(lag <= 0)
            msg{end+1} = sprintf('%d outcome marker(s) not just after a peck', ...
                                 sum(lag > 0.5 | lag <= 0));
        end
    end
    % Tagging stimulus codes are project codes (>= 16) and must appear.
    nTag = sum(ismember(codes, 8001:8006));
    if nTag == 0 && any(strcmp({S.truth.blocks.task}, 'nft'))
        msg{end+1} = 'no tagging stimulus codes in the record';
    end
    known = [16 17 18 19 20];                       % block markers
    if any(codes > 15 & codes < 8001 & ~ismember(codes, known))
        msg{end+1} = 'unexpected project codes';
    end
    if isempty(msg)
        r = localResult(true, ['%d events, %d trials paired, %d tagging codes; ' ...
            'outcomes follow the peck by %.0f ms'], numel(codes), nStart, nTag, ...
            1000 * median(lag));
    else
        r = localResult(false, '%s', strjoin(msg, '; '));
    end
end

function r = localCheckBlocks(S)
% The block markers are read once per session and must land on the blocks that
% were actually run: right label, right times, right trials. This is the table
% every later analysis consults to know which trials it may touch.
    tbl = S.blocks;
    want = S.truth.blocks;
    msg = {};
    if numel(tbl) ~= numel(want)
        msg{end+1} = sprintf('%d blocks found, %d run', numel(tbl), numel(want));
    else
        for k = 1:numel(want)
            if ~strcmp(tbl(k).label, want(k).task)
                msg{end+1} = sprintf('block %d labelled %s, ran %s', k, ...
                                     tbl(k).label, want(k).task); %#ok<AGROW>
            end
            % The markers sit half a second outside the block, so the recovered
            % range must contain the block and not drift far from it.
            if tbl(k).tStart > want(k).tStart || tbl(k).tEnd < want(k).tEnd ...
                    || want(k).tStart - tbl(k).tStart > 1 ...
                    || tbl(k).tEnd - want(k).tEnd > 1
                msg{end+1} = sprintf('block %d spans %.1f-%.1f s, ran %.1f-%.1f s', ...
                    k, tbl(k).tStart, tbl(k).tEnd, want(k).tStart, want(k).tEnd); %#ok<AGROW>
            end
            if tbl(k).nTrials ~= numel(want(k).trials)
                msg{end+1} = sprintf('block %d holds %d trials, ran %d', k, ...
                    tbl(k).nTrials, numel(want(k).trials)); %#ok<AGROW>
            end
        end
    end
    % Blocks must tile the session without overlapping: a trial belongs to one
    % block, and an analysis that selects by block must not see it twice.
    if numel(tbl) > 1 && any([tbl(2:end).tStart] < [tbl(1:end-1).tEnd])
        msg{end+1} = 'blocks overlap in time';
    end
    idx = [tbl.trialIdx];
    if numel(unique(idx)) ~= numel(idx)
        msg{end+1} = 'a trial belongs to more than one block';
    end
    if numel(idx) ~= S.truth.nTrials
        msg{end+1} = sprintf('%d of %d trials fall in a block', ...
                             numel(idx), S.truth.nTrials);
    end
    % The table's job: an analysis asks it which trials it may use, and gets
    % only the trials of tasks that carry that analysis - none for a task that
    % does not, and none at all for an analysis nobody declared.
    gate = {'contrast', 'dms'; 'nft', 'nft'; 'csd', 'dms arena'};
    for g = 1:size(gate, 1)
        got = localBlocksFor(S, gate{g, 1});
        wantBlocks = strsplit(gate{g, 2});
        expect = [];
        for wb = wantBlocks
            j = find(strcmp({tbl.label}, wb{1}), 1);
            if ~isempty(j), expect = [expect, tbl(j).trialIdx]; end %#ok<AGROW>
        end
        if ~isequal(sort(got(:)), sort(expect(:)))
            msg{end+1} = sprintf('%s gated to %d trials, expected %d', ...
                                 gate{g, 1}, numel(got), numel(expect)); %#ok<AGROW>
        end
    end
    if ~isempty(localBlocksFor(S, 'noSuchAnalysis'))
        msg{end+1} = 'an undeclared analysis was granted trials';
    end
    if isempty(msg)
        parts = arrayfun(@(b) sprintf('%s %.0f-%.0f s (%d)', b.label, ...
                         b.tStart, b.tEnd, b.nTrials), tbl, 'UniformOutput', false);
        r = localResult(true, '%s; gating holds', strjoin(parts, ', '));
    else
        r = localResult(false, '%s', strjoin(msg, '; '));
    end
end

function r = localCheckChannels(S)
% The bad channels have to be bad in the data, not just in the notes. A dead
% channel that still carries the planted response is the worst case of all: its
% noise floor is the lowest on the probe, so every signal-to-noise measure
% ranks it first and a figure picks it to show.
    FT = localLoadFT(S, '');                        % continuous
    x = FT.trial{1};
    v = var(double(x), 0, 2);
    iDead  = find(strcmp(FT.label, S.truth.lfp.deadChannel), 1);
    iNoisy = find(strcmp(FT.label, S.truth.lfp.noisyChannel), 1);
    others = setdiff(1:numel(v), [iDead, iNoisy]);
    typical = median(v(others));
    msg = {};
    if v(iDead) > 0.05 * typical
        msg{end+1} = sprintf('%s is not dead (%.0f%% of a typical channel)', ...
                             S.truth.lfp.deadChannel, 100 * v(iDead) / typical);
    end
    if v(iNoisy) < 5 * typical
        msg{end+1} = sprintf('%s is not noisy (%.1fx a typical channel)', ...
                             S.truth.lfp.noisyChannel, v(iNoisy) / typical);
    end
    % And it must carry no planted response either: the tagging block drives
    % every other NCL channel, so this is where a leak would show.
    blocks = S.truth.blocks;
    b = find(strcmp({blocks.task}, 'nft'), 1);
    if ~isempty(b) && ~isempty(iDead)
        trials = blocks(b).trials;
        in = [trials.rate] == min([trials.rate]);
        win = [min([trials(in).tITI]), max([trials(in).tEnd])];
        epochs = nftEpochs(FT, win, min([trials.rate]));
        spc = computeTaggingSpectrum(epochs, S.opt);
        resp = taggingResponse(spc, min([trials.rate]), S.opt);
        zDead = resp.z(iDead, 1);
        if zDead > 3
            msg{end+1} = sprintf('the dead channel answers the tagging (z=%.1f)', zDead);
        end
    else
        zDead = NaN;
    end
    if isempty(msg)
        r = localResult(true, ['%s dead (%.1f%% of typical variance, tagging ' ...
            'z=%.1f), %s noisy (%.0fx)'], S.truth.lfp.deadChannel, ...
            100 * v(iDead) / typical, zDead, S.truth.lfp.noisyChannel, ...
            v(iNoisy) / typical);
    else
        r = localResult(false, '%s', strjoin(msg, '; '));
    end
end

function r = localCheckSpikes(S)
% The units are planted with specific jobs: some fire after stimOn1, one after
% reward, one drifts in amplitude, one is labelled noise. Each is checked
% against a unit that has no such job, so a general rate increase cannot pass.
    spike = localField(fullfile(S.spike, 'spike.mat'), 'spike');
    truth = S.truth.spikes;
    ids = cellfun(@str2double, spike.label);
    msg = {};
    if numel(ids) ~= truth.nUnits || ~isequal(sort(ids(:))', sort(truth.ids(:))')
        msg{end+1} = 'unit ids do not match the plan';
    end
    % Curation labels must survive: a fixture without a noise cluster would let
    % code that ignores curation pass.
    if ~any(strcmp(spike.HumanLabel, 'noise')), msg{end+1} = 'no noise cluster'; end
    if ~any(strcmp(spike.bc_unitType, 'GOOD')),  msg{end+1} = 'no GOOD unit';    end

    trials = [S.truth.blocks.trials];
    stim1 = [trials.tStim1]; stim1 = stim1(~isnan(stim1));
    % Driven units: more spikes 40-70 ms after stimOn1 than in the 100 ms
    % before it. The undriven ones must not show the same.
    drivenRatio = @(u) localPeriRatio(spike.timestamp{u}, stim1, [0.04 0.07], [-0.1 0]);
    isDriven = ismember(ids, truth.drivenUnits);
    rDriven = arrayfun(drivenRatio, find(isDriven));
    rOther  = arrayfun(drivenRatio, find(~isDriven));
    if ~(min(rDriven) > 1.5 && max(rOther) < 1.5)
        msg{end+1} = sprintf('stimOn1 drive: driven %.2f-%.2f, others up to %.2f', ...
                             min(rDriven), max(rDriven), max(rOther));
    end
    % The drifting unit loses amplitude across the session; nobody else does.
    uD = find(ids == truth.driftingUnit, 1);
    amp = spike.ampl{uD}; ts = spike.timestamp{uD};
    half = ts > median(ts);
    lost = 1 - mean(amp(half)) / mean(amp(~half));
    if lost < 0.15
        msg{end+1} = sprintf('drifting unit lost only %.0f%% of its amplitude', 100*lost);
    end
    if isempty(msg)
        r = localResult(true, ['%d units; stimOn1 drive %.1fx (others %.1fx), ' ...
            'drift -%.0f%% over the session'], numel(ids), mean(rDriven), ...
            max(rOther), 100*lost);
    else
        r = localResult(false, '%s', strjoin(msg, '; '));
    end
end

function ratio = localPeriRatio(ts, events, win, base)
% Spikes per second in a window after each event, over the same in a baseline
% window before it.
    ts = ts(:);
    n = @(w) sum(arrayfun(@(e) sum(ts >= e + w(1) & ts < e + w(2)), events(:)));
    rWin  = n(win)  / (numel(events) * diff(win));
    rBase = n(base) / (numel(events) * diff(base));
    ratio = rWin / max(rBase, eps);
end

function r = localCheckPSTH(S)
% The spike side end to end, through the pipeline's own path: sort2trials cuts
% the trains into trials, calcFireRate turns them into rates. The driven unit
% must rise after stimOn1 in the two tasks that present a stimulus to respond
% to, and the undriven one must not - in every task.
    spike = localField(fullfile(S.spike, 'spike.mat'), 'spike');
    neurons = sort2trials(spike, S.trialdef, S.opt);

    ids = cellfun(@str2double, spike.label);
    uDriven = find(ids == S.truth.spikes.drivenUnits(1), 1);
    % The control pools EVERY unit with no stimulus job. One of them alone is
    % a handful of spikes in the arena's 24 trials, and a ratio built on that
    % wanders by chance - which would make the control fail on its own noise
    % rather than on anything the pipeline did.
    uQuiet  = find(~ismember(ids, S.truth.spikes.drivenUnits) & ...
                   ids ~= S.truth.spikes.rewardUnit & ...
                   ids ~= S.truth.spikes.taggedUnit);

    param = struct('binSize', 50, 'stepSz', 10, 'interval', [-500 500], ...
                   'smpRate', 1000, 'baseline', 500, 'plot', false);
    lines = {}; ok = true;
    for task = {'dms', 'arena', 'nft'}
        trialIdx = localBlockTrials(S, task{1});
        if isempty(trialIdx), continue; end
        [tAxis, rDriven] = localPSTH(neurons.stimOn1{uDriven}, trialIdx, param);
        rQ = zeros(numel(uQuiet), numel(rDriven));
        for q = 1:numel(uQuiet)
            [~, rQ(q, :)] = localPSTH(neurons.stimOn1{uQuiet(q)}, trialIdx, param);
        end
        rQuiet = sum(rQ, 1);
        % The planted drive is 40-70 ms after the stimulus; compare it with
        % this unit's own pre-stimulus rate rather than with another unit's.
        gain  = @(v) mean(v(tAxis >= 20 & tAxis <= 120)) / ...
                     max(mean(v(tAxis >= -400 & tAxis <= -50)), eps);
        gD = gain(rDriven); gQ = gain(rQuiet);
        % Every task presents a stimulus, and this unit answers a stimulus -
        % including in the passive stream, where there is nothing to do about
        % it. The control is the other unit, which answers in none of them.
        ok = ok && gD > 1.5 && gQ < 1.5;
        lines{end+1} = sprintf('%s %.1fx (quiet %.1fx)', task{1}, gD, gQ); %#ok<AGROW>
    end
    r = localResult(ok, 'unit %d after stimOn1: %s (control = %d pooled units)', ...
                    S.truth.spikes.drivenUnits(1), strjoin(lines, ', '), numel(uQuiet));
end

function [t, rate] = localPSTH(trialSpikes, trialIdx, param)
% calcFireRate over one set of trials, averaged. Spike times arrive in ms
% relative to the alignment, which is what it expects.
    sel = trialSpikes(trialIdx);
    sel(cellfun(@isempty, sel)) = {NaN};
    fr = calcFireRate(sel, param, struct());
    if iscell(fr), fr = cell2mat(fr); end
    rate = mean(fr, 1, 'omitnan');
    t = param.interval(1) + param.binSize/2 + ...
        (0:size(fr, 2)-1) * param.stepSz;
end

function idx = localBlockTrials(S, label)
% The trials of one block, from the session's block table.
    k = find(strcmp({S.blocks.label}, label), 1);
    if isempty(k), idx = []; else, idx = S.blocks(k).trialIdx; end
end

function r = localCheckIMU(S)
% The preprocessing estimate_pecking does - high-pass the acceleration,
% differentiate to jerk, take the magnitude - and then a plain threshold on
% that magnitude. What is being tested is the fixture: that each planted peck
% leaves a jerk transient far above the noise and nothing else does. The
% jerk magnitude is what the pecking detectors consume, so if this check
% passes, they have something to find; if it fails, the IMU generator broke.
    L = load(fullfile(S.preproc, 'MotionData_raw.mat'), '-mat');
    raw = L.raw;
    fs = raw.fs;
    acc = [raw.acc.X(:), raw.acc.Y(:), raw.acc.Z(:)];
    [b, a] = butter(3, 20 / (fs/2), 'high');        % opt.hpass = 20 Hz
    aHP = filtfilt(b, a, acc);
    jerk = zeros(size(aHP));
    jerk(2:end-1, :) = (aHP(3:end, :) - aHP(1:end-2, :)) / (2 / fs);
    mag = sqrt(sum(jerk.^2, 2));

    % The pecks sit ~50 MAD above the noise, so where the threshold goes
    % between those two is not a tuning question; 15 is well inside the gap.
    thr = median(mag) + 15 * mad(mag, 1);
    refr = round(0.1 * fs);
    found = localThresholdEvents(mag, thr, refr) / fs;
    want = S.truth.imu.peckTimes(:);
    tol = 0.05;                                     % the jerk is 30 ms long
    hit = arrayfun(@(x) any(abs(found - x) <= tol), want);
    matched = arrayfun(@(x) any(abs(want - x) <= tol), found);
    recall = mean(hit); prec = mean(matched);
    sep = median(arrayfun(@(x) max(mag(max(1, round(x*fs)-20) : ...
                    min(numel(mag), round(x*fs)+60))), want)) / mad(mag, 1);
    ok = recall >= 0.95 && prec >= 0.90;
    r = localResult(ok, ['recall %.0f%% of %d pecks, precision %.0f%%; ' ...
        'peck jerk is %.0fx the noise MAD'], 100*recall, numel(want), 100*prec, sep);
end

function idx = localThresholdEvents(x, thr, refr)
% One event per above-threshold excursion, then a refractory sweep keeping the
% largest. Deliberately the simplest detector that can work: a fixture check
% should not depend on a detector's tuning.
    above = x(:) > thr;
    d = diff([false; above; false]);
    starts = find(d == 1); stops = find(d == -1) - 1;
    idx = zeros(numel(starts), 1);
    for k = 1:numel(starts)
        [~, rel] = max(x(starts(k):stops(k)));
        idx(k) = starts(k) + rel - 1;
    end
    keep = true(size(idx));
    last = -Inf;
    for k = 1:numel(idx)
        if idx(k) - last < refr, keep(k) = false; else, last = idx(k); end
    end
    idx = idx(keep);
end

function r = localCheckCSD(S)
% The planted sink is at a known depth on shank 1. computeCSD must put its
% strongest sink there, shortly after stimOn1.
    FT = localLoadFT(S, 'stimOn1');
    % Trials differ in length by design, and a CSD needs one time axis, so the
    % evoked window is cut out first - what NGL07 would have to do as well.
    [FT, keep] = localCommonWindow(FT, [-0.2 0.4]);
    shanks = lfpChannelGeometry(FT, S.opt.fixtureInput, S.opt, ...
                                'chanMap', S.truth.chanMapFile);
    sh = shanks([shanks.shank] == S.truth.lfp.csdShank);
    assert(~isempty(sh), 'no shank %d in the geometry', S.truth.lfp.csdShank);
    % Only the trials that have a stimulus-evoked response: the tagging block
    % has no stimOn1-locked laminar profile planted.
    % Which trials may carry a CSD is a property of the task, so it is asked of
    % the block table rather than assumed of the session.
    idx = FT.trialinfo(keep);
    allowed = localBlocksFor(S, 'csd');
    useTrial = ismember(idx(:), allowed(:));
    assert(any(useTrial), 'no trials in a block that allows a CSD');
    csd = computeCSD(FT, sh, S.opt, 'trials', useTrial);

    win = csd.time >= 0.02 & csd.time <= 0.15;     % the evoked window
    [~, iT] = min(min(csd.csd(:, win), [], 1));    % strongest sink in time
    tWin = csd.time(win); tSink = tWin(iT);
    prof = csd.csd(:, win);
    [~, iD] = min(prof(:, iT));                    % ... and its depth
    got = csd.depth(iD);
    want = S.truth.lfp.csdSinkDepth;
    spacing = median(diff(csd.depth));
    ok = abs(got - want) <= 1.5 * spacing;         % within a contact and a half
    r = localResult(ok, 'sink at %g um (planted %g, spacing %g) at %+.0f ms', ...
                    got, want, spacing, 1000*tSink);
end

function r = localCheckContrast(S)
% The 20 Hz burst was planted on CORRECT dms trials in NCL only. The contrast
% must find it in NCL and must NOT find it in STR - a test that only checks
% the positive half would pass on a function that flags everything.
    FT = localLoadFT(S, 'stimOn2');
    cond = localSubsetCondition(S.condition, FT.trialinfo);
    % dms trials only: the arena block shares this alignment but carries theta,
    % not the beta burst, and mixing tasks would dilute the contrast.
    allowed = localBlocksFor(S, 'contrast');
    isDMS = ismember(FT.trialinfo(:), allowed(:));
    TFR = computeTrialparsedTFR(FT, cond, struct(), S.opt, 'stimOn2');
    band = TFR{1};

    out = struct('area', {}, 'p', {}, 'f', {}, 'effect', {});
    for area = {'NCL', 'STR'}
        labels = FT.label(strcmp(FT.chanArea, area{1}));
        sub = ft_selectdata(struct('channel', {labels}), band);
        spec = parseTrialContrast(struct('A', cond.correct(:) > 0 & isDMS, ...
                                         'B', cond.incorrect(:) > 0 & isDMS, ...
                                         'labelA', 'correct', 'labelB', 'incorrect'), ...
                                  [], size(sub.powspctrm, 1));
        res = computeTFRcontrast(sub, spec, S.opt, 'area', area{1}, 'align', 'stimOn2');
        out(end+1) = localClusterPeak(res, area{1}); %#ok<AGROW>
    end
    ncl = out(strcmp({out.area}, 'NCL'));
    str = out(strcmp({out.area}, 'STR'));
    beta = S.truth.lfp.betaBand;
    % The negative control is about SIZE, not only about significance. A
    % cluster test on 16 channels finds something below 0.05 in the undriven
    % area every so often - that is what 0.05 means - and a fixture that fails
    % when it does would be testing the weather. What must hold is that
    % anything STR shows is a fraction of what NCL shows.
    inBand = ncl.f >= beta(1) && ncl.f <= beta(2);
    if str.effect > 0
        ratio = ncl.effect / str.effect;
        tail = sprintf('%s, %.0fx smaller', localNone(str.f, str.p), ratio);
    else
        ratio = Inf;                                % nothing there to compare
        tail = localNone(str.f, str.p);
    end
    ok = ncl.p <= 0.01 && inBand && ratio >= 5;
    r = localResult(ok, 'NCL cluster p=%.3f at %.0f Hz (planted 20, band %g-%g); STR %s', ...
                    ncl.p, ncl.f, beta(1), beta(2), tail);
end

function r = localCheckNFT(S)
% The tagging block drives NCL at the stimulation rate and leaves STR alone.
    FTc = localLoadFT(S, '');                       % continuous
    [~, blk] = localBlocksFor(S, 'nft');
    assert(~isempty(blk), 'no block in this session allows a tagging analysis');
    blocks = S.truth.blocks;
    b = find(strcmp({blocks.task}, blk(1).label), 1);
    trials = blocks(b).trials;
    rates = unique([trials.rate]);
    lines = {}; ok = true;
    for rr = rates(:)'
        in = [trials.rate] == rr;
        win = [min([trials(in).tITI]), max([trials(in).tEnd])];
        epochs = nftEpochs(FTc, win, rr);
        spec = computeTaggingSpectrum(epochs, S.opt);
        resp = taggingResponse(spec, rr, S.opt);
        isNCL = strcmp(FTc.chanArea, 'NCL');
        % Compare the driven area against the undriven one at the base
        % frequency; the dead and noisy channels are deliberately included,
        % so the median is what is read, not the maximum.
        zNCL = median(resp.z(isNCL, 1));
        zSTR = median(resp.z(~isNCL, 1));
        % The response spreads over harmonics, and how much of it sits in the
        % fundamental depends on the pulse's duty cycle - so the test is that
        % the driven area has significant harmonics and the undriven one has
        % none, with the fundamental's z reported alongside.
        nNCL = median(resp.nSignificant(isNCL));
        nSTR = median(resp.nSignificant(~isNCL));
        % The driven channels reach nearly every harmonic; the undriven ones
        % clip one now and then, which is what 16 channels x 8 harmonics of
        % testing costs - and one of those channels is deliberately noisy.
        % Demanding zero there would be testing the noise, not the response.
        ok = ok && zNCL > 3 && nNCL >= 4 && zSTR < 3 && nSTR <= 1;
        lines{end+1} = sprintf('%.1f Hz: NCL z %.1f, %.1f harmonics; STR z %.1f, %.1f', ...
                               rr, zNCL, nNCL, zSTR, nSTR); %#ok<AGROW>
    end
    r = localResult(ok, '%s', strjoin(lines, '; '));
end

% ======================= check helpers =======================
function FT = localLoadFT(S, alignName)
    if isempty(alignName)
        f = fullfile(S.preproc, [S.session '_FTcont.mat']);
    else
        f = fullfile(S.trial, [S.session '_' alignName '.mat']);
    end
    FT = localField(f, 'FT_data');
    if isfield(FT, 'FT_data'), FT = FT.FT_data; end
end

function [FT, keep] = localCommonWindow(FT, win)
% Cut every trial to the same window and drop the ones that do not span it.
% Variable-length trials are real (a response can come at any time); anything
% that averages across trials in time needs them squared off first.
    FT = ft_redefinetrial(struct('toilim', win), FT);
    % Keep only trials that cover the FULL window. The most common length is
    % not the criterion here: the tagging block has the most trials and the
    % shortest ones, so the majority would pick exactly the wrong set.
    want = floor(diff(win) * FT.fsample);
    n = cellfun(@numel, FT.time);
    keep = n >= want;
    cut = min(n(keep));                             % one length for all of them
    for k = find(keep(:)')
        FT.trial{k} = FT.trial{k}(:, 1:cut);
        FT.time{k}  = FT.time{k}(1:cut);
    end
    FT.trial = FT.trial(keep);
    FT.time  = FT.time(keep);
    if isfield(FT, 'sampleinfo') && size(FT.sampleinfo, 1) == numel(keep)
        FT.sampleinfo = FT.sampleinfo(keep, :);
    end
end

function [trialIdx, blk] = localBlocksFor(S, analysis)
% The trials this analysis is allowed to use, taken from the block table and
% the study's declaration of which task carries which analysis. An analysis
% that asks for trials it has no business in gets none, which is the point:
% a tagging spectrum on a delay-match block is not a weaker result, it is a
% meaningless one.
    trialIdx = []; blk = S.blocks([]);
    for k = 1:numel(S.blocks)
        label = S.blocks(k).label;
        if ~isfield(S.truth.taskAnalyses, label), continue; end
        if any(strcmp(S.truth.taskAnalyses.(label), analysis))
            trialIdx = [trialIdx, S.blocks(k).trialIdx]; %#ok<AGROW>
            blk(end+1) = S.blocks(k);                    %#ok<AGROW>
        end
    end
end

function cond = localSubsetCondition(cond, idx)
% Trials the parser dropped are gone from the FT data but still present in the
% condition struct; .trialinfo says which ones survived.
    f = fieldnames(cond);
    for k = 1:numel(f)
        v = cond.(f{k});
        if isnumeric(v) && numel(v) >= max(idx), cond.(f{k}) = v(idx); end
    end
end

function out = localClusterPeak(res, area)
% The most significant cluster: its frequency, its p, and how big the
% difference inside it actually is. NaN / zero when the test found nothing,
% which for STR is the expected answer.
    out = struct('area', area, 'p', NaN, 'f', NaN, 'effect', 0);
    if isempty(res.stat) || ~isfield(res.stats, 'anySignificant') ...
            || ~res.stats.anySignificant || ~any(res.mask(:))
        return
    end
    pos = res.diff; pos(~res.mask) = NaN;
    [out.effect, i] = max(abs(pos(:)));
    [iF, ~] = ind2sub(size(pos), i);
    out.f = res.freq(iF);
    out.p = min([res.stats.pPos, res.stats.pNeg, Inf]);
    if ~isfinite(out.p), out.p = NaN; end
end

function s = localNone(f, p)
    if isnan(f), s = 'no cluster (correct)';
    else, s = sprintf('cluster p=%.3f at %.0f Hz (should be none)', p, f);
    end
end
