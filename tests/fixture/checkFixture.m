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
%                {'layout','trials','events','spikes','imu','csd','contrast','nft'}
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
    p.addParameter('checks',  {'layout','trials','events','spikes','imu','csd', ...
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
                case 'trials',   r = localCheckTrials(S);
                case 'events',   r = localCheckEvents(S);
                case 'spikes',   r = localCheckSpikes(S);
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
    opt.toi        = [-1 0.8];
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
        if any(rt > cap + 1e-9) || any(rt <= 0)
            msg{end+1} = sprintf('%s rt outside (0, %g]', task{1}, cap); %#ok<AGROW>
        end
    end
    % Trial lengths must vary: a fixture with fixed-length trials would pass
    % code that assumes fixed-length trials, which is what it must not do.
    len = [trials.tEnd] - [trials.tITI];
    if std(len(strcmp({trials.task}, 'dms'))) < 0.05
        msg{end+1} = 'dms trial lengths are too uniform';
    end
    if isempty(msg)
        r = localResult(true, ['%d trials, %d correct / %d incorrect / %d omitted / ' ...
            '%d aborted / %d passive; dms length %.2f-%.2f s'], n, ...
            sum(S.condition.correct), sum(S.condition.incorrect), ...
            sum(S.condition.omission), sum(S.condition.aborted), ...
            sum(S.condition.passive), min(len), max(len));
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
    % Tagging stimulus codes are project codes (>= 16) and must appear.
    nTag = sum(ismember(codes, 8001:8006));
    if nTag == 0 && any(strcmp({S.truth.blocks.task}, 'nft'))
        msg{end+1} = 'no tagging stimulus codes in the record';
    end
    if any(codes > 15 & codes < 8001), msg{end+1} = 'unexpected project codes'; end
    if isempty(msg)
        r = localResult(true, '%d events, %d trials paired, %d tagging codes', ...
                        numel(codes), nStart, nTag);
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
    idx = FT.trialinfo(keep);
    trials = [S.truth.blocks.trials];
    useTrial = ~strcmp({trials(idx).task}, 'nft');
    csd = computeCSD(FT, sh, S.opt, 'trials', useTrial(:));

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
    trials = [S.truth.blocks.trials];
    isDMS = strcmp({trials(FT.trialinfo).task}, 'dms')';
    TFR = computeTrialparsedTFR(FT, cond, struct(), S.opt, 'stimOn2');
    band = TFR{1};

    out = struct('area', {}, 'p', {}, 'f', {});
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
    ok = ncl.p < 0.05 && ncl.f >= beta(1) && ncl.f <= beta(2) && isnan(str.f);
    r = localResult(ok, 'NCL cluster p=%.3f at %.0f Hz (planted 20, band %g-%g); STR %s', ...
                    ncl.p, ncl.f, beta(1), beta(2), ...
                    localNone(str.f, str.p));
end

function r = localCheckNFT(S)
% The tagging block drives NCL at the stimulation rate and leaves STR alone.
    FTc = localLoadFT(S, '');                       % continuous
    blocks = S.truth.blocks;
    b = find(strcmp({blocks.task}, 'nft'), 1);
    assert(~isempty(b), 'this session has no tagging block');
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
% The frequency of the most significant cluster, and its p; NaN when the test
% found nothing, which for STR is the expected answer.
    out = struct('area', area, 'p', NaN, 'f', NaN);
    if isempty(res.stat) || ~isfield(res.stats, 'anySignificant') ...
            || ~res.stats.anySignificant || ~any(res.mask(:))
        return
    end
    pos = res.diff; pos(~res.mask) = NaN;
    [~, i] = max(abs(pos(:)));
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
