function [lfp, spike, motion, truth] = makeFixtureSignals(blocks, geom, duration, cfg)
%MAKEFIXTURESIGNALS  The signals of a synthetic session, and what was planted in them.
%
% PURPOSE:
%   Everything a test can check lives here: what goes into the LFP, the spike
%   trains and the IMU, and nothing else does. The returned `truth` says in
%   words what each analysis should find, so a test asserts against the design
%   rather than against a previous run's output.
%
% USAGE:
%   [lfp, spike, motion, truth] = makeFixtureSignals(blocks, geom, duration, cfg)
%   Called by makeFakeSession; not meant to be called directly.
%
% INPUTS:
%   blocks   - struct array from makeFakeSession's block builder
%   geom     - channel geometry (label, ycoords, kcoords, area)
%   duration - session length in seconds
%   cfg      - fixture configuration (fs, bands, codes)
%
% OUTPUT:
%   lfp    - [nChan x nSamples] single, at cfg.fs
%   spike  - the struct NGL02_postPhy writes after Kilosort + Phy
%   motion - MotionData_raw contents (3-axis accelerometer)
%   truth  - what was planted: effects, units, pecks, and notes
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).

    [lfp, truth.lfp] = localBuildLFP(blocks, geom, duration, cfg);
    [spike, truth.spikes] = localBuildSpikes(blocks, geom, duration, cfg);
    [motion, truth.imu] = localBuildIMU(blocks, duration);
end

% ======================= LFP =======================
function [lfp, truth] = localBuildLFP(blocks, geom, duration, cfg)
    fs = cfg.fs;
    n = round(duration * fs);
    t = (0:n-1) / fs;
    nCh = numel(geom.label);
    isNCL = strcmp(geom.area, 'NCL');

    % 1/f-ish background, partly shared across channels (volume conduction).
    common = localPink(1, n);
    lfp = 25 * localPink(nCh, n) + 18 * repmat(common, nCh, 1);

    truth = struct();
    truth.notes = {};

    % A dead channel and a noisy one: real recordings have them, and code that
    % cannot cope with them should fail here rather than on real data.
    truth.deadChannel  = geom.label{7};
    truth.noisyChannel = geom.label{20};
    lfp(7, :) = 0.5 * randn(1, n);
    lfp(20, :) = lfp(20, :) + 120 * randn(1, n);

    % Laminar evoked response with a KNOWN current source density: a sink at
    % 325 um with flanking sources. The potential is its double spatial
    % integral, so computeCSD recovers what was planted.
    sinkDepth = 325; sigma = 0.3;
    zf = (-800:2:1600)';
    csdTrue = -exp(-0.5*((zf - sinkDepth)/55).^2) ...
              + 0.5*exp(-0.5*((zf - (sinkDepth-170))/85).^2) ...
              + 0.5*exp(-0.5*((zf - (sinkDepth+170))/85).^2);
    csdTrue = csdTrue - mean(csdTrue);
    phiF = cumtrapz(cumtrapz(-csdTrue / sigma)) * median(diff(zf))^2;
    phiF = phiF - interp1([zf(1) zf(end)], [phiF(1) phiF(end)], zf);
    phiF = phiF / max(abs(phiF));
    phiByDepth = interp1(zf, phiF, geom.ycoords, 'linear', 0);
    phiByDepth(geom.kcoords ~= 1) = 0;             % shank 1 only
    truth.csdSinkDepth = sinkDepth;
    truth.csdShank = 1;
    truth.notes{end+1} = sprintf(['CSD: sink at %g um on shank 1, ~60 ms after ' ...
        'stimOn1 of dms and arena trials'], sinkDepth);

    truth.betaBand = cfg.bands.beta;
    truth.thetaBand = cfg.bands.theta;
    truth.tagRates = [];
    for b = 1:numel(blocks)
        trials = blocks(b).trials;
        switch blocks(b).task
            case 'dms'
                for k = 1:numel(trials)
                    tr = trials(k);
                    lfp = localAddEvoked(lfp, t, fs, tr.tStim1 + 0.06, phiByDepth * 140);
                    if strcmp(tr.outcome, 'correct')
                        % Beta only on correct trials and only in NCL: the
                        % contrast must find it there and nowhere else.
                        lfp = localAddBurst(lfp, t, fs, tr.tStim2 + 0.35, 0.22, ...
                                            20, 55, isNCL);
                    end
                end
                truth.notes{end+1} = ['dms: 20 Hz burst 0.35 s after stimOn2 on ' ...
                    'CORRECT trials, NCL only (STR has none)'];

            case 'arena'
                for k = 1:numel(trials)
                    tr = trials(k);
                    lfp = localAddEvoked(lfp, t, fs, tr.tStim1 + 0.06, phiByDepth * 110);
                    lfp = localAddBurst(lfp, t, fs, (tr.tStim1 + tr.tStim2)/2, ...
                                        max(tr.travel/2, 0.2), 6, 40, true(nCh,1));
                    if ~isempty(tr.pecks)
                        lfp = localAddArtifact(lfp, t, fs, tr.pecks(1) + 0.02);
                    end
                end
                truth.notes{end+1} = ['arena: 6 Hz theta while walking, both ' ...
                    'areas; a movement artifact at each peck'];

            case 'nft'
                rates = unique([trials.rate]);
                truth.tagRates = rates;
                gain = 30 * double(isNCL) .* (0.6 + 0.4 * double(geom.kcoords == 1));
                for k = 1:numel(trials)
                    lfp = localAddEvoked(lfp, t, fs, trials(k).tStim1, gain);
                end
                truth.notes{end+1} = sprintf(['nft: stimulus-locked response at ' ...
                    '%s Hz in NCL; STR not driven'], mat2str(rates));
        end
    end
    lfp = single(lfp);
end

function x = localPink(nCh, n)
% Cheap 1/f-ish noise: white noise integrated, slowest drift removed.
    x = cumsum(randn(nCh, n), 2);
    x = x - movmean(x, max(3, round(n/20)), 2);
    x = x ./ std(x, 0, 2);
end

function lfp = localAddEvoked(lfp, t, fs, tEvent, ampByChan)
    w = round(0.12 * fs);
    i0 = round((tEvent - t(1)) * fs);
    if i0 < 1 || i0 + w > size(lfp, 2), return; end
    shape = sin(pi * (0:w-1) / w).^2;
    lfp(:, i0:i0+w-1) = lfp(:, i0:i0+w-1) + ampByChan(:) * shape;
end

function lfp = localAddBurst(lfp, t, fs, tCentre, halfWidth, freq, amp, chanMask)
    i0 = max(1, round((tCentre - 3*halfWidth - t(1)) * fs));
    i1 = min(size(lfp, 2), round((tCentre + 3*halfWidth - t(1)) * fs));
    if i1 <= i0, return; end
    tt = t(i0:i1);
    env = exp(-0.5 * ((tt - tCentre) / max(halfWidth, 1e-3)).^2);
    wave = amp * sin(2*pi*freq*tt + 2*pi*rand) .* env;
    gain = chanMask(:) .* (0.7 + 0.6 * rand(size(lfp, 1), 1));
    lfp(:, i0:i1) = lfp(:, i0:i1) + gain * wave;
end

function lfp = localAddArtifact(lfp, t, fs, tAt)
% A movement artifact: large, brief, on every channel at once.
    w = round(0.03 * fs);
    i0 = round((tAt - t(1)) * fs);
    if i0 < 1 || i0 + w > size(lfp, 2), return; end
    shape = (1 - cos(2*pi*(0:w-1)/w)) / 2;
    lfp(:, i0:i0+w-1) = lfp(:, i0:i0+w-1) + (300 + 150*rand) * shape;
end

% ======================= spikes =======================
function [spike, truth] = localBuildSpikes(blocks, geom, duration, cfg) %#ok<INUSD>
    units = localUnitPlan(geom);
    nU = numel(units);
    spike = struct();
    [spike.label, spike.timestamp, spike.depth, spike.ampl, spike.templampl, ...
     spike.ch, spike.shank, spike.roi, spike.KSLabel, spike.bc_unitType, ...
     spike.HumanLabel, spike.phyLabel] = deal(cell(1, nU));

    for u = 1:nU
        U = units(u);
        st = localPoisson(U.baseRate, duration);
        for b = 1:numel(blocks)
            for k = 1:numel(blocks(b).trials)
                st = [st; localDrivenSpikes(U, blocks(b).trials(k), blocks(b).task)]; %#ok<AGROW>
            end
        end
        st = sort(st(st > 0 & st < duration));
        st(diff([0; st]) < 0.0015) = [];            % refractory period
        spike.label{u}      = num2str(U.id);
        spike.timestamp{u}  = st;
        spike.depth{u}      = repmat(geom.ycoords(U.chan), numel(st), 1);
        amp = U.amp * ones(numel(st), 1);
        if U.drifts
            amp = amp .* (1 - 0.45 * st / duration);   % electrode drift
        end
        spike.ampl{u}       = amp;
        spike.templampl{u}  = amp / U.amp;
        spike.ch{u}         = U.chan;
        spike.shank{u}      = geom.kcoords(U.chan);
        spike.roi{u}        = geom.area{U.chan};
        spike.KSLabel{u}    = U.ksLabel;
        spike.bc_unitType{u}= U.bcType;
        spike.HumanLabel{u} = U.human;
        spike.phyLabel{u}   = U.human;
    end

    truth = struct();
    truth.nUnits       = nU;
    truth.ids          = [units.id];
    truth.drivenUnits  = [units([units.stim1Gain] > 0).id];
    truth.rewardUnit   = units(find([units.rewardGain] > 0, 1)).id;
    truth.taggedUnit   = units(find([units.tagLocked], 1)).id;
    truth.driftingUnit = units(find([units.drifts], 1)).id;
    truth.noiseUnit    = units(find(strcmp({units.human}, 'noise'), 1)).id;
    truth.notes = {'driven units fire 40-70 ms after stimOn1', ...
                   'one unit follows reward, one is locked to the tagging stream', ...
                   'one unit loses ~45% of its amplitude across the session', ...
                   'one cluster is labelled noise, as Phy curation would leave it'};
end

function units = localUnitPlan(geom)
    shank1 = find(geom.kcoords == 1);
    shank2 = find(geom.kcoords == 2);
    mk = @(id, ch, rate, amp, s1, rw, tag, drift, ks, bc, hum) struct( ...
        'id', id, 'chan', ch, 'baseRate', rate, 'amp', amp, 'stim1Gain', s1, ...
        'rewardGain', rw, 'tagLocked', tag, 'drifts', drift, 'ksLabel', ks, ...
        'bcType', bc, 'human', hum);
    units = [ ...
        mk(101, shank1(6),  4.0,  90, 12, 0, false, false, 'good', 'GOOD',  'good'), ...
        mk(102, shank1(7),  6.5, 120, 18, 0, false, false, 'good', 'GOOD',  'good'), ...
        mk(103, shank1(11), 2.0,  70,  0, 9, false, false, 'good', 'GOOD',  'good'), ...
        mk(104, shank1(9),  3.0,  80,  6, 0, true,  false, 'good', 'GOOD',  'good'), ...
        mk(105, shank1(13), 5.0, 100,  0, 0, false, true,  'good', 'MUA',   'mua'), ...
        mk(201, shank2(5),  8.0,  60,  4, 0, false, false, 'good', 'GOOD',  'good'), ...
        mk(202, shank2(10), 1.5,  55,  0, 5, false, false, 'mua',  'MUA',   'mua'), ...
        mk(203, shank2(14), 22.0, 30,  0, 0, false, false, 'mua',  'NOISE', 'noise')];
end

function st = localPoisson(rate, duration)
    st = cumsum(-log(rand(ceil(rate * duration * 1.4) + 20, 1)) / rate);
    st = st(st < duration);
end

function st = localDrivenSpikes(U, tr, task)
% The extra spikes events drive, on top of the background rate.
    st = zeros(0, 1);
    if U.stim1Gain > 0 && ~isnan(tr.tStim1)
        n = localPoissonCount(U.stim1Gain * 0.05);
        st = [st; tr.tStim1 + 0.04 + 0.03 * rand(n, 1)];
    end
    if U.rewardGain > 0 && strcmp(tr.outcome, 'correct')
        n = localPoissonCount(U.rewardGain * 0.15);
        st = [st; tr.tEnd - 0.2 + 0.15 * rand(n, 1)];
    end
    if U.tagLocked && strcmp(task, 'nft')
        n = localPoissonCount(3);
        st = [st; tr.tStim1 + 0.03 + 0.02 * randn(n, 1)];
    end
end

function n = localPoissonCount(lambda)
% Knuth's method: avoids a Statistics Toolbox dependency for a few counts.
    L = exp(-lambda); k = 0; p = 1;
    while true
        p = p * rand;
        if p <= L, break; end
        k = k + 1;
        if k > 50 + 10 * lambda, break; end
    end
    n = k;
end

% ======================= IMU =======================
function [motion, truth] = localBuildIMU(blocks, duration)
% Three-axis accelerometer in the shape GetMotionSensors writes for INTAN:
% raw.acc.X/.Y/.Z in m/s^2 at 1 kHz, no gyroscope or magnetometer (the
% ADXL335 has neither). Pecks are the sharp events detect_jerkEvents is built
% to find; walking is slow sway.
    fs = 1000;                       % par.targetFs in GetMotionSensors
    g0 = 9.81;
    n = round(duration * fs);
    t = (0:n-1) / fs;
    acc = 0.15 * randn(3, n);        % m/s^2 of sensor noise
    acc(3, :) = acc(3, :) + g0;      % gravity on z, as when the head is level

    peckTimes = [];
    for b = 1:numel(blocks)
        for k = 1:numel(blocks(b).trials)
            tr = blocks(b).trials(k);
            for pk = tr.pecks(:)'
                acc = localJerk(acc, t, fs, pk, g0);
                peckTimes(end+1) = pk; %#ok<AGROW>
            end
            if strcmp(blocks(b).task, 'arena') && ~isnan(tr.tStim2)
                acc = localSway(acc, t, fs, tr.tStim1, tr.tStim2, g0);
            end
        end
    end

    raw = struct();
    raw.acc = struct('X', acc(1, :), 'Y', acc(2, :), 'Z', acc(3, :));
    raw.gyr = []; raw.mag = [];
    raw.tsec = t;
    raw.fs = fs;
    raw.source = 'fixture-adxl335-like';
    raw.dof = 3;
    raw.chipFrame = struct('acc', 'Z-up when level, as an ADXL335 mounted flat.', ...
                           'note', 'Synthetic: no mounting rotation needed.');
    raw.meta = struct('fixture', true, 'units', 'm/s^2');
    motion = struct('raw', raw);

    truth = struct('fs', fs, 'nPecks', numel(peckTimes), ...
                   'peckTimes', sort(peckTimes(:)));
    truth.notes = {['pecks are sharp 30 ms jerks; detect_jerkEvents on the ' ...
                    'jerk magnitude should find them'], ...
                   'arena trials add slow sway while the animal walks'};
end

function acc = localJerk(acc, t, fs, tAt, g0)
% A peck: one sharp cycle, ~1 g on the forward axis and half that on z. Short
% enough that the jerk (the derivative) is what stands out, which is what the
% detector works on.
    w = round(0.03 * fs);
    i0 = round((tAt - t(1)) * fs);
    if i0 < 1 || i0 + w > size(acc, 2), return; end
    shape = sin(2*pi*(0:w-1)/w);
    acc(1, i0:i0+w-1) = acc(1, i0:i0+w-1) + (0.9 + 0.4*rand) * g0 * shape;
    acc(3, i0:i0+w-1) = acc(3, i0:i0+w-1) + (0.5 + 0.3*rand) * g0 * shape;
end

function acc = localSway(acc, t, fs, t0, t1, g0)
% Walking: slow, large and smooth, so a detector tuned to sharp events must
% ignore it. If it does not, that is the fixture doing its job.
    i0 = max(1, round((t0 - t(1)) * fs));
    i1 = min(size(acc, 2), round((t1 - t(1)) * fs));
    if i1 <= i0, return; end
    tt = t(i0:i1);
    % Ramp in and out over 200 ms. Starting the sway as a step would put a
    % discontinuity in the acceleration, and a jerk detector would be right to
    % call that an event - an animal does not begin walking instantaneously.
    ramp = ones(size(tt));
    w = min(round(0.2 * fs), floor(numel(tt) / 2));
    if w > 1
        edge = (1 - cos(pi * (0:w-1) / (w-1))) / 2;
        ramp(1:w) = edge;
        ramp(end-w+1:end) = fliplr(edge);
    end
    acc(2, i0:i1) = acc(2, i0:i1) + 0.25 * g0 * ramp .* sin(2*pi*2.5*tt + 2*pi*rand);
end
