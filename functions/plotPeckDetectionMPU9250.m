function h = plotPeckDetectionMPU9250(data, fs, out, opts)
% plotPeckDetectionMPU9250  Plot raw + filtered IMU signals, thresholds, detections,
%                           a SMALL sample of candidates, and feature visualizations.
%
% h = plotPeckDetectionMPU9250(data, fs, out)
% h = plotPeckDetectionMPU9250(data, fs, out, opts)
%
% INPUTS
%   data: struct used as input to detectPecksMPU9250 (raw axes in data.acc.*, data.gyr.*)
%   fs:   sampling rate (Hz)
%   out:  output struct from detectPecksMPU9250 (expects out.debug.signals + out.debug.thresholds)
%   opts: (optional) struct
%       .gt = []                         % (Nx1) logical/0-1. Marker that happens AFTER a true peck
%       .titlePrefix = ''
%       .showCandidates = true
%
%       % Marker downsampling (keep plots readable)
%       .maxCandMarkersOverview = 40     % max candidate markers drawn in overview
%       .maxGTMarkersOverview   = 40     % max GT markers drawn in overview
%       .maxDetMarkersOverview  = 200    % max detected markers drawn in overview
%       .maxCandMarkersPerEvent = 10     % max candidate markers in each zoomed event panel
%
%       % Event panels
%       .maxEventPanels = 6
%       .eventWinSec = 0.8
%
%       % Feature visualization for the sampled candidates
%       .featureWinSec = 0.25            % window around candidate for feature extraction
%       .maxFeatureCandidates = 15       % number of candidates shown in feature plots (subset of overview sample)
%
% OUTPUT
%   h.overview : overview figure
%   h.events   : zoomed detection panels figure (empty if no detections)
%   h.features : features figure (empty if candidates not shown)
%
% Notes
% - Candidate markers are sampled by highest composite score among candidates.
% - Feature visualization uses windowed peaks of |jerk|, |gyro|, hfRMS and jerk half-max duration.

if nargin < 4 || isempty(opts), opts = struct(); end
opts = setDefault(opts, 'gt', []);
opts = setDefault(opts, 'titlePrefix', '');
opts = setDefault(opts, 'showCandidates', true);

opts = setDefault(opts, 'maxCandMarkersOverview', 40);
opts = setDefault(opts, 'maxGTMarkersOverview', 40);
opts = setDefault(opts, 'maxDetMarkersOverview', 200);
opts = setDefault(opts, 'maxCandMarkersPerEvent', 10);

opts = setDefault(opts, 'maxEventPanels', 6);
opts = setDefault(opts, 'eventWinSec', 0.8);

opts = setDefault(opts, 'featureWinSec', 0.25);
opts = setDefault(opts, 'maxFeatureCandidates', 15);

% ---------- Validate output struct ----------
if ~isfield(out,'debug') || ~isfield(out.debug,'signals') || ~isfield(out.debug,'thresholds')
    error('out.debug.signals and out.debug.thresholds are required. Pass the output of detectPecksMPU9250.');
end

sig = out.debug.signals;
thr = out.debug.thresholds;

t = sig.t(:);
N = numel(t);

% Raw accel magnitude from input (pre-alignment)
axr = data.acc.X(:); ayr = data.acc.Y(:); azr = data.acc.Z(:);
aMagRaw = sqrt(axr.^2 + ayr.^2 + azr.^2);

% Aligned signals from detector
jerk = sig.jerk(:);
gMag = sig.gMag(:);
hfRMS = sig.hfRMS(:);

scoreSignal = zeros(N,1);
if isfield(sig,'scoreSignal') && numel(sig.scoreSignal)==N
    scoreSignal = sig.scoreSignal(:);
end

% Candidates (same logical used by detector)
isCandidate = (abs(jerk) > thr.jerk) & (gMag > thr.gyro) & (hfRMS > thr.hfEnergy);
candIdxAll = find(isCandidate);

% Downsampled candidate indices for display
candIdx = candIdxAll;
if opts.showCandidates && ~isempty(candIdxAll)
    candIdx = sampleTopByScore(candIdxAll, scoreSignal, opts.maxCandMarkersOverview);
else
    candIdx = [];
end

% Detected events (may be many; cap for overview)
detIdxAll = [];
if isfield(out,'idx') && ~isempty(out.idx)
    detIdxAll = out.idx(:);
end
detIdx = detIdxAll;
if numel(detIdxAll) > opts.maxDetMarkersOverview
    detIdx = detIdxAll( round(linspace(1, numel(detIdxAll), opts.maxDetMarkersOverview)) );
end

% Ground truth markers (optional; cap)
gtIdxAll = [];
if ~isempty(opts.gt)
    gt = logical(opts.gt(:));
    gt = gt(1:min(end,N));
    gtIdxAll = find(gt);
end
gtIdx = gtIdxAll;
if numel(gtIdxAll) > opts.maxGTMarkersOverview
    gtIdx = gtIdxAll( round(linspace(1, numel(gtIdxAll), opts.maxGTMarkersOverview)) );
end

% ---------- Marker colors ----------
mk = struct();
mk.detColor  = [0 0 0];
mk.candColor = [0.6 0.6 0.6];
mk.gtColor   = [0.85 0 0];

titlePrefix = opts.titlePrefix;
if ~isempty(titlePrefix) && titlePrefix(end) ~= ' '
    titlePrefix = [titlePrefix ' '];
end

% ======================= Figure 1: Overview =======================
h = struct();
h.overview = figure('Name','Peck detection overview','Color','w');
tiledlayout(6,1,'TileSpacing','compact','Padding','compact');

% 1) Raw accel magnitude
nexttile;
plot(t, aMagRaw, 'LineWidth', 1); grid on;
ylabel('|a| raw');
title([titlePrefix 'Raw accel magnitude (input axes)']);
hold on; plotEventMarkers(t, detIdx, candIdx, gtIdx, mk, opts.showCandidates); hold off;

% 2) Aligned accel magnitude
nexttile;
plot(t, sig.aMag(:), 'LineWidth', 1); grid on;
ylabel('|a| aligned');
title('Accel magnitude (aligned frame, includes gravity)');
hold on; plotEventMarkers(t, detIdx, candIdx, gtIdx, mk, opts.showCandidates); hold off;

% 3) High-pass accel magnitude
nexttile;
plot(t, sig.aHPmag(:), 'LineWidth', 1); grid on;
ylabel('aHP');
title('High-pass accel magnitude');
hold on; plotEventMarkers(t, detIdx, candIdx, gtIdx, mk, opts.showCandidates); hold off;

% 4) Bandpassed RMS (hfRMS) + threshold
nexttile;
plot(t, hfRMS, 'LineWidth', 1); grid on;
yline(thr.hfEnergy, '--', 'thr.hf', 'LabelHorizontalAlignment','left');
ylabel('hfRMS');
title('Bandpassed accel RMS (impact energy proxy)');
hold on; plotEventMarkers(t, detIdx, candIdx, gtIdx, mk, opts.showCandidates); hold off;

% 5) |jerk| + threshold
nexttile;
plot(t, abs(jerk), 'LineWidth', 1); grid on;
yline(thr.jerk, '--', 'thr.jerk', 'LabelHorizontalAlignment','left');
ylabel('|jerk|');
title('Jerk magnitude (|d/dt aHP|)');
hold on; plotEventMarkers(t, detIdx, candIdx, gtIdx, mk, opts.showCandidates); hold off;

% 6) |gyro| + score + gyro threshold
nexttile;
yyaxis left;
plot(t, gMag, 'LineWidth', 1); grid on;
yline(thr.gyro, '--', 'thr.gyro', 'LabelHorizontalAlignment','left');
ylabel('|gyro|');
yyaxis right;
plot(t, scoreSignal, 'LineWidth', 1);
ylabel('score');
title('Gyro magnitude and composite score');
hold on; plotEventMarkers(t, detIdx, candIdx, gtIdx, mk, opts.showCandidates); hold off;

ax = gca;
legend(ax, composeLegend(opts.showCandidates && ~isempty(candIdx), ~isempty(gtIdx), ~isempty(detIdx)), 'Location','best');

% ======================= Figure 2: Zoomed event panels =======================
h.events = [];
K = min(opts.maxEventPanels, numel(detIdxAll));
if K > 0
    h.events = figure('Name','Peck detection event panels','Color','w');
    tiledlayout(K,1,'TileSpacing','compact','Padding','compact');

    halfW = round((opts.eventWinSec/2) * fs);

    % pick K detections spread across time (not necessarily the first K)
    pick = unique(round(linspace(1, numel(detIdxAll), K)));
    detForPanels = detIdxAll(pick);

    for k = 1:K
        i0 = detForPanels(k);
        i1 = max(1, i0 - halfW);
        i2 = min(N, i0 + halfW);

        nexttile;
        plot(t(i1:i2), abs(jerk(i1:i2)), 'LineWidth', 1); grid on; hold on;
        yline(thr.jerk, '--');

        plot(t(i1:i2), hfRMS(i1:i2), 'LineWidth', 1);
        yline(thr.hfEnergy, '--');

        plot(t(i1:i2), gMag(i1:i2), 'LineWidth', 1);
        yline(thr.gyro, '--');

        xline(t(i0), '-', 'Detected', 'LabelVerticalAlignment','bottom', ...
            'Color', mk.detColor, 'LineWidth', 1.2);

        % candidates inside window: sample a few by score
        if opts.showCandidates && ~isempty(candIdxAll)
            cwin = candIdxAll(candIdxAll>=i1 & candIdxAll<=i2);
            if ~isempty(cwin)
                cwinS = sampleTopByScore(cwin, scoreSignal, opts.maxCandMarkersPerEvent);
                xline(t(cwinS), ':', 'Color', mk.candColor);
            end
        end

        % GT inside window (do not explode)
        if ~isempty(gtIdxAll)
            gwin = gtIdxAll(gtIdxAll>=i1 & gtIdxAll<=i2);
            if ~isempty(gwin)
                gwin = gwin(1:min(end, 10)); % hard cap
                xline(t(gwin), '-', 'GT', 'Color', mk.gtColor, 'LineWidth', 1.0);
            end
        end

        % annotate with detector features if present
        s = '';
        if isfield(out,'features') && ~isempty(out.features) && any(out.features.idx == i0)
            row = out.features(out.features.idx == i0, :);
            s = sprintf(' jerkPeak=%.3g, gyroPeak=%.3g, hfRMSPeak=%.3g, durHalf=%.3gs', ...
                row.jerkPeak, row.gyroPeak, row.hfRMSPeak, row.durAboveHalf);
        end
        title(sprintf('Detection @ t=%.3fs.%s', t(i0), s));
        if k == K, xlabel('time (s)'); end
        hold off;
    end
end

% ======================= Figure 3: Feature visualization (sample candidates) =======================
h.features = [];
if opts.showCandidates && ~isempty(candIdxAll)
    % Choose feature candidates: top by score, limited
    candFeatIdx = sampleTopByScore(candIdxAll, scoreSignal, min(opts.maxFeatureCandidates, numel(candIdxAll)));

    % Compute features for these candidates from aligned signals
    halfWF = round((opts.featureWinSec/2) * fs);
    featCand = computeCandidateFeatures(candFeatIdx, halfWF, jerk, gMag, hfRMS, fs);

    % Also show detections as reference (limited)
    detFeatIdx = detIdxAll;
    if numel(detFeatIdx) > 40
        detFeatIdx = detIdxAll(round(linspace(1, numel(detIdxAll), 40)));
    end
    featDet = computeCandidateFeatures(detFeatIdx, halfWF, jerk, gMag, hfRMS, fs);

    h.features = figure('Name','Candidate feature visualization','Color','w');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % (1) jerkPeak vs gyroPeak with thresholds
    nexttile;
    scatter(featDet.jerkPeak, featDet.gyroPeak, 'filled'); grid on; hold on;
    scatter(featCand.jerkPeak, featCand.gyroPeak, 'o', 'LineWidth', 1.2);
    xline(thr.jerk, '--', 'thr.jerk');
    yline(thr.gyro, '--', 'thr.gyro');
    xlabel('jerkPeak');
    ylabel('gyroPeak');
    title('Feature space: jerkPeak vs gyroPeak');
    legend({'detections (sample)','candidates (sample)'}, 'Location','best');
    hold off;

    % (2) hfRMSPeak vs jerkPeak with thresholds
    nexttile;
    scatter(featDet.jerkPeak, featDet.hfRMSPeak, 'filled'); grid on; hold on;
    scatter(featCand.jerkPeak, featCand.hfRMSPeak, 'o', 'LineWidth', 1.2);
    xline(thr.jerk, '--', 'thr.jerk');
    yline(thr.hfEnergy, '--', 'thr.hf');
    xlabel('jerkPeak');
    ylabel('hfRMSPeak');
    title('Feature space: jerkPeak vs hfRMSPeak');
    legend({'detections (sample)','candidates (sample)'}, 'Location','best');
    hold off;

    % (3) Normalized feature bars for candidate sample
    nexttile([1 2]);
    X = [featCand.jerkPeak featCand.gyroPeak featCand.hfRMSPeak featCand.durAboveHalf];
    Xn = normalizeCols01(X);
    bar(Xn); grid on;
    xlabel('candidate # (sorted by score)');
    ylabel('normalized feature value');
    title('Candidate sample: normalized feature vectors');
    legend({'jerkPeak','gyroPeak','hfRMSPeak','durAboveHalf'}, 'Location','best');

    % Sort the bars by score (descending)
    [~,ord] = sort(featCand.score,'descend');
    reorderBarOrd(ord);
end

end

% =================== Local helper functions ===================

function s = setDefault(s, field, value)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = value;
end
end

function idxS = sampleTopByScore(idx, scoreSignal, K)
% Take top-K indices by scoreSignal, preserving time order.
if isempty(idx), idxS = idx; return; end
K = min(K, numel(idx));
sc = scoreSignal(idx);
[~,ord] = sort(sc, 'descend');
idxS = sort(idx(ord(1:K)));
end

function plotEventMarkers(t, detIdx, candIdx, gtIdx, mk, showCandidates)
yl = ylim;

if showCandidates && ~isempty(candIdx)
    x = t(candIdx);
    for i = 1:numel(x)
        line([x(i) x(i)], yl, 'Color', mk.candColor, 'LineStyle', ':', 'LineWidth', 0.8);
    end
end

if ~isempty(detIdx)
    x = t(detIdx);
    for i = 1:numel(x)
        line([x(i) x(i)], yl, 'Color', mk.detColor, 'LineStyle', '-', 'LineWidth', 1.2);
    end
end

if ~isempty(gtIdx)
    x = t(gtIdx);
    for i = 1:numel(x)
        line([x(i) x(i)], yl, 'Color', mk.gtColor, 'LineStyle', '-', 'LineWidth', 1.0);
    end
end

ylim(yl);
end

function L = composeLegend(hasCand, hasGT, hasDet)
L = {};
L{end+1} = 'signal';
if hasCand, L{end+1} = 'candidate sample'; end
if hasDet,  L{end+1} = 'detected'; end
if hasGT,   L{end+1} = 'ground truth sample'; end
end

function feat = computeCandidateFeatures(idx, halfW, jerk, gMag, hfRMS, fs)
% Windowed feature extraction around indices, using already-aligned signals.
n = numel(idx);
jerkPeak = zeros(n,1);
gyroPeak = zeros(n,1);
hfRMSPeak = zeros(n,1);
durHalf = zeros(n,1);
hfEnergyApprox = zeros(n,1);
score = zeros(n,1);

N = numel(jerk);
for k = 1:n
    i0 = idx(k);
    i1 = max(1, i0-halfW);
    i2 = min(N, i0+halfW);

    jw = abs(jerk(i1:i2));
    gw = gMag(i1:i2);
    hw = hfRMS(i1:i2);

    jerkPeak(k) = max(jw);
    gyroPeak(k) = max(gw);
    hfRMSPeak(k)= max(hw);

    th = 0.5*max(jw);
    durHalf(k) = sum(jw > th)/fs;

    hfEnergyApprox(k) = sum(hw.^2)/fs;

    % score proxy (same structure as detector)
    score(k) = mean(jw) + mean(gw) + mean(hw);
end

feat = table(idx(:), jerkPeak, gyroPeak, hfRMSPeak, durHalf, hfEnergyApprox, score, ...
    'VariableNames', {'idx','jerkPeak','gyroPeak','hfRMSPeak','durAboveHalf','hfEnergyApprox','score'});
end

function Xn = normalizeCols01(X)
Xn = X;
for j = 1:size(X,2)
    col = X(:,j);
    col = col - min(col);
    d = max(col);
    if d > 0
        Xn(:,j) = col/d;
    else
        Xn(:,j) = zeros(size(col));
    end
end
end

function reorderBarOrd(ord)
% Reorder bar groups by ord (works by reordering XData/YData on current axes)
ax = gca;
ch = ax.Children;
% Find Bar object (usually first child, but not guaranteed)
barObj = [];
for k = 1:numel(ch)
    if isa(ch(k),'matlab.graphics.chart.primitive.Bar')
        barObj = ch(k);
        break;
    end
end
if isempty(barObj), return; end
Y = barObj.YData;
barObj.YData = Y(:,ord);
end
