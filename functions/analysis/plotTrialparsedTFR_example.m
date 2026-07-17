function figFile = plotTrialparsedTFR_example(TFR, condition, opt, spec)
% plotTrialparsedTFR_example  Quick-look condition-mean spectrogram for NGL02_LFP.
%
% PURPOSE:
%   Renders one PNG per (area, alignment): the trial-mean, channel-mean
%   power spectrogram for a trial subset (default: `condition.correct == 1`),
%   stacked as one tiled panel per band from opt.freqInterest. Meant as
%   the trial-parsed sibling of plotContinuousTFR - a fast per-session
%   look so you can eyeball whether the LFP looks reasonable.
%
% USAGE:
%   figFile = plotTrialparsedTFR_example(TFR, condition, opt, spec);
%
% INPUTS:
%   TFR       - cell {nBands x 1} of ft_freqanalysis structs (per-band
%               output of computeTrialparsedTFR; each with
%               .powspctrm [nTrials x nChan x nFreq x nTime], .freq,
%               .time, .label).
%   condition - condition struct (from NGL02_postPhy). Consumed to build
%               the trial mask: default keeps trials with
%                 condition.(spec.trialFilter) == 1
%               where spec.trialFilter defaults to 'correct'.
%   opt       - resolved options. Consumes:
%                 .analysis            output root
%                 .SavFileName         session name (for filename)
%                 .toi                 [t0 t1] x-axis bounds (seconds)
%                 .lfp.plot.visible    figure Visible; default 'off'
%                 .lfp.plot.Resolution PNG DPI; default 300
%                 .lfp.plot.zlim       [zmin zmax] shared clim; [] -> auto
%                                       symmetric across ALL bands.
%                 .lfp.plot.colormap   colormap name/matrix; default hot
%                 .lfp.plot.interp     'bilinear' (default) | 'nearest' | 'none'
%   spec      - (optional) struct describing what to plot:
%                 .alignName    alignment tag (for title / filename)
%                 .areaTag      area tag       (for title / filename)
%                 .trialFilter  condition field forming the trial mask;
%                               default 'correct' (all trials with
%                               condition.correct == 1). Special value
%                               '' keeps every trial with valid TFR data.
%                 .baseline     [t0 t1] seconds -> baseline-corrected
%                               (dB) plot; default [-1 0] (or
%                               opt.lfp.plot.baseline via schema).
%                 .titlePrefix  extra prefix for the sgtitle
%
% OUTPUT:
%   figFile - absolute path to the PNG written.
%
% FILE:
%   <opt.analysis>/plots/TFR/<SavFileName>_<align>[_<area>]_TFR_<trialFilter>.png
%
% LAYOUT DETAILS (26.06.2026 polish):
%   * Panels stacked highest-freq on top, lowest on bottom.
%   * One shared symmetric color scale across ALL bands; one colorbar
%     centered on the tiledlayout's east side.
%   * One Freq (Hz) label centered on the tiledlayout's west side.
%   * X-axis extents forced to opt.toi([1 end]) exactly.
%   * Y-tick labels de-duplicated at shared tile boundaries.
%   * X-axis and X ticks shown only on the bottom tile so the stack
%     reads as a single figure.
%   * imagesc uses opt.lfp.plot.interp (default 'bilinear').
%
% NOTES:
%   * TFR{f} must carry per-trial data (keeptrials='yes' at compute time).
%     computeTrialparsedTFR ensures this.
%   * Trial-mask semantics: values of condition.(trialFilter) are
%     coerced to logical(v ~= 0). Length mismatch between the mask and
%     size(TFR{f}.powspctrm,1) -> warn and skip THIS band.
%
% SEE ALSO:
%   computeTrialparsedTFR (compute half), plotContinuousTFR (sibling
%   quick-look for the continuous mode), NGL02_LFP (caller).
%
% Last modified 26.06.2026 (Jesus) - shared color scale, tiledlayout east
%                                     colorbar / west ylabel, exact toi
%                                     x-range, tick-boundary de-duplication,
%                                     opt.lfp.plot.interp.

    if nargin < 4, spec = struct(); end
    trialFilter = getf(spec, 'trialFilter', 'correct');
    alignName   = getf(spec, 'alignName',   '');
    areaTag     = getf(spec, 'areaTag',     '');
    baseline    = getf(spec, 'baseline',    [-1 0]);
    titlePref   = getf(spec, 'titlePrefix', '');
    visible     = localOptField(opt, {'lfp','plot','visible'},    'off');
    resPNG      = localOptField(opt, {'lfp','plot','Resolution'}, 300);
    zLimUser    = localOptField(opt, {'lfp','plot','zlim'},       []);
    cmap        = localOptField(opt, {'lfp','plot','colormap'},   hot);
    interpMode  = localOptField(opt, {'lfp','plot','interp'},     'bilinear');

    assert(iscell(TFR) && ~isempty(TFR), 'NGL:plotTrialparsedTFR:badTFR', ...
        'TFR must be a non-empty cell of ft_freqanalysis structs.');

    nBands = numel(TFR);

    %% Trial mask.
    if isempty(trialFilter)
        makeMask = @(nTr) true(nTr, 1);
        maskName = 'allTrials';
    else
        if ~isstruct(condition) || ~isfield(condition, trialFilter)
            warning('NGL:plotTrialparsedTFR:noFilter', ...
                ['condition.%s not found; falling back to all trials for the ', ...
                 'quick-look plot.'], trialFilter);
            makeMask = @(nTr) true(nTr, 1);
            maskName = 'allTrials';
        else
            v = condition.(trialFilter);
            makeMask = @(nTr) localAlignMask(v, nTr);
            maskName = trialFilter;
        end
    end

    %% Panel ordering: highest-frequency band on TOP, lowest on bottom.
    bandFmin = nan(nBands, 1);
    for b = 1:nBands
        pow = TFR{b};
        if isstruct(pow) && isfield(pow, 'freq') && ~isempty(pow.freq)
            bandFmin(b) = min(pow.freq);
        end
    end
    [~, tileOrder] = sort(bandFmin, 'descend', 'MissingPlacement', 'last');

    %% PASS 1 - precompute per-band M so we can share a color scale
    %% across ALL bands (single colorbar at the end).
    Ms            = cell(nBands, 1);   % {slot} -> reduced [nFreq x nTime]
    nUsed         = zeros(nBands, 1);
    firstNonEmpty = 0;
    cbLbl         = 'dB';
    hasBaseline   = ~isempty(baseline) && numel(baseline) == 2;

    for slot = 1:nBands
        b   = tileOrder(slot);
        pow = TFR{b};
        if ~isstruct(pow) || ~isfield(pow, 'powspctrm') || ~isfield(pow, 'freq') || ~isfield(pow, 'time')
            continue
        end
        P = pow.powspctrm;
        if ndims(P) >= 4, nTr = size(P, 1); else, nTr = 1; end
        mask = makeMask(nTr);
        if numel(mask) ~= nTr
            warning('NGL:plotTrialparsedTFR:maskShape', ...
                'Trial mask (%d) does not match TFR{%d} trial dim (%d); skipping band.', ...
                numel(mask), b, nTr);
            continue
        end
        Pmask = P(mask, :, :, :);
        nUsed(b) = size(Pmask, 1);
        if nUsed(b) == 0, continue, end
        if firstNonEmpty == 0, firstNonEmpty = nUsed(b); end

        M = squeeze(mean(mean(Pmask, 1, 'omitnan'), 2, 'omitnan'));
        if hasBaseline
            tIdx = pow.time >= baseline(1) & pow.time <= baseline(2);
            if any(tIdx)
                base = mean(M(:, tIdx), 2, 'omitnan');
                M    = 10 * log10(M ./ max(base, eps));
                cbLbl = 'dB (re. baseline)';
            end
        end
        Ms{slot} = M;
    end

    %% Shared symmetric color scale across all bands.
    if ~isempty(zLimUser) && numel(zLimUser) == 2
        zLim = zLimUser;
    else
        allV = [];
        for slot = 1:nBands
            if isempty(Ms{slot}), continue, end
            v = Ms{slot}(:);
            v = v(isfinite(v));
            allV = [allV; v];
        end
        if isempty(allV)
            zLim = [-1 1];
        else
            m = max(abs(allV));
            zLim = [-m m];
        end
    end

    %% Figure + tiled layout.
    fig = figure('Visible', visible, 'Position', [50 50 900 max(180*nBands, 260)]);
    tl  = tiledlayout(fig, nBands, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    axs = gobjects(nBands, 1);

    for slot = 1:nBands
        b   = tileOrder(slot);
        ax  = nexttile(tl);
        axs(slot) = ax;
        pow = TFR{b};
        M   = Ms{slot};
        if isempty(M) || ~isstruct(pow) || ~isfield(pow, 'freq') || ~isfield(pow, 'time')
            axis(ax, 'off');
            continue
        end

        % imagesc with optional bilinear interpolation.
        if strcmpi(interpMode, 'bilinear')
            imagesc(ax, pow.time, pow.freq, M, 'Interpolation', 'bilinear');
        else
            imagesc(ax, pow.time, pow.freq, M);
        end
        set(ax, 'YDir', 'normal');

        % Exact axis extents. X = opt.toi (trims wavelet edges that
        % extend past the requested toi). Y = band bounds (prevents
        % MATLAB extending YLim to a "round" tick, which was causing
        % the shared-boundary duplicate freq numbers).
        if isfield(opt,'toi') && numel(opt.toi) >= 2
            ax.XLim = [opt.toi(1) opt.toi(end)];
        end
        ax.YLim = [min(pow.freq) max(pow.freq)];

        clim(ax, zLim);
        colormap(ax, cmap);

        xline(ax, 0, '--k', 'LineWidth', 2);
        box off

        % Bottom tile keeps the time axis + label; others hide theirs
        % so the stack reads as a single figure. Per-tile ylabel is
        % suppressed - the tiledlayout carries a single centered one.
        if slot < nBands
            ax.XTickLabel     = {};
            xlabel(ax, '');
            ax.XAxis.Visible  = 'off';
        else
            xlabel(ax, 'Time (s)');
        end
        ylabel(ax, '');
    end

    %% Dedupe freq ticks at shared tile boundaries. For each adjacent
    %% pair (upper = slot, lower = slot+1), if the upper's LOWEST tick
    %% equals the lower's HIGHEST tick (within 0.5 Hz tolerance), blank
    %% the lower tile's top label so the number is not printed twice.
    for slot = 1:(nBands - 1)
        up = axs(slot);       % higher-freq band (upper tile)
        lo = axs(slot + 1);   % lower-freq band  (lower tile)
        if ~isgraphics(up) || ~isgraphics(lo),   continue, end
        if isempty(up.YTick) || isempty(lo.YTick), continue, end
        if abs(up.YTick(1) - lo.YTick(end)) < 0.5
            lbl = lo.YTickLabel;
            if ~isempty(lbl)
                lbl{end} = '';
                lo.YTickLabel = lbl;
            end
        end
    end

    %% One shared colorbar centered on the east side of the tiledlayout.
    lastLive = find(arrayfun(@(a) isgraphics(a) && ~strcmp(a.Visible,'off'), axs), 1, 'last');
    if isempty(lastLive), lastLive = numel(axs); end
    cb = colorbar(axs(lastLive));
    cb.Layout.Tile   = 'east';
    cb.Label.String  = cbLbl;

    %% Single Freq (Hz) label centered on the west side of the tiledlayout.
    tl.YLabel.String     = 'Freq (Hz)';
    tl.YLabel.FontWeight = 'normal';

    %% Suptitle.
    stbits = {};
    if ~isempty(titlePref),          stbits{end+1} = titlePref;                end
    if isfield(opt,'SavFileName'),   stbits{end+1} = opt.SavFileName;          end
    if ~isempty(areaTag),            stbits{end+1} = ['area ' areaTag];        end
    if ~isempty(alignName),          stbits{end+1} = ['align ' alignName];     end
    if firstNonEmpty > 0
        stbits{end+1} = sprintf('n=%d %s', firstNonEmpty, maskName);
    end
    title(tl, strjoin(stbits, ' | '), 'Interpreter', 'none');

    %% Save.
    outDir = fullfile(opt.analysis, 'plots', 'TFR');
    if ~isfolder(outDir), mkdir(outDir); end
    stem = opt.SavFileName;
    if ~isempty(alignName), stem = [stem '_' alignName]; end
    if ~isempty(areaTag),   stem = [stem '_' areaTag];   end
    figFile = fullfile(outDir, [stem '_TFR_' trialFilter '.png']);
    exportgraphics(fig, figFile, 'Resolution', resPNG);
    close(fig);
    fprintf('plotTrialparsedTFR_example: wrote %s\n', figFile);
end


% =======================================================================
function v = getf(s, f, dflt)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function m = localAlignMask(v, nTr)
% Coerce a per-trial condition vector to a length-nTr logical mask.
    v = v(:);
    if numel(v) == nTr
        m = logical(v ~= 0);
    elseif numel(v) > nTr
        m = logical(v(1:nTr) ~= 0);
    else
        m = false(nTr, 1);
        m(1:numel(v)) = logical(v ~= 0);
    end
end

function v = localOptField(opt, path, dflt)
    v = dflt;  cursor = opt;
    for k = 1:numel(path)
        if isstruct(cursor) && isfield(cursor, path{k})
            cursor = cursor.(path{k});
        else
            return
        end
    end
    v = cursor;
end
