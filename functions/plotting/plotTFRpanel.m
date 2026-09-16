function h = plotTFRpanel(ax, t, f, M, st, varargin)
% plotTFRpanel  One time x frequency panel: the unit every LFP figure is built from.
%
% PURPOSE:
%   Power maps, contrasts, comodulograms and tagging spectra are all a matrix
%   over two axes with an event at zero and, often, a significance mask. Drawn
%   here once, they share colour limits, interpolation, the event line and the
%   way significance is marked, so panels from different analyses can be put
%   in one figure without re-tuning anything.
%
% USAGE:
%   plotTFRpanel(ax, freq.time, freq.freq, squeeze(mean(P, 1)), lfpStyle(opt))
%   plotTFRpanel(ax, t, f, M, st, 'mask', stat.mask, 'title', 'NCL correct')
%
% INPUTS:
%   ax  - target axes (gca if []).
%   t   - time vector, seconds, length nT.
%   f   - frequency vector, Hz, length nF.
%   M   - [nF x nT] matrix to draw (already normalised; see normalizeTFR).
%   st  - style struct from lfpStyle.
%   Name/value pairs:
%     'mask'      [nF x nT] logical; significant bins. Drawn as an outline,
%                 with everything outside it faded (see 'maskStyle').
%     'maskStyle' 'outline' (default) | 'fade' | 'both' | 'none'
%     'title'     panel title
%     'xlabel' / 'ylabel'   axis labels ('' to suppress, for shared labels)
%     'event'     time(s) of vertical marker lines (default 0; [] for none)
%     'clim'      explicit [lo hi], overriding the style's rule
%     'logf'      true = logarithmic frequency axis (default false)
%     'ydir'      'normal' (default) or 'reverse', for a depth axis that
%                 should run downwards
%     'colorbar'  true = attach one to this panel (default false; a shared
%                 colorbar per figure is usually what you want)
%
% OUTPUT:
%   h - struct with .image, .maskLine, .colorbar handles for further tweaking.
%
% NOTES:
%   * Colour limits: 'clim' if given, else st.zlim, else the data's robust
%     range ('climPercentile', 2nd to 98th by default) so a single artifact
%     bin cannot flatten the whole map. Diverging styles force the range
%     symmetric about zero. Widen the percentiles when the interesting part
%     of the map is a small fraction of it - an evoked CSD occupies a few
%     tens of milliseconds and saturates at the default.
%   * The mask outlines the cluster rather than blanking what is outside it:
%     a non-significant trend stays visible and honest, instead of being
%     hidden by the test.
%   * NaN bins are drawn transparent, so a wavelet's edge cone reads as "no
%     data" instead of taking the colour map's extreme.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 0).

    p = inputParser;
    p.addParameter('mask',      []);
    p.addParameter('maskStyle', 'outline');
    p.addParameter('title',     '');
    p.addParameter('xlabel',    'Time (s)');
    p.addParameter('ylabel',    'Frequency (Hz)');
    p.addParameter('event',     0);
    p.addParameter('clim',      []);
    p.addParameter('logf',      false);
    p.addParameter('colorbar',  false);
    p.addParameter('ydir',      'normal');   % 'reverse' for depth downwards
    p.addParameter('climPercentile', [2 98]);
    p.parse(varargin{:});
    a = p.Results;

    if isempty(ax), ax = gca; end
    M = squeeze(M);
    assert(ismatrix(M) && isequal(size(M), [numel(f) numel(t)]), 'plotTFRpanel:size', ...
        'M must be [numel(f) x numel(t)] = [%d x %d]; got [%s].', ...
        numel(f), numel(t), strjoin(string(size(M)), ' x '));

    cl = a.clim;
    if isempty(cl), cl = st.zlim; end
    if isempty(cl), cl = localRobustLimits(M, st.diverging, a.climPercentile); end

    h.image = imagesc(ax, t, f, M);
    if strcmpi(st.interp, 'bilinear')
        set(h.image, 'Interpolation', 'bilinear');   % R2022a+; ignored if absent
    end
    set(ax, 'YDir', a.ydir, 'FontSize', st.fontSize, 'Layer', 'top', ...
            'TickDir', 'out', 'Box', 'off');
    clim(ax, cl);
    colormap(ax, st.colormap);
    if a.logf, set(ax, 'YScale', 'log'); end
    xlim(ax, [t(1) t(end)]);
    ylim(ax, [f(1) f(end)]);

    % NaN bins (a wavelet's edge cone, a dropped channel) must read as absent.
    % Left opaque they take the colormap's end colour, which on a diverging map
    % looks exactly like a strong effect.
    alpha = double(isfinite(M));
    set(ax, 'Color', 'w');

    h.maskLine = gobjects(0);
    if ~isempty(a.mask) && ~strcmpi(a.maskStyle, 'none')
        mask = squeeze(logical(a.mask));
        assert(isequal(size(mask), size(M)), 'plotTFRpanel:maskSize', ...
            'mask must match M ([%s] vs [%s]).', ...
            strjoin(string(size(mask)), ' x '), strjoin(string(size(M)), ' x '));
        if any(strcmpi(a.maskStyle, {'fade', 'both'}))
            alpha = alpha .* (0.35 + 0.65 * double(mask));
        end
        if any(strcmpi(a.maskStyle, {'outline', 'both'})) && any(mask(:))
            hold(ax, 'on');
            [~, h.maskLine] = contour(ax, t, f, double(mask), [0.5 0.5], ...
                                      'LineColor', 'k', 'LineWidth', 1.2);
            hold(ax, 'off');
        end
    end
    set(h.image, 'AlphaData', alpha);

    if ~isempty(a.event)
        hold(ax, 'on');
        for e = a.event(:)'
            xline(ax, e, 'k-', 'LineWidth', 1, 'Alpha', 0.6);
        end
        hold(ax, 'off');
    end

    % Interpreter 'none': channel and condition labels carry underscores,
    % which TeX would silently turn into subscripts.
    if ~isempty(a.title)
        title(ax, a.title, 'FontWeight', 'normal', 'Interpreter', 'none');
    end
    if ~isempty(a.xlabel), xlabel(ax, a.xlabel); end
    if ~isempty(a.ylabel), ylabel(ax, a.ylabel); end
    h.colorbar = gobjects(0);
    if a.colorbar, h.colorbar = colorbar(ax); end
end

function cl = localRobustLimits(M, diverging, pct)
    v = M(isfinite(M));
    if isempty(v), cl = [0 1]; return; end
    lo = prctile(v, pct(1)); hi = prctile(v, pct(2));
    if diverging
        m = max(abs([lo hi]));
        if m == 0, m = 1; end
        cl = [-m m];
    else
        if hi <= lo, hi = lo + eps(lo) + 1; end
        cl = [lo hi];
    end
end
