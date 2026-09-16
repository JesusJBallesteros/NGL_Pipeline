function [fig, figFile] = plotTFRcontrast(res, opt, varargin)
% plotTFRcontrast  The standard event-centered power figure: A, B, and A - B.
%
% PURPOSE:
%   Three panels per band - each condition against baseline, then their
%   difference with the significant cluster outlined. The two condition panels
%   share one colour scale so they can be compared by eye; the difference gets
%   its own, symmetric about zero.
%
% USAGE:
%   [fig, f] = plotTFRcontrast(res, opt)                 % one band
%   [fig, f] = plotTFRcontrast({resLow, resHigh}, opt)   % one row per band
%   plotTFRcontrast(res, opt, 'save', false)             % keep it open
%
% INPUTS:
%   res - result struct from computeTFRcontrast, or a cell of them (one per
%         band). Bands are drawn highest-frequency first, matching the
%         quick-look TFR plot.
%   opt - options; uses opt.lfp.plot.* through lfpStyle, and opt.analysis /
%         opt.SavFileName for the file name.
%   Name/value pairs:
%     'save'    true (default) writes the PNG via saveLFPfigure
%     'visible' override for the figure's visibility
%     'title'   override for the figure title
%
% OUTPUT:
%   fig     - figure handle (closed by the caller; still open on return)
%   figFile - path written, '' when 'save' is false
%
% NOTES:
%   Bands keep separate colour scales. Power falls steeply with frequency, so
%   a scale shared across a 4 Hz and a 60 Hz band would flatten the faster one
%   to a single colour - the opposite of the quick-look plot's choice, which
%   shares one scale to compare bands' magnitudes directly. Here the question
%   is the shape of each band's response, not which band is larger.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 1).

    p = inputParser;
    p.addParameter('save',    true);
    p.addParameter('visible', '');
    p.addParameter('title',   '');
    p.parse(varargin{:});
    a = p.Results;

    if ~iscell(res), res = {res}; end
    res = res(~cellfun(@isempty, res));
    assert(~isempty(res), 'plotTFRcontrast:empty', 'nothing to plot.');

    stSeq = lfpStyle(opt);
    stDiv = lfpStyle(opt, 'diverging');
    visible = a.visible; if isempty(visible), visible = stSeq.visible; end

    % Highest-frequency band on top, as in plotTrialparsedTFR_example.
    fLow = cellfun(@(r) min(r.freq), res);
    [~, order] = sort(fLow, 'descend');
    res = res(order);
    nB = numel(res);

    figSize = stSeq.figSize;
    fig = figure('Visible', visible, 'Color', 'w', ...
                 'Position', [80 80 figSize(1) max(figSize(2), 260 * nB)]);
    tl = tiledlayout(fig, nB, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for b = 1:nB
        r = res{b};
        % One scale for the two conditions, so their panels are comparable.
        condLim = localLimits([r.meanA(:); r.meanB(:)], false);
        bottom  = (b == nB);
        xlab = ''; if bottom, xlab = 'Time (s)'; end

        plotTFRpanel(nexttile(tl), r.time, r.freq, r.meanA, stSeq, ...
            'clim', condLim, 'title', localTitle(r, r.spec.labelA, r.nA), ...
            'xlabel', xlab);
        % The pair shares one scale, so one colorbar serves both; it sits on
        % the second panel, next to the difference's own.
        h2 = plotTFRpanel(nexttile(tl), r.time, r.freq, r.meanB, stSeq, ...
            'clim', condLim, 'title', localTitle(r, r.spec.labelB, r.nB), ...
            'xlabel', xlab, 'ylabel', '', 'colorbar', true);
        h2.colorbar.Label.String = r.norm.units;

        h3 = plotTFRpanel(nexttile(tl), r.time, r.freq, r.diff, stDiv, ...
            'mask', r.mask, 'title', localDiffTitle(r), 'xlabel', xlab, ...
            'ylabel', '', 'colorbar', true);
        h3.colorbar.Label.String = ['\Delta ' r.norm.units];
    end

    if isempty(a.title)
        a.title = localFigTitle(res{1});
    end
    title(tl, a.title, 'FontWeight', 'normal');

    figFile = '';
    if a.save
        figFile = saveLFPfigure(fig, 'TFRcontrast', opt, 'area', res{1}.area, ...
                                'align', res{1}.align, 'tags', {res{1}.spec.label});
    end
end

% ---------------- helpers ----------------
function s = localTitle(r, label, n)
    s = sprintf('%s (n = %d)', label, n);
    if numel(r.freq) > 0
        % ft_freqanalysis returns the frequencies it could resolve, not the
        % ones asked for, so round rather than printing 12.9395-44.9219 Hz.
        s = sprintf('%s  |  %g-%g Hz', s, round(min(r.freq)), round(max(r.freq)));
    end
end

function s = localDiffTitle(r)
    s = sprintf('%s - %s', r.spec.labelA, r.spec.labelB);
    if isempty(r.stat)
        s = [s '  (not tested)'];
        return
    end
    p = [r.stats.pPos(:); r.stats.pNeg(:)];
    p = p(isfinite(p));
    if isempty(p) || min(p) >= r.stats.alpha
        s = sprintf('%s  (no cluster, min p = %.3f)', s, min([p; 1]));
    elseif min(p) <= r.stats.pFloor
        s = sprintf('%s  (cluster p < %.3f)', s, r.stats.pFloor);
    else
        s = sprintf('%s  (cluster p = %.3f)', s, min(p));
    end
end

function s = localFigTitle(r)
    bits = {};
    if ~isempty(r.area),  bits{end+1} = r.area;  end
    if ~isempty(r.align), bits{end+1} = ['aligned to ' r.align]; end
    bits{end+1} = r.spec.request;
    s = strjoin(bits, '  |  ');
end

function cl = localLimits(v, diverging)
    v = v(isfinite(v));
    if isempty(v), cl = [0 1]; return; end
    lo = prctile(v, 2); hi = prctile(v, 98);
    if diverging
        m = max(abs([lo hi])); if m == 0, m = 1; end
        cl = [-m m];
    else
        if hi <= lo, hi = lo + 1; end
        cl = [lo hi];
    end
end
