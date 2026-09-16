function [fig, figFile] = plotCSD(csd, opt, varargin)
% plotCSD  Sink-source map down the shank, with the LFP it came from.
%
% PURPOSE:
%   Depth against time, sinks blue and sources red, with the laminar LFP
%   traces beside it. The two panels answer different questions - where
%   current moved, and what the raw signal looked like - and disagreeing with
%   each other is informative, so both are drawn.
%
% USAGE:
%   [fig, f] = plotCSD(csd, opt)
%   [fig, f] = plotCSD(csd, opt, 'align', 'stim', 'traces', false)
%
% INPUTS:
%   csd - from computeCSD.
%   opt - options; uses opt.lfp.plot.* through lfpStyle.
%   Name/value pairs:
%     'align'  alignment name, for the title and file name
%     'traces' draw the LFP traces panel (default true)
%     'save'   true (default)
%     'title'  override
%
% OUTPUT:
%   fig, figFile
%
% NOTES:
%   * Depth runs downwards, as the probe does. The colour scale is symmetric
%     about zero, so white is "no net current" and the sign is readable; with
%     an asymmetric scale a weak source and a strong sink can take the same
%     colour.
%   * When the CSD was computed with Vaknin's extension the top and bottom
%     depths rest on an assumption about the potential beyond the probe, not
%     on a measurement. They are marked on the axis.
%
% Last modified 16.09.2026 (Jesus) - new (LFP analysis Phase 2, CSD).

    p = inputParser;
    p.addParameter('align',  '');
    p.addParameter('traces', true);
    p.addParameter('save',   true);
    p.addParameter('title',  '');
    p.parse(varargin{:});
    a = p.Results;

    st = lfpStyle(opt, 'diverging');
    nPanels = 1 + logical(a.traces);
    fig = figure('Visible', st.visible, 'Color', 'w', ...
                 'Position', [80 80 st.figSize(1) st.figSize(2)]);
    tl = tiledlayout(fig, 1, nPanels, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(tl);
    % Wide percentiles: the response lasts tens of milliseconds out of a whole
    % epoch, so the default robust limits would put the baseline in view and
    % saturate the sink that is the point of the figure.
    h = plotTFRpanel(ax, csd.time, csd.depth, csd.csd, st, ...
        'climPercentile', [0.2 99.8], ...
        'ylabel', 'Depth (\mum)', 'ydir', 'reverse', 'colorbar', true, ...
        'title', sprintf('CSD  (%d trials, %d smoothing pass(es))', ...
                         csd.info.nTrials, csd.info.smoothPasses));
    h.colorbar.Label.String = sprintf('%s   (blue = sink)', csd.info.units);
    if csd.info.edgeEstimated
        localMarkEdges(ax, csd.depth);
    end

    if a.traces
        ax2 = nexttile(tl);
        localTraces(ax2, csd, st);
    end

    if isempty(a.title)
        a.title = localTitle(csd, a.align);
    end
    title(tl, a.title, 'FontWeight', 'normal', 'Interpreter', 'none');

    figFile = '';
    if a.save
        tags = {sprintf('shank%g', csd.info.shank)};
        figFile = saveLFPfigure(fig, 'CSD', opt, 'area', csd.info.area, ...
                                'align', a.align, 'tags', tags);
    end
end

% ---------------- helpers ----------------
function localTraces(ax, csd, st)
% LFP per contact, stacked at its own depth so the panel shares the y axis
% with the CSD map. Scaled to the largest deflection, printed in the label.
    span = max(abs(csd.lfp(:)));
    if span == 0 || ~isfinite(span), span = 1; end
    gap = median(diff(csd.depth));
    if isempty(gap) || gap == 0, gap = 1; end
    scale = 0.45 * gap / span;
    hold(ax, 'on');
    for k = 1:size(csd.lfp, 1)
        plot(ax, csd.time, csd.depth(k) - scale * csd.lfp(k, :), ...
             'Color', [0.25 0.28 0.35], 'LineWidth', 0.7);
    end
    xline(ax, 0, 'k-', 'Alpha', 0.6);
    hold(ax, 'off');
    set(ax, 'YDir', 'reverse', 'FontSize', st.fontSize, 'TickDir', 'out', 'Box', 'off');
    xlim(ax, [csd.time(1) csd.time(end)]);
    ylim(ax, [min(csd.depth) - gap, max(csd.depth) + gap]);
    xlabel(ax, 'Time (s)');
    title(ax, sprintf('LFP (peak %.0f \\muV, up = negative)', span), ...
          'FontWeight', 'normal');
end

function localMarkEdges(ax, depth)
% The Vaknin depths are extrapolated, not measured; say so on the axis.
    yl = [depth(1) depth(end)];
    txt = {'(edge: assumed flat)', '(edge: assumed flat)'};
    for k = 1:2
        text(ax, min(get(ax, 'XLim')), yl(k), txt{k}, 'FontSize', 7, ...
             'Color', [0.35 0.35 0.4], 'VerticalAlignment', 'middle', ...
             'HorizontalAlignment', 'left', 'Margin', 1);
    end
end

function s = localTitle(csd, align)
    bits = {};
    if ~isempty(csd.info.area), bits{end+1} = csd.info.area; end
    bits{end+1} = sprintf('shank %g', csd.info.shank);
    bits{end+1} = sprintf('%g um spacing', csd.info.spacing_um);
    if ~isempty(align), bits{end+1} = ['aligned to ' align]; end
    s = strjoin(bits, '  |  ');
end
