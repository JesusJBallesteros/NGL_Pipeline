function [fig, figFile] = plotTaggingSpectrum(resp, opt, varargin)
% plotTaggingSpectrum  The tagging figure: where the response sits in the spectrum.
%
% PURPOSE:
%   One panel showing the z-scored spectrum with the tagging frequency and its
%   harmonics marked, and one showing the response per harmonic. Read together
%   they answer the two questions a tagging block raises: is there a peak where
%   the experiment put one, and how far up the harmonics does it go.
%
% USAGE:
%   [fig, f] = plotTaggingSpectrum(resp, opt)
%   [fig, f] = plotTaggingSpectrum(resp, opt, 'itpc', itpc, 'block', 'stream1.3')
%
% INPUTS:
%   resp - from taggingResponse.
%   opt  - options; uses opt.lfp.plot.* through lfpStyle, opt.analysis and
%          opt.SavFileName for the file name.
%   Name/value pairs:
%     'itpc'    from computeITPC; adds a phase-consistency panel
%     'block'   block name, for the title and file name
%     'area'    area name, for the title and file name
%     'channel' index or label to draw; default is the channel with the
%               largest summed response (see NOTES)
%     'save'    true (default)
%
% OUTPUT:
%   fig, figFile
%
% NOTES:
%   Defaulting to the strongest channel is a convenience for looking, not a
%   result: picking the best of many channels and then reading its z as a
%   p-value is circular. The panel labels which channel it drew, and the
%   per-channel numbers live in the saved .mat.
%
% Last modified 16.09.2026 (Jesus) - new (NFT project).

    p = inputParser;
    p.addParameter('itpc',    []);
    p.addParameter('block',   '');
    p.addParameter('area',    '');
    p.addParameter('channel', []);
    p.addParameter('save',    true);
    p.parse(varargin{:});
    a = p.Results;

    st = lfpStyle(opt);
    ch = localPickChannel(resp, a.channel);
    hasITPC = ~isempty(a.itpc);
    nPanels = 2 + hasITPC;

    fig = figure('Visible', st.visible, 'Color', 'w', ...
                 'Position', [80 80 st.figSize(1) 300 + 120 * nPanels]);
    tl = tiledlayout(fig, nPanels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- z-scored spectrum with the harmonics marked ----------------------
    ax = nexttile(tl);
    plot(ax, resp.freq, resp.zSpectrum(ch, :), 'Color', [0.30 0.35 0.45], 'LineWidth', 0.9);
    hold(ax, 'on');
    yl = [min(-2, min(resp.zSpectrum(ch, :))), ...
          max(4, 1.15 * max(resp.zSpectrum(ch, :)))];
    for k = 1:numel(resp.harmonics)
        col = [0.75 0.75 0.78];
        if resp.significant(ch, k), col = [0.84 0.38 0.30]; end
        xline(ax, resp.harmonics(k), '-', 'Color', col, 'LineWidth', 1.1, ...
              'Alpha', 0.9);
    end
    yline(ax, resp.info.zThreshold, ':', 'Color', [0.45 0.45 0.5]);
    plot(ax, resp.freq, resp.zSpectrum(ch, :), 'Color', [0.30 0.35 0.45], 'LineWidth', 0.9);
    hold(ax, 'off');
    ylim(ax, yl); xlim(ax, [0 max(resp.freq)]);
    set(ax, 'FontSize', st.fontSize, 'TickDir', 'out', 'Box', 'off');
    ylabel(ax, 'z vs neighbours');
    % Interpreter 'none': channel labels carry underscores, which TeX would
    % turn into subscripts (NCL_01 -> NCL with a subscript 0).
    title(ax, sprintf('%s   (harmonics of %g Hz marked; red = z > %.2f)', ...
          localChanName(resp, ch), resp.info.baseFreq, resp.info.zThreshold), ...
          'FontWeight', 'normal', 'Interpreter', 'none');

    % --- response per harmonic --------------------------------------------
    ax = nexttile(tl);
    b = bar(ax, 1:numel(resp.harmonics), resp.corrected(ch, :), 0.6);
    b.FaceColor = 'flat';
    for k = 1:numel(resp.harmonics)
        if resp.significant(ch, k)
            b.CData(k, :) = [0.84 0.38 0.30];
        else
            b.CData(k, :) = [0.80 0.80 0.83];
        end
    end
    set(ax, 'XTick', 1:numel(resp.harmonics), 'FontSize', st.fontSize, ...
            'TickDir', 'out', 'Box', 'off', 'XTickLabel', ...
            arrayfun(@(x) sprintf('%.2g', x), resp.harmonics, 'UniformOutput', false));
    xlabel(ax, 'Harmonic (Hz)');
    ylabel(ax, 'Amplitude - neighbours');
    title(ax, sprintf('summed over significant harmonics = %.3g (%d of %d)', ...
          resp.sumAmp(ch), resp.nSignificant(ch), numel(resp.harmonics)), ...
          'FontWeight', 'normal');

    % --- phase consistency -------------------------------------------------
    if hasITPC
        ax = nexttile(tl);
        it = a.itpc;
        b2 = bar(ax, 1:numel(it.freq), it.itpc(ch, :), 0.6);
        b2.FaceColor = 'flat';
        for k = 1:numel(it.freq)
            if it.significant(ch, k)
                b2.CData(k, :) = [0.26 0.52 0.78];
            else
                b2.CData(k, :) = [0.80 0.80 0.83];
            end
        end
        set(ax, 'XTick', 1:numel(it.freq), 'FontSize', st.fontSize, ...
                'TickDir', 'out', 'Box', 'off', 'XTickLabel', ...
                arrayfun(@(x) sprintf('%.2g', x), it.freq, 'UniformOutput', false));
        ylim(ax, [0 1]);
        xlabel(ax, 'Harmonic (Hz)'); ylabel(ax, 'ITPC');
        title(ax, sprintf('phase consistency over %d epochs (blue = Rayleigh p < %.3g)', ...
              it.nEpochs, it.alpha), 'FontWeight', 'normal');
    end

    title(tl, localTitle(resp, a), 'FontWeight', 'normal', 'Interpreter', 'none');

    figFile = '';
    if a.save
        tags = {sprintf('%gHz', resp.info.baseFreq)};
        if ~isempty(a.block), tags = [{a.block}, tags]; end
        figFile = saveLFPfigure(fig, 'NFT', opt, 'area', a.area, 'tags', tags);
    end
end

% ---------------- helpers ----------------
function ch = localPickChannel(resp, want)
    if isempty(want)
        [~, ch] = max(resp.sumAmp);
        return
    end
    if isnumeric(want)
        ch = want;
    else
        ch = find(strcmp(resp.label, want), 1);
        assert(~isempty(ch), 'plotTaggingSpectrum:channel', ...
            'channel ''%s'' is not in this spectrum.', want);
    end
end

function s = localChanName(resp, ch)
    if isfield(resp, 'label') && numel(resp.label) >= ch
        s = resp.label{ch};
    else
        s = sprintf('channel %d', ch);
    end
end

function s = localTitle(resp, a)
    bits = {};
    if ~isempty(a.block), bits{end+1} = a.block; end
    if ~isempty(a.area),  bits{end+1} = a.area;  end
    bits{end+1} = sprintf('tagged at %g Hz', resp.info.baseFreq);
    s = strjoin(bits, '  |  ');
end
