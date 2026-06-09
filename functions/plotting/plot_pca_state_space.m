function plot_pca_state_space(result, plotCfg)
% plot_pca_state_space  2D + 3D neural-state-space figure for one PCA result.
%
% PURPOSE:
%   Render a side-by-side 2D (PC1 vs PC2) and 3D (PC1 vs PC2 vs PC3)
%   trajectory figure from a calculate_pca_from_pool output. One trace
%   per condition (coloured), with the chosen uncertainty variant
%   overlaid:
%     - 'singleTrials' : thin per-session marginal traces in the
%                        condition's colour at low alpha. Answers
%                        "does the mean generalize across sessions?"
%     - 'ciTube'       : trial-bootstrap 5-95% envelope rendered as a
%                        shaded ribbon (2D) and per-axis CI crosshairs
%                        at subsampled timepoints (3D).
%   Filled circle at t=0 (alignment) on every mean trajectory; open
%   diamond at the end of the window. Variance-explained percentages
%   appear in each axis label.
%
% USAGE:
%   plot_pca_state_space(result, plotCfg);
%
% INPUTS:
%   result   - struct from calculate_pca_from_pool. Required fields:
%                .traj_mean, .explained, .timeAxis, .conditions,
%                .traj_session (singleTrials fallback), .ciLo/.ciHi (ciTube),
%                .sessionKeys, .nClustPerSess.
%              Optional:
%                .traj_trial   {Ncond x 1}, [Nbins x K x Ntrials_c]. When
%                              populated (single-session case), takes
%                              precedence over traj_session for the
%                              'singleTrials' overlay.
%   plotCfg  - struct:
%                .variant      'singleTrials' | 'ciTube'    (default 'singleTrials')
%                .titlePrefix  char (default '')
%                .outFile      full PNG path (required)
%                .palette      [Ncond x 3] (default lines(Ncond))
%                .sessionAlpha 0..1 line alpha for per-session traces (default 0.18)
%                .ciAlpha      0..1 fill alpha for CI band (default 0.20)
%                .ciStride     integer, draw 3D crosshair every N bins (default 10)
%                .visible      logical, figure 'Visible' (default false)
%
% Last modified 09.06.2026 (Jesus)

    if nargin < 2, plotCfg = struct(); end
    if ~isfield(plotCfg,'variant'),      plotCfg.variant      = 'singleTrials'; end
    if ~isfield(plotCfg,'titlePrefix'),  plotCfg.titlePrefix  = '';              end
    if ~isfield(plotCfg,'palette'),      plotCfg.palette      = [];              end
    if ~isfield(plotCfg,'sessionAlpha'), plotCfg.sessionAlpha = 0.18;            end
    if ~isfield(plotCfg,'ciAlpha'),      plotCfg.ciAlpha      = 0.20;            end
    if ~isfield(plotCfg,'ciStride'),     plotCfg.ciStride     = 10;              end
    if ~isfield(plotCfg,'visible'),      plotCfg.visible      = true;           end
    assert(isfield(plotCfg,'outFile') && ~isempty(plotCfg.outFile), ...
        'NGL:plot_pca_state_space:noOut', 'plotCfg.outFile is required.');

    traj_mean    = result.traj_mean;     % [Nbins x K x Ncond]
    explained    = result.explained;     % [K x 1]
    conditions   = result.conditions;
    [Nbins, K, Ncond] = size(traj_mean);

    palette = plotCfg.palette;
    if isempty(palette), palette = lines(max(2, Ncond)); end

    has3D = K >= 3;

    fig = figure('Visible', ternaryVis(plotCfg.visible), 'Position', [80 80 1300 600]);
    tl  = tiledlayout(fig, 1, ternaryInt(has3D, 2, 1), 'Padding','compact', 'TileSpacing','compact');

    %% 2D panel: PC1 vs PC2.
    ax2 = nexttile(tl);
    hold(ax2,'on');
    handles2D = gobjects(Ncond, 1);
    for c = 1:Ncond
        col = palette(c, :);
        xMean = traj_mean(:, 1, c);
        yMean = traj_mean(:, 2, c);

        switch lower(plotCfg.variant)
            case 'singletrials'
                traces = localPickTraces(result, c);   % [Nbins x K x N]
                N      = size(traces, 3);
                for s = 1:N
                    xs = traces(:, 1, s);
                    ys = traces(:, 2, s);
                    if all(isnan(xs)) || all(xs == 0 & ys == 0), continue; end
                    plot(ax2, xs, ys, '-', 'Color', [col, plotCfg.sessionAlpha], ...
                        'LineWidth', 0.8);
                end
            case 'citube'
                localFillCITube2D(ax2, result, c, palette(c,:), plotCfg.ciAlpha);
            otherwise
                error('NGL:plot_pca_state_space:badVariant', ...
                    'Unknown variant ''%s'' (expected ''singleTrials'' or ''ciTube'').', ...
                    plotCfg.variant);
        end

        handles2D(c) = plot(ax2, xMean, yMean, '-', 'Color', col, 'LineWidth', 2.0);
        % t=0 marker (alignment event)
        [~, iZero] = min(abs(result.timeAxis));
        plot(ax2, xMean(iZero), yMean(iZero), 'o', ...
             'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        % end marker
        plot(ax2, xMean(end), yMean(end), 'd', ...
             'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
    end
    xlabel(ax2, sprintf('PC1 (%.1f%%)', explained(1)));
    ylabel(ax2, sprintf('PC2 (%.1f%%)', explained(min(2, K))));
    grid(ax2,'on'); axis(ax2,'square'); box(ax2,'off');
    title(ax2, '2D state space  |  filled o = t=0, diamond = end');

    if Ncond >= 1
        legend(ax2, handles2D, conditions, 'Location', 'bestoutside');
        legend(ax2, 'boxoff');
    end

    %% 3D panel: PC1 vs PC2 vs PC3.
    if has3D
        ax3 = nexttile(tl);
        hold(ax3,'on');
        view(ax3, 3);
        for c = 1:Ncond
            col = palette(c, :);
            xMean = traj_mean(:, 1, c);
            yMean = traj_mean(:, 2, c);
            zMean = traj_mean(:, 3, c);

            switch lower(plotCfg.variant)
                case 'singletrials'
                    traces = localPickTraces(result, c);
                    N      = size(traces, 3);
                    for s = 1:N
                        xs = traces(:, 1, s);
                        ys = traces(:, 2, s);
                        zs = traces(:, 3, s);
                        if all(isnan(xs)) || all(xs == 0 & ys == 0 & zs == 0), continue; end
                        plot3(ax3, xs, ys, zs, '-', ...
                              'Color', [col, plotCfg.sessionAlpha], 'LineWidth', 0.8);
                    end
                case 'citube'
                    localCICrosshair3D(ax3, result, c, col, plotCfg.ciStride);
            end

            plot3(ax3, xMean, yMean, zMean, '-', 'Color', col, 'LineWidth', 2.0);
            [~, iZero] = min(abs(result.timeAxis));
            plot3(ax3, xMean(iZero), yMean(iZero), zMean(iZero), 'o', ...
                  'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
            plot3(ax3, xMean(end),   yMean(end),   zMean(end),   'd', ...
                  'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
        end
        xlabel(ax3, sprintf('PC1 (%.1f%%)', explained(1)));
        ylabel(ax3, sprintf('PC2 (%.1f%%)', explained(2)));
        zlabel(ax3, sprintf('PC3 (%.1f%%)', explained(3)));
        grid(ax3,'on'); box(ax3,'off');
        title(ax3, '3D state space');
    end

    %% Suptitle with variant label and pool composition.
    sessSummary = sprintf('%d sessions (%d clust each on avg)', ...
                          numel(result.sessionKeys), round(mean(result.nClustPerSess)));
    if isempty(plotCfg.titlePrefix)
        sgtitle(fig, sprintf('%s  |  %s', upper(plotCfg.variant), sessSummary));
    else
        sgtitle(fig, sprintf('%s  |  %s  |  %s', ...
            plotCfg.titlePrefix, upper(plotCfg.variant), sessSummary));
    end

    %% Save.
    outDir = fileparts(plotCfg.outFile);
    if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
    exportgraphics(fig, plotCfg.outFile, 'Resolution', 300);
    close(fig);
end

% ------------------------------------------------------------------------
function traces = localPickTraces(result, c)
% Prefer real single-trial projections (traj_trial) when present and
% non-empty for this condition. Otherwise fall back to per-session
% marginal trajectories (traj_session). Both share the same
% [Nbins x K x N] shape so the caller can iterate uniformly.
    if isfield(result, 'traj_trial') && numel(result.traj_trial) >= c ...
            && ~isempty(result.traj_trial{c}) && size(result.traj_trial{c}, 3) > 0
        traces = result.traj_trial{c};
        return
    end
    traces = result.traj_session{c};
end

% ------------------------------------------------------------------------
function localFillCITube2D(ax, result, c, col, fillAlpha)
% Draw a filled 2-D ribbon between the per-axis CI bounds in PC1/PC2.
% Approximation: at each timepoint, the corners are
%   (ciLo_PC1(t), ciLo_PC2(t)) -> (ciHi_PC1(t), ciLo_PC2(t))
%   (ciHi_PC1(t), ciHi_PC2(t)) -> (ciLo_PC1(t), ciHi_PC2(t))
% We render two patches: the "x-uncertainty" band uses the mean PC2 path,
% and the "y-uncertainty" band uses the mean PC1 path. Together they
% give a sense of the trajectory's reliability.
    if ~isfield(result,'ciLo') || all(isnan(result.ciLo(:))), return; end
    xMean = result.traj_mean(:, 1, c);
    yMean = result.traj_mean(:, 2, c);
    xLo   = result.ciLo(:, 1, c);  xHi = result.ciHi(:, 1, c);
    yLo   = result.ciLo(:, 2, c);  yHi = result.ciHi(:, 2, c);

    bad = isnan(xLo) | isnan(xHi) | isnan(yLo) | isnan(yHi);
    if all(bad), return; end

    % X-band (PC1 uncertainty): polygon spanning xLo..xHi, mean PC2.
    px = [xLo(~bad); flipud(xHi(~bad))];
    py = [yMean(~bad); flipud(yMean(~bad))];
    patch(ax, px, py, col, 'FaceAlpha', fillAlpha, 'EdgeColor', 'none');

    % Y-band (PC2 uncertainty): polygon spanning yLo..yHi, mean PC1.
    px = [xMean(~bad); flipud(xMean(~bad))];
    py = [yLo(~bad);   flipud(yHi(~bad))];
    patch(ax, px, py, col, 'FaceAlpha', fillAlpha, 'EdgeColor', 'none');
end

function localCICrosshair3D(ax, result, c, col, stride)
% At every `stride`-th timepoint, draw thin per-axis CI segments centred
% on the mean trajectory point. Three short segments (one per PC axis).
    if ~isfield(result,'ciLo') || all(isnan(result.ciLo(:))), return; end
    xMean = result.traj_mean(:, 1, c);
    yMean = result.traj_mean(:, 2, c);
    zMean = result.traj_mean(:, 3, c);
    xLo   = result.ciLo(:, 1, c);  xHi = result.ciHi(:, 1, c);
    yLo   = result.ciLo(:, 2, c);  yHi = result.ciHi(:, 2, c);
    zLo   = result.ciLo(:, 3, c);  zHi = result.ciHi(:, 3, c);

    Nbins = numel(xMean);
    idx   = 1:stride:Nbins;
    for t = idx
        if any(isnan([xLo(t) xHi(t) yLo(t) yHi(t) zLo(t) zHi(t)])), continue; end
        plot3(ax, [xLo(t) xHi(t)], [yMean(t) yMean(t)], [zMean(t) zMean(t)], ...
              '-', 'Color', [col 0.6], 'LineWidth', 0.8);
        plot3(ax, [xMean(t) xMean(t)], [yLo(t) yHi(t)], [zMean(t) zMean(t)], ...
              '-', 'Color', [col 0.6], 'LineWidth', 0.8);
        plot3(ax, [xMean(t) xMean(t)], [yMean(t) yMean(t)], [zLo(t) zHi(t)], ...
              '-', 'Color', [col 0.6], 'LineWidth', 0.8);
    end
end

function s = ternaryVis(b)
    if b, s = 'on'; else, s = 'off'; end
end

function v = ternaryInt(b, a, c)
    if b, v = a; else, v = c; end
end
