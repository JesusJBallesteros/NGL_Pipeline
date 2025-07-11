function plot_single_fireRate(fireRate, normFireRate, param, opt)
% Funtion to create two plots, side by side, showing the raw firing rate
% and the normalized firing rate for a single cluster, across all trials.
% INPUTS: 'fireRate': [trials, nbins] array as calculated from 'calcFireRate.m'
%         'normFireRate': [trials, nbins] array as calculated from 'calcFireRate.m'
%         'param': set of parameters given for the processing pipeline
%         'opt': options collected for current session
%
% OUTPUTS: saved plots in the session folder
%
% Jesus. 14.03.1985

%% Defaults
if ~isfield(param,'xtick')
    param.xtick = (0:(0-param.interval(1)):diff([param.interval(1) param.interval(2)]))/param.stepSz; % time ticks
end
if ~isfield(param,'xticklabels'),   param.xticklabels = {((param.xtick*param.stepSz)+param.interval(1))/1000}; end
if ~isfield(param,'timelabel'),     param.timelabel = 'time (s)';                               end
if ~isfield(param,'inibin'),        param.inibin = param.xtick(find(param.xticklabels{1}==0));  end
if ~isfield(param,'blockchange'),   param.blockchange    = [];                                  end

%% Plot simple heatmaps for each calculated rate.
hp = figure('Visible', param.visible, 'Position', [1100 200 700 400]);
tlo = tiledlayout(hp, 1, 2, "TileSpacing","compact");

if ~isempty(fireRate)
    fr = nexttile;
      imagesc(fireRate);
        % X axis 
        fr.XLim = [1 size(fireRate,2)+.5];
        fr.XAxis.Label.String = param.timelabel;
        fr.XAxis.TickValues = param.xtick;
        fr.XAxis.TickLabels = param.xticklabels;
        fr.XAxis.TickLabelRotation = 0;
        % Y axis 
        fr.YLim = [0 size(fireRate,1)+.5];
        fr.YAxis.Label.String = 'trial';
        fr.YDir = "normal";
        % Plot
        fr.Box = 'off';
        fr.TickDir = "out";
        % Ini
        xline(param.inibin, LineWidth=2, Color='w', LineStyle='--')
        % Color scale
        fr.CLim = [0 max(max(fireRate))*0.9];
        % Color Bar
        cb1 = colorbar;
        cb1.Position = [0.49 0.35 0.017 0.30];
        cb1.Label.String = 'sp/s';
        cb1.Label.Position = [1 (max(max(fireRate)))/2 0];
        cb1.Label.HorizontalAlignment = "center";
        cb1.Label.VerticalAlignment = "bottom";
        cb1.TickDirection = 'out';
        % Color map
        colormap(fr, "hot")
        % Blocks
        if ~isempty(param.blockchange)
            for b=1:size(param.blockchange,2)
                yline(param.blockchange(1,b),LineWidth=1, Color='w', LineStyle='-')
            end
        end
end

if ~isempty(normFireRate)
    Nfr = nexttile;
      imagesc(normFireRate);
        % X axis
        Nfr.XLim = [1 size(normFireRate,2)+.5];
        Nfr.XAxis.Label.String = param.timelabel;
        Nfr.XAxis.TickValues = param.xtick;
        Nfr.XAxis.TickLabels = param.xticklabels;
        Nfr.XAxis.TickLabelRotation = 0;
        Nfr.YAxis.Visible = 'off';
        % Y axis        
        Nfr.YLim = [0 size(normFireRate,1)+.5];
        Nfr.YDir = "normal";
        % Plot
        Nfr.Box = 'off';
        Nfr.TickDir = "out";
        % Ini
        xline(param.inibin, LineWidth=2, Color='k', LineStyle='--')
        % Color Scale
        Nfr.CLim = [-4 4];
        % Color bar
        cb2 = colorbar;
        cb2.Position = [0.91 0.35 0.017 0.30];
        cb2.Label.String = 'Z';
        cb2.Label.Position = [1.1 0 0];
        cb2.Label.HorizontalAlignment = "center";
        cb2.Label.VerticalAlignment = "bottom";
        cb2.TickDirection = 'out';
        % Color map
        cRedBlue = redwhiteblue(Nfr.CLim(1), Nfr.CLim(2), size(colormap(Nfr),1));
        colormap(Nfr, cRedBlue)
        % Blocks
        if ~isempty(param.blockchange)
            for b=1:size(param.blockchange,2)
                yline(param.blockchange(1,b),LineWidth=1, Color='k', LineStyle='-')
            end
        end
end
tlo.Title.String = [opt.alignto{param.cl(1)}, ' c', num2str(param.cl(2))];

% Save figure per alignment&cluster    
if ~exist(fullfile(opt.analysis,'plots','single_fr'),"dir")
    mkdir(fullfile(opt.analysis,'plots','single_fr'))
end
exportgraphics(tlo, fullfile(opt.analysis,'plots','single_fr', ...
                    [opt.SavFileName,'_', tlo.Title.String, '.png']), ...
                    'Resolution', param.Resolution);

close all hidden
close all force
end