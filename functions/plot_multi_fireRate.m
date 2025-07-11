function plot_multi_fireRate(fireRate, param, opt)
% firerate comes Normalized, as {cluster}. We want to plot time vs cluster 
% (all clusters piled up), per phase. With some delineation of trial structure

% Is it normalized
% norm = true;

%% Default params
if ~isfield(param,'xtick')
    param.xtick = (0:(0-param.interval(1)):diff([param.interval(1) param.interval(2)]))/param.stepSz; % time ticks
end
if ~isfield(param,'xticklabels'),   param.xticklabels = {((param.xtick*param.stepSz)+param.interval(1))/1000}; end
if ~isfield(param,'timelabel'),     param.timelabel     = 'time (s)';                           end
if ~isfield(param,'blockchange'),   param.blockchange   = [];                                   end
if ~isfield(param,'ylabel'),        param.ylabel        = {'cluster #'};                          end % rate label
if ~isfield(param,'ytick'),         param.ytick         = 0:10:1000;                            end % fire rate ticks
if ~isfield(param,'yticklabels'),   param.yticklabels   = {mat2cell(param.ytick,1)};            end % cluster labels
if ~isfield(param,'inibin'),        param.inibin = param.xtick(find(param.xticklabels{1}==0));  end
if ~isfield(param,'blockchange'),   param.blockchange    = [];                                  end

%% Figure size
if strcmpi('adaptive', param.size)
    param.screen.size   = get(0, 'ScreenSize');  
    param.screen.width  = param.screen.size(3);
    param.screen.height = param.screen.size(4);
else
    param.screen.width  = param.size(1); 
    param.screen.height = param.size(2);
end

%% Create plotting matrices from data
plotmatrix = cell2mat(fireRate);

% Remove invalid elements
plotmatrix(all(plotmatrix==Inf,2),:) = [];
plotmatrix(all(isnan(plotmatrix),2),:) = [];

% Sort matrices by max firing rate, at shorter time bin, increasingly
for i = 1:size(plotmatrix,1)
    maxfr(i,1) = max(plotmatrix(i,:));
    maxtb(i,1) = find(plotmatrix(i,:)==maxfr(i,1),1);
end

[~, sorted] = sort(maxtb);
sortplotmatrix = plotmatrix(sorted,:);

%% Plots
figure('Visible', param.visible, 'Position', [0 0 700 400]);
fig = imagesc(sortplotmatrix);
    fig.Parent.Title.String = 'allAvClusters';
    % X
    fig.Parent.XLim = [param.xtick(1) size(sortplotmatrix,2)];
    fig.Parent.XAxis.Label.String = param.timelabel;
    fig.Parent.XAxis.TickValues = param.xtick;
    fig.Parent.XAxis.TickLabels = param.xticklabels{1};
    fig.Parent.XAxis.TickLabelRotation = 0;
    % Y
    fig.Parent.YLim = [.5 size(sortplotmatrix,1)+.5];
    fig.Parent.YAxis.Label.String = param.ylabel;
    fig.Parent.YDir = "normal";
    % Box and scale
    fig.Parent.CLim = [-4 4]; 
    fig.Parent.Box = 'off';
    fig.Parent.TickDir = "out";
    % Zero Aligment line
    xline(param.inibin, LineWidth=2, Color='k', LineStyle='--')
    % colorbar
    cb = colorbar;
    cb.Position = [0.91 0.35 0.017 0.30];
    cb.Label.String = 'Z';
    cb.Label.Position = [1.1 0 0];
    cb.Label.HorizontalAlignment = "center";
    cb.Label.VerticalAlignment = "bottom";
    cb.TickDirection = 'out';
    % Color map
    cRedBlue = redwhiteblue(fig.Parent.CLim(1), fig.Parent.CLim(2), size(colormap(fig.Parent),1));
    colormap(fig.Parent, cRedBlue)
    
    % Save figure per alignment&cluster    
    if ~exist(fullfile(opt.analysis,'plots','multi_fr'),"dir")
        mkdir(fullfile(opt.analysis,'plots','multi_fr'))
    end
    exportgraphics(fig.Parent, fullfile(opt.analysis,'plots','multi_fr', ...
                        [opt.SavFileName,'_', fig.Parent.Title.String, '.png']), ...
                        'Resolution', param.Resolution);
    close all hidden
    close all force
end