function plot_fireRateColMap(fireRate, conditions, param, opt)
% firerate comes as {cluster, phase}. We want to plot time vs cluster 
% in each panel (all clusters piled up), per phase. With a line to mark 
% the baseline end.

% Is it normalized or not?
norm = true;
if iscell(fireRate{1,1}), norm = false; end

%% Default params
param.ylabel        = {'cluster'}; % rate label
param.ytick         = 0:10:1000; % fire rate ticks
param.yticklabels   = {mat2cell(param.ytick,1)}; % cluster labels
param.xlabel        = {'time (s)'}; % time label
% % whole x-axis
% param.xtick         = (0:1000:diff([param.interval(1) param.interval(2)]))/param.stepSz; % time ticks
% param.xticklabels   = {mat2cell((param.xtick-40)*param.stepSz/1000,1)}; % time labels
% minimal x-axis
param.xtick          = (0:2000:diff([param.interval(1) param.interval(2)]))/param.stepSz; % time ticks
param.xticklabels{1} = {'-2', 'Ini','StimOn', '4', '6', '8', '10'};

%% Figure size
if strcmpi('adaptive', param.size)
    param.screen.size   = get(0, 'ScreenSize');  
    param.screen.width  = param.screen.size(3);
    param.screen.height = param.screen.size(4);
else
    param.screen.width  = param.size(1); 
    param.screen.height = param.size(2);
end

%% Blocks description
blocks.Pos = {'A' 'B' 'A' 'B' 'A' 'C' 'A' 'C'};
blocks.Amb = {'A' 'B' 'A' 'B' 'C' 'B' 'A' 'C'};
blocks.Neg = {'A' 'B' 'A' 'B' 'C' 'B' 'C' 'A'};

% Current application
blks{1} = blocks.(conditions.type);
blks{2} = {'NS' 'FS'};

%% Create plotting matrices from data
if norm
    for p = 1:size(fireRate,2)
        for i = 1:size(fireRate{1,p},1)
            for c = 1:size(fireRate,1)
                if all(isnan(fireRate{c,p}(i,:)))
                    try plotmatrix{p,i}(c,:) = zeros(size(plotmatrix{p,i}(c-1,:)));
                    catch plotmatrix{p,i}(c,:) = zeros(1,param.xtick(end)); end
                else
                    plotmatrix{p,i}(c,:) = fireRate{c,p}(i,:);
                end
            end
        end
    end
else
    for p = 1:size(fireRate,2)
        for i = 1:size(fireRate{1,p},2)
            for c = 1:size(fireRate,1)
                meanfr = mean(fireRate{c,p}{i},1,'includenan');
                if all(isnan(meanfr))
                    try plotmatrix{p,i}(c,:) = zeros(size(plotmatrix{p,i}(c-1,:)));
                    catch plotmatrix{p,i}(c,:) = zeros(1,param.xtick(end)); end
                else
                    plotmatrix{p,i}(c,:) = meanfr;
                end
            end
        end
    end
end    
% minmax{p,i} = [min(min(plotmatrix{p,i})) max(max(plotmatrix{p,i}))];

% OPTIONAL. Sort matrices as lower ordinal for higher lower column
for p = 1:size(fireRate,2)
    [plotmatrix{p,1}, sortidx{p}] = sortrows(plotmatrix{p,1});
    for i = 2:size(fireRate{1,p},1)
        plotmatrix{p,i} = sortrows(plotmatrix{p,i}(sortidx{p},:));
    end
end

%% Plots
for p = 1:size(fireRate,2)
    nblocks{p} = ~cellfun(@isempty, plotmatrix(p,:));
    fig = figure('Visible', param.visible);
    tlo{p} = tiledlayout(fig, sum(nblocks{p})/2, 2, "TileSpacing","tight");    
        if p == 1, tlo{p}.Parent.Position = [0, 0, round(param.screen.height), round(param.screen.width)]; % fig vertical
        else,      tlo{p}.Parent.Position = [0, 0, round(param.screen.width), round(param.screen.height)]; % fig landscape
        end

    for i = 1:sum(nblocks{p})
        sp{i} = nexttile;
            imagesc(sp{i}, plotmatrix{p,i});
            
            sp{i}.XLim = [param.xtick(1)-.5 param.xtick(end)+.5];
            sp{i}.YLim = [.5 size(plotmatrix{p,i},1)+.5];
                sp{i}.YDir = "normal";
            if ~norm, sp{i}.CLim = [0 max(max(plotmatrix{p,i}))*0.75];
            else,     sp{i}.CLim = [-5 5]; end
            sp{i}.Box = 'off';
            sp{i}.TickDir = "out";
            
            if p == 1 && i > 6
                sp{i}.XAxis.Label.String = param.xlabel;
                    sp{i}.XAxis.TickValues = param.xtick;
                    sp{i}.XAxis.TickLabels = param.xticklabels{1};
                    sp{i}.XAxis.TickLabelRotation = 0;
            elseif p == 1 && i <= 6
                sp{i}.XAxis.Visible = 'off';
            end

            if p == 2
                sp{i}.XAxis.Label.String = param.xlabel;
                    sp{i}.XAxis.TickValues = param.xtick;
                    sp{i}.XAxis.TickLabels = param.xticklabels{1};
                    sp{i}.XAxis.TickLabelRotation = 0;
            end

            if ismember(i,[1 3 5 7])
                sp{i}.YAxis.Label.String = param.ylabel;
                sp{i}.YAxis.TickValues = param.ytick;
            else
                sp{i}.YAxis.Visible = 'off';
            end
            
            xline(40, LineWidth=2, Color='k', LineStyle='--')
            xline(80, LineWidth=1, Color='k', LineStyle='--')

            sp{i}.Subtitle.String = blks{p}{i}; % TODO
    end
    
    cb = colorbar;
    if norm
        strtitl = 'NormFR';
        colormap('cool');
    else     
        strtitl = 'FR'; 
        colormap('parula');
    end

    cb.Label.String = strtitl;
    if p == 1, cb.Position = [0.925 0.35 0.017 0.30];
    else, cb.Position = [0.94 0.35 0.017 0.30]; end

    if p == 1,    tlo{p}.Title.String = ['Cluster ', strtitl, ' Block'];
    else,         tlo{p}.Title.String = ['Cluster ', strtitl, ' Stim']; end

    % Save figure per alignment&cluster    
    exportgraphics(tlo{p}, fullfile(param.save,'fireRate', ...
                        [opt.SavFileName,'_', tlo{p}.Title.String, '.png']), ...
                        'Resolution', param.Resolution);
    close all hidden
    close all force
end
end