function plot_BarErrAllData(data_all, nCats, categories)
%
%

%% Set Visualization Parameters
alpha = 0.4;
grey = [97,97,97]./255;
okabe_ito = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66]]./255;
vmin = round(min(data_all, [], 'all'), 1)-0.1;
vmax = round(max(data_all, [], 'all'), 1)+0.1;
if mod(vmin,0.2) ~= 0, vmin = vmin - 0.1; end
if mod(vmax,0.2) ~= 0, vmax = vmax + 0.1; end

%% Data
x_axis = 1:nCats;
data_mean = mean(data_all, 1, 'omitnan');
data_std  = std(data_all, 1, 'omitnan');
nSelTrial = size(data_all, 1);

%% Plot 
figure(); hold on;

% Bars
b = bar(x_axis, data_mean);
b.FaceColor = 'flat';
% color
for n = 1:nCats
    b.CData(n,:) = okabe_ito(n,:);
end

% Data points
s = cell(size(x_axis));
for n = 1:nCats
    rnd = unifrnd(-0.2,0.2,1,nSelTrial);
    s{n} = scatter(ones(1,nSelTrial)*n+rnd, data_all(:,n), 80, grey, 'filled');
    s{n}.MarkerFaceAlpha = 0.65;
    s{n}.LineWidth = 1.5;
end

% Error bars
er = cell(size(x_axis));
for n = 1:nCats
    er{n} = errorbar(x_axis(n), data_mean(n), data_std(n), data_std(n));
    er{n}.Color = okabe_ito(n,:);
    if n == nCats
        er{n}.Color = [209,206,23]./255; % darker yellow for visibility
    end
end

% Make pretty
xticks(x_axis);
xticklabels(categories);
yticks(linspace(vmin,vmax,5));
xlim([0.2,nCats+0.8]);
ylim([vmin,vmax]);
xlabel('Categories');
ylabel({'Data magnitude'; '(Units)'});
title('Data Title');

% Set properties
set(b,'EdgeColor','none','BaseValue',vmin+0.004,'ShowBaseLine','off', ...
    'FaceAlpha',alpha);
set([er{:}],'CapSize',20,'LineStyle','none','LineWidth',5);
set(gca,'FontSize',20,'LineWidth',3,'TickDir','out','TickLength',[0.015, 0.025], ...
    'XTickLabelRotation',0,'Box','off','FontWeight','bold',...
    'Position',[0.2216,0.1545,0.7362,0.7784]);
set(gcf,'Color','w','Position',[453,327,582,509]);
end