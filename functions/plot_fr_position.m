function plot_fr_position(tstamps, ts_t, pos, varargin)
%
%
%
%
%

%% Check if variable input are provided; if not, set a default value
if nargin < 4
    param.sr        = 60; % Sampling rate in Hz (1000 samples per second)
    param.sig       = 1; % Gaussian kernel width (500 ms)
    param.title     = 'Sequence ';  %
    param.nbins     = 50;        % Define grid resolution (higher = finer)
    param.subtitle  = '';           % Except if interesting for saving, empty
    param.shape     = [0 256; 256 256; 256 0; 530 0; 530 256; 768 256; 768 530; 530 530; 530 768; 256 768; 256 530; 0 530];
    param.visible   = 'on';         % see the figure
    param.AXvisible   = 'off';      % hide axes
    param.colmap    = 'hot';        % colormap
    param.flipY     = 0;        % Do not Reverse Y axis
    param.barunit   = 'sp/s';
    param.savepath  = pwd;
else
    param = varargin{1}; %
end

%% Call a function to calculate firing rate with a Gaussian smooth,
% outputting also a time vector that matches the whole range of timestamps
[fr_t, fr] = calculateFiringRate(tstamps, param.sr, param.sig);

%% Interpolate spatial positions sampling to the firing rate time vector.
if size(ts_t,1) > size(unique(ts_t),1)
    z = find(diff(ts_t)==0);
    ts_t(z) = [];
    pos(z,:) = [];
end

x_pos = interp1(ts_t, pos(:,1), fr_t, 'linear');
y_pos = interp1(ts_t, pos(:,2), fr_t, 'linear');

%% Define spatial binning
num_bins = param.nbins;
x_edges = linspace(min(x_pos), max(x_pos), num_bins+1); % X bins
y_edges = linspace(min(y_pos), max(y_pos), num_bins+1); % Y bins

% Assign each position to a spatial bin
[x_bin_idx, ~] = discretize(x_pos, x_edges);
[y_bin_idx, ~] = discretize(y_pos, y_edges);

%% Create a 2D matrix for accumulating firing rate
rate_matrix = zeros(num_bins, num_bins);
count_matrix = zeros(num_bins, num_bins);

for i = 1:length(fr)
    if ~isnan(x_bin_idx(i)) && ~isnan(y_bin_idx(i))
        rate_matrix(x_bin_idx(i), y_bin_idx(i)) = rate_matrix(x_bin_idx(i), y_bin_idx(i)) + fr(i);
        count_matrix(x_bin_idx(i), y_bin_idx(i)) = count_matrix(x_bin_idx(i), y_bin_idx(i)) + 1;
    end
end

% Compute average firing rate per spatial bin
avg_fr = rate_matrix ./ (count_matrix + (count_matrix == 0)); % Avoid division by zero

%% Plot the results and save
% Check/Create folder.
if ~exist(param.savepath,"dir"), mkdir(param.savepath), end

% % Continuous firing rate 
% fig = figure('visible', param.visible); % visibility; 'painters' allow for transparency
%     scatter(x_pos, y_pos, 10, fr', 'filled'); % Color-coded scatter plot
%     colorbar; % Add color bar to show firing rate intensity
%     colormap('hot'); % Use a colormap for better visualization
%     xlim([0 max(max(param.shape))]*1.1);
%     ylim([0 max(max(param.shape))]*1.1);
% 
%     ax = gca;  % get scatter properties
%     ax.Color = [0 0 0]; % black background
%     ax.XAxis.Visible = param.AXvisible; % show no show axis
%     ax.YAxis.Visible = param.AXvisible;  % show/ no show axis
%     if param.flipY % Set (0,0) on top left
%         set(ax,'XDir','normal')
%         set(ax,'YDir','reverse')
%     end
% 
%     % this creates the shape of the arena 
%     v = param.shape; % In this experiment, we removed one half
%     f = 1:size(v,1);
%     patch('Faces', f , 'Vertices' , v, ...
%           'EdgeColor', 'white', 'FaceColor', 'none', 'LineWidth', 1);
%     
%     % Get CLim
%     ax.CLim; % obtain Max value for colorbar
%     
%     % Color bar labels and limits
%     c = colorbar('Ticks', ax.CLim, 'TickLabels', {ceil(ax.CLim(1)), floor(ax.CLim(2))}, 'TickLength', 0, 'FontSize', 10);
%         c.Label.String = 'sp/s';
%         c.Label.Position = [1.25 ax.CLim(2)/2 0];    
%     title('Firing Rate vs Position');
%     
%     exportgraphics(fig, fullfile(param.savepath, [param.title, param.subtitle, '.png']), 'Resolution', 300);
%     close all
%     clear ax c fig

% Average firing rate
fig = figure('visible', param.visible); % visibility; 'painters' allow for transparency
    imagesc(x_edges, y_edges, avg_fr'); % Use imagesc for 2D plotting
    colorbar; % Add color bar
    colormap('hot'); % Use colormap for visual clarity
    xlim([0 max(max(param.shape))]*1.1);
    ylim([0 max(max(param.shape))]*1.1);

    ax = gca;  % get scatter properties
    ax.Color = [0 0 0]; % black background
    ax.XAxis.Visible = param.AXvisible; % show no show axis
    ax.YAxis.Visible = param.AXvisible;  % show/ no show axis

    if param.flipY % Set (0,0) on top left
        set(ax,'XDir','normal')
        set(ax,'YDir','reverse')
    end

    % this creates the shape of the arena 
    v = param.shape; % In this experiment, we removed one half
    f = 1:size(v,1);
    patch('Faces', f , 'Vertices' , v, ...
          'EdgeColor', 'white', 'FaceColor', 'none', 'LineWidth', 1);
    
    % Get CLim
    ax.CLim; % obtain Max value/2 for colorbar
    
    % Color bar labels and limits
    c = colorbar('Ticks', ax.CLim, 'TickLabels', {ceil(ax.CLim(1)), floor(ax.CLim(2))}, 'TickLength', 0, 'FontSize', 10);
        c.Label.String = 'Average sp/s';
        c.Label.Position = [1.25 ax.CLim(2)/2 0];    
    title('Average Firing Rate');

    exportgraphics(fig, fullfile(param.savepath, [param.title, param.subtitle, '_avg.png']), 'Resolution', 300);
    close all

end