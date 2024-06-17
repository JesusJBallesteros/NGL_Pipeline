function plot_pooledstats(neurons, opt, param)
    %% TODO

if ~isfield(param,'pooledstats'),   param.pooledstats    = true;         end
if ~isfield(param,'res'),           param.res            = false;        end
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = [600 1000];   end
if ~isfield(param,'treatment'),     param.treatment      = false;        end
if ~isfield(param,'plotcol'),       param.plotcol        = zeros(1,3);   end
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'square';     end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 2;            end
if ~isfield(param,'lineLength'),    param.lineLength     = .5;           end
if ~isfield(param,'timelim'),       param.timelim        = [-500 2000];  end
if ~isfield(param,'timelimItiOn'),  param.timelimItiOn   = [0 2500];     end
if ~isfield(param,'timelimRwd'),    param.timelimRwd     = [-2000 500];  end

% Pooled raster
param.poolraster.ylabel = {'Cluster # (all trials)'}; % trial label
param.poolraster.xlabel = {'ms'};   % time label
param.poolraster.ytick = 0:length(neurons.itiOn{1,1}):length(neurons.itiOn{1,1})*length(neurons.(opt.alignto{1})); % trial ticks
param.poolraster.xtick = param.timelim(1):250:param.timelim(2); % time ticks
param.poolraster.yticklabels = {mat2cell((param.poolraster.ytick/length(neurons.itiOn{1,1}))+1,1)}; % ticks label
param.poolraster.xticklabels = {mat2cell(param.poolraster.xtick,1)}; % ticks label

%% Get details
nalign = numel(opt.alignto);

% Figure size
if strcmpi('adaptive', param.size)
    screen.size = get(0, 'ScreenSize');  
    screen.width =  screen.size(3);
    screen.height = screen.size(4);
else
    screen.width =  param.size(1); 
    screen.height = param.size(2);
end

%% For ALL clusters
% For each alignment
for a = 1:nalign
    if param.pooledstats
        % Special params to change color every cluster, only for this specific plot
        param.pooled_levels = length(neurons.(opt.alignto{a})); % number clusters
        trialperclus = length(neurons.(opt.alignto{a}){1}); % trials per cluster
        param.pooled_trial_change = (1:trialperclus:(param.pooled_levels*trialperclus)+1);
        param.pooled_plotcol = [0 0 0; 0.6350 0.0780 0.1840]; % two colors to alternate
    
        % Allocate and concatenate all cells from neurons cell array (pile
        % up all trials along all clusters)
    %         poolneurons = cell(param.levels*trialperclus,1);
        poolneurons = cat(1, neurons.(opt.alignto{a}){:});
    
        % Fig title and subtitle
        param.poolraster.title = ['Aligned to: ', opt.alignto{a}];
        param.poolraster.subtitle = ['ALL clusters'];
    
        % Initialize figure
        fig = figure('visible', param.visible); % switch visibility
        set(fig, 'Position', [0, 0, round(screen.width), round(screen.height)]); % Set fig size as screen
            
        % Event-aligned Spike Raster
        trialCounter = 1;
        for lvl = 1:param.pooled_levels
            trialCounter = plotRaster(poolneurons(param.pooled_trial_change(lvl):param.pooled_trial_change(lvl+1)-1), ... % spikes
                                      trialCounter,               ... % trialCounter
                                      'plotcol',    param.pooled_plotcol(mod(lvl,2)+1,:), ... % alternate color per cluster
                                      'spkwidth',   2,       ...
                                      'linelength', 1,     ...
                                      'plotstyle',  param.plotStyle,      ...
                                      'timelim',    param.timelim);
        end
        % Prettify
        prettify(param.poolraster)
            ylim([0 length(neurons.(opt.alignto{a}))+1])
            clear poolneurons

        %TODO Plot a 'driftmap' basically a time vs depth plot colored for
        % spike amplitude, for all clusters.
        plot_driftmap(sp.st, sp.amps, sp.depths);
    
        %TODO A probe deph amp pdf
        [pdfs, cdfs, ampBins, depthBins] = computeAndPlotWFampHist(ksDir, varargin);
    end
end

% Save figure 
exportgraphics(fig,fullfile(opt.trialSorted,'genplots',['raster_',opt.alignto{a},'_pooled.png']),'Resolution',300);
close all
end
