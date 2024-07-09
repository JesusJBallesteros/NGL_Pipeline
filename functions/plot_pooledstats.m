function plot_pooledstats(neurons, opt, param)
    %% TODO

if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'treatment'),     param.treatment      = true;         end
if ~isfield(param,'pooled_plotcol'),param.pooled_plotcol = [0.0 0.0 0.0; 0.6 0.1 0.2]; end
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'lines';      end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 1;            end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
if ~isfield(param,'baseline'),      param.baseline       = 500;          end
if ~isfield(param,'post'),          param.post           = 2500;         end
if ~isfield(param,'timelim'),       param.timelim        = [-param.baseline param.post]; end
if ~isfield(param,'timelimItiOn'),  param.timelimItiOn   = [   0 8000];  end

% Pooled raster
param.poolraster.ylabel = {'Cluster # (all trials)'}; % cluster/trial label
param.poolraster.ytick = 0:length(neurons.itiOn{1,1}):length(neurons.itiOn{1,1})*length(neurons.(opt.alignto{1})); % trial ticks
param.poolraster.yticklabels = {mat2cell((param.poolraster.ytick/length(neurons.itiOn{1,1}))+1,1)}; % ticks label

param.poolraster.xlabel = {'ms'}; % time label
param.poolraster.xtick = param.timelim(1):param.baseline:param.timelim(2); % time ticks
param.poolraster.xtickiti = param.timelimItiOn(1):1000:param.timelimItiOn(2); % time ticks iti
param.poolraster.xticklabels = {mat2cell(param.poolraster.xtick,1)}; % ticks label
param.poolraster.xticklabelsiti = {mat2cell(param.poolraster.xtickiti,1)}; % ticks label iti

%% For ALL clusters
% For each alignment
for a = 1:numel(opt.alignto)
    
    % Cases
    if strcmpi(opt.alignto{a}, 'itiOn'),     timelim = param.timelimItiOn; param.size = [800 1000]; % itiOn
    elseif strcmpi(opt.alignto{a}, 'rwd'),   timelim = param.timelim;      param.size = [300 1000]; % needed?
    else,                                    timelim = param.timelim;      param.size = [300 1000]; % all others
    end

    % To change color every cluster, only for this specific kind of figure
    param.pooled_levels = length(neurons.(opt.alignto{a})); % number clusters
    trialperclus = length(neurons.(opt.alignto{a}){1}); % trials per cluster
    param.pooled_trial_change = (1:trialperclus:(param.pooled_levels*trialperclus)+1);

    % Allocate and concatenate all cells from neurons cell array (pile
    % up all trials along all clusters)
    poolneurons = cat(1, neurons.(opt.alignto{a}){:});
    % poolneurons = cellfun(@(x) x*1000, poolneurons, 'UniformOutput', false); % In case time units need to be changed

    % prepare title and subtitle
    param.poolraster.title = ['Aligned to: ', opt.alignto{a}];
    param.poolraster.subtitle = ['ALL clusters'];

    % Initialize
    fig = figure('visible', param.visible); % switch visibility
    set(fig, 'Position', [0, 0, round(param.size(1)), round(param.size(2))]); % Set fig size as screen
        
    % Event-aligned Spike Raster for all clusters in session
    trialCounter = 1;
    for lvl = 1:param.pooled_levels
        trialCounter = plotRaster(poolneurons(param.pooled_trial_change(lvl):param.pooled_trial_change(lvl+1)-1), ... % spikes
                                  trialCounter,               ... % trialCounter
                                  'plotcol',    param.pooled_plotcol(mod(lvl,2)+1,:), ... % alternate color per cluster
                                  'spkwidth',   1,       ...
                                  'linelength', 1,     ...
                                  'plotstyle',  param.plotStyle,      ...
                                  'timelim',    timelim);
    end
    clear poolneurons % Free some memory

    % Prettify
    prettify(param.poolraster, a)
        ylim([0 param.poolraster.ytick(end)])

    % Save figure 
    exportgraphics(fig,fullfile(opt.analysis,'genplots',['raster_',opt.alignto{a},'_pooled.png']),'Resolution',600);
    close all
end

end
