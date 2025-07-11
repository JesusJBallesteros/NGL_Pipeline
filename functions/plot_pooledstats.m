function plot_pooledstats(neurons, opt, param, input)
    %% TODO

if ~isfield(param,'visible'),       param.visible        = 'off';        end
% if ~isfield(param,'treatment'),     param.treatment      = true;         end
if ~isfield(param,'pooled_plotcol'),param.pooled_plotcol = [0.0 0.0 0.0; 0.6 0.1 0.2]; end
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'lines';      end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 1;            end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
if ~isfield(param,'baseline'),      param.baseline       = 500;          end
if ~isfield(param,'post'),          param.post           = 2500;         end
if ~isfield(param,'timelim'),       param.timelim        = [-param.baseline param.post]; end
if ~isfield(param,'timelimItiOn'),  param.timelimItiOn   = [   0 8000];  end
if ~isfield(param,'size'),          param.size           = [1900 1000];   end

% overrule para.treatment, to deprecate (TODO)
% param.treatment = opt.treatment;

%% Set
toalignto = opt.alignto;

% Add case for social interactions
if opt.plotSocial
    toalignto = [toalignto , 'interactions'];
end

% How many alignments
nalign = numel(toalignto); 

%% For ALL clusters
% For each alignment
for a = 1:nalign
    % Pooled raster
    param.poolraster.ylabel = {'Cluster #'}; % cluster/trial label
    param.poolraster.ytick = 0:length(neurons.(toalignto{a}){1,1}):length(neurons.(toalignto{a}){1,1})*length(neurons.(toalignto{a})); % trial ticks
    param.poolraster.yticklabels = {mat2cell((param.poolraster.ytick/length(neurons.(toalignto{a}){1,1}))+1,1)}; % ticks label
    
    param.poolraster.xlabel = {'time (ms)'}; % time label
    param.poolraster.xtick = param.timelim(1):param.baseline:param.timelim(2); % time ticks
    param.poolraster.xticklabels = {mat2cell(param.poolraster.xtick,1)}; % ticks label

    % Cases
    timelim = param.timelim; % General, and rwd
    if strcmpi(toalignto{a}, 'itiOn') % itiOn
        timelim = param.timelimItiOn;
        param.size = [800 1000]; 
        param.poolraster.xtick = timelim(1):1000:timelim(2); % time ticks iti
        param.poolraster.xticklabel = {mat2cell(param.poolraster.xtick,1)}; % ticks label iti
    elseif strcmpi(toalignto{a}, 'interactions') % social interactions
        longst_interact = 10000;
        timelim = [-2000 longst_interact];
        param.poolraster.xlabel = {'time (s)'}; % time label
        param.poolraster.xtick = timelim(1):2000:timelim(2); % time ticks
        param.poolraster.xticklabels = {mat2cell(param.poolraster.xtick/1000,1)}; % ticks label
    end

    % To change color every cluster, only for this specific kind of figure
    param.pooled_levels = length(neurons.(toalignto{a})); % number clusters
    trialperclus = length(neurons.(toalignto{a}){1}); % trials per cluster
    param.pooled_trial_change = (1:trialperclus:(param.pooled_levels*trialperclus)+1);

    % Allocate and concatenate all cells from neurons cell array (pile
    % up all trials along all clusters)
    poolneurons = cat(1, neurons.(toalignto{a}){:});
    % poolneurons = cellfun(@(x) x*1000, poolneurons, 'UniformOutput', false); % In case time units need to be changed

    % prepare title and subtitle
    param.poolraster.title = ['Aligned to: ', toalignto{a}];
    param.poolraster.subtitle = ['ALL clusters'];

    %% Figure
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

    %% Save figure 
    if ~exist(fullfile(input.analysis,'genplots'),'dir')
        mkdir(fullfile(input.analysis,'genplots'));
    end
    exportgraphics(fig,fullfile(input.analysis,'genplots',['raster_',toalignto{a},'_pooled.png']),'Resolution',600);
    close all
end

end
