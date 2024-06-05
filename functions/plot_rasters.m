function plot_rasters(neurons, events, spike, opt, param)
% Will take neuron-trial data and plot a series of basic rasters,
% histograms and other statistics to inspect clusters in relation to task
% events. A variable number of options can be given to modify plots without
% need to go to low level functions.

%% Default options.
if ~isfield(param,'res'),           param.res            = false;        end
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'treatment'),     param.treatment      = false;        end
if ~isfield(param,'genstats'),      param.genstats       = true;         end
if ~isfield(param,'pooledstats'),   param.pooledstats    = true;         end
if ~isfield(param,'plotcol'),       param.plotcol        = zeros(1,3);   end
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'square';     end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 3;            end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
if ~isfield(param,'timelim'),       param.timelim        = [-500 2000];  end
if ~isfield(param,'timelimItiOn'),  param.timelimItiOn   = [0 2500];    end
if ~isfield(param,'binSize'),       param.binSize        = 100;  end
if ~isfield(param,'stepSz'),        param.stepSz         = 10;  end
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;  end
if ~isfield(param,'interval'),      param.interval       = [0 (param.timelim(2)-param.timelim(1))/10];  end

%% Default figure attributes. 
% Raster plots
% if ~isfield(param,'raster')
    param.raster.ylabel = {'Trial # (single cluster)'}; % trial label
    param.raster.xlabel = {'msec'};   % time label
    param.raster.ytick = [0:50:1000]; % trial ticks
    param.raster.xtick = [param.timelim(1):250:param.timelim(2)]; % time ticks
    param.raster.yticklabels = {mat2cell(param.raster.ytick,1)}; % ticks label
    param.raster.xticklabels = {mat2cell(param.raster.xtick,1)}; % ticks label
% end
% PSH
% if ~isfield(param,'psh')
    param.psh.ylabel = {'spikes/s'}; % trial label
    param.psh.xlabel = {'msec'};   % time label
    param.psh.ytick = [0:5:60]; % trial ticks
    param.psh.xtick = [0:param.stepSz*5:param.interval(2)/param.stepSz]; % time ticks
    param.psh.yticklabels = {mat2cell(param.psh.ytick,1)}; % ticks label
    param.psh.xticklabels = {mat2cell(param.psh.xtick*param.stepSz,1)}; % ticks label
% end
% Pooled raster
% if ~isfield(param,'poolraster')
    param.poolraster.ylabel = {'Trial # (all clusters)'}; % trial label
    param.poolraster.xlabel = {'msec'};   % time label
    param.poolraster.ytick = [0:200:10000]; % trial ticks
    param.poolraster.xtick = [param.timelim(1):250:param.timelim(2)]; % time ticks
    param.poolraster.yticklabels = {mat2cell(param.poolraster.ytick,1)}; % ticks label
    param.poolraster.xticklabels = {mat2cell(param.poolraster.xtick,1)}; % ticks label
% end

%% Load needs
% if param.res
%     load(fullfile(opt.behavFiles, "res.mat"), 'par')
% end

%% Get important options and particular details
nalign = numel(opt.alignto);

% screen adaptive
screen.size = get(0, 'ScreenSize');   
screen.width = screen.size(3);
screen.height = screen.size(4);

%% Initialize levels
if param.treatment == true
    param.levels = 1; 
    param.trial_change = 1;
    % check existence of defined treatments
    if isfield(opt,'trEvents')
        % how many
        ntreatments = numel(opt.trEvents);
        % have they start/end or just a single application time?
        % add levels accordingly, e.g. basal/treatment_present/post (+2) or basal/post (+1)
        for i=1:ntreatments
            param.levels = param.levels + numel(events.(opt.trEvents{i}).time{1});
            param.trial_change = [param.trial_change events.(opt.trEvents{i}).trial{1}];
        end
    end
else
    param.levels = 1; % force to integer
    param.trial_change = 1;
end

% add last trial to the level-limits vector
param.trial_change = [param.trial_change length(neurons.(opt.alignto{1}){1})];

%% Individual stats, cluster by cluster
% TODO Add ISI, Autocorr, trace, others

% For each alignment
for a = 1:nalign
    % for each cluster
    for c = 1:length(neurons.(opt.alignto{a}))
         if param.genstats
            % For each cluster
            param.raster.title = ['Aligned to: ', opt.alignto{a}];
            param.raster.subtitle = ['cluster: ', spike.label{c}];
            
            % Special case for itiON, since there is no -time
            if a == 1, tlim = param.timelimItiOn;
            else,      tlim = param.timelim;
            end

            % Initialize figure
            fig = figure('visible', param.visible); % switch visibility
            set(fig, 'Position', [0, 0, round(screen.width), round(screen.height)]); % Set fig size as screen
            
            % Event-aligned Spike Raster
            subplot(2,1,1)
            trialCounter = 1;
            for lvl=1:param.levels
                % raster plot
                trialCounter = plotRaster(neurons.(opt.alignto{a}){c}(param.trial_change(lvl):param.trial_change(lvl+1)), ... % spikes
                                        trialCounter,               ... % trialCounter
                                       'plotcol',    param.plotcol(lvl,:), ... % color per align (for now)
                                       'spkwidth',   param.spkWidth,       ...
                                       'linelength', param.lineLength,     ...
                                       'plotstyle',  param.plotStyle,      ...
                                       'timelim',    tlim);
            end
            % Prettify
            if a == 1, prettify(param.raster, 1) % Special case for itiOn only
            else,      prettify(param.raster),  end

            % Event-aligned PSH
            subplot(2,1,2)
            for lvl=1:param.levels
                upperY2(lvl) = plotPSTH(neurons.(opt.alignto{a}){c}(param.trial_change(lvl):param.trial_change(lvl+1)), ... % spikes
                                param.stepSz,   ... % stepSz
                                param.binSize,  ...  % binSize
                                param.interval,  ... % interval
                                param.smpRate,  ...  % samples per second in the feeded data
                                'plotcol',      param.plotcol(lvl,:),...
                                'meanline',     '-',...
                                'smoothplot',   true);
            end
            % Prettify
            prettify(param.psh)
            
         end

         % Save figure per alignment&cluster
         exportgraphics(fig,fullfile(opt.trialSorted,'genplots',['raster_',opt.alignto{a},'_',spike.label{c},'.png']),'Resolution',300);
         close all

    end

    %% For ALL clusters
    if param.pooledstats
        % Trick the params to change color every cluster
        param.levels = length(neurons.(opt.alignto{a})); % number clusters
        trialperclus = length(neurons.(opt.alignto{a}){1}); % trials per cluster
        param.trial_change = (1:trialperclus:(param.levels*trialperclus)+1);
        param.plotcol = [0 0 0; 0.6350 0.0780 0.1840]; % two colors to alternate

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
        for lvl=1:param.levels
            trialCounter = plotRaster(poolneurons(param.trial_change(lvl):param.trial_change(lvl+1)-1), ... % spikes
                                      trialCounter,               ... % trialCounter
                                      'plotcol',    param.plotcol(mod(lvl,2)+1,:), ... % alternate color per cluster
                                      'spkwidth',   2,       ...
                                      'linelength', 1,     ...
                                      'plotstyle',  param.plotStyle,      ...
                                      'timelim',    param.timelim);
        end
        % Prettify
        prettify(param.poolraster)

    end
    
    % Save figure 
    exportgraphics(fig,fullfile(opt.trialSorted,'genplots',['raster_',opt.alignto{a},'_pooled.png']),'Resolution',300);
    close all

end %function