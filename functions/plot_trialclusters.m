function plot_trialclusters(neurons, events, spike, opt, param)
% Will take neuron-trial data and plot a series of basic rasters,
% histograms and other statistics to inspect clusters in relation to task
% events. A variable number of options can be given to modify plots without
% need to go to low level functions.

%% Default options.
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = [1000 600];   end
if ~isfield(param,'treatment'),     param.treatment      = false;        end
if ~isfield(param,'plotcol'),       param.plotcol        = [0.4 0.4 0.4; 0.6 0.1 0.2; 0.0 0.0 0.0]; end
% Rasters
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'square';     end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 1;            end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
if ~isfield(param,'timelim'),       param.timelim        = [0 8000];     end % Hard coded, TODO
% PSH
if ~isfield(param,'binSize'),       param.binSize        = 200;          end
if ~isfield(param,'stepSz'),        param.stepSz         = 20;           end
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end
if ~isfield(param,'interval'),      param.interval       = param.timelim;end

%% Default figure attributes. 
% Raster plot
param.raster.ylabel = {'Trial #'}; % trial label
param.raster.xlabel = {'ms'};   % time label
param.raster.ytick = 0:50:1000; % trial ticks
param.raster.xtick = param.timelim(1):1000:param.timelim(2); % time ticks
param.raster.yticklabels = {mat2cell(param.raster.ytick,1)}; % ticks label
param.raster.xticklabels = {mat2cell(param.raster.xtick,1)}; % ticks label

% PSH
param.psh.ylabel = {'spikes/s'}; % rate label
param.psh.xlabel = {'ms'}; % time label
param.psh.ytick = 0:5:60; % fire rate ticks
param.psh.xtick = (param.interval(1):1000:param.interval(2))/param.stepSz; % time ticks
param.psh.yticklabels = {mat2cell(param.psh.ytick,1)}; % rate labels
param.psh.xticklabels = {mat2cell(param.psh.xtick*param.stepSz,1)}; % time labels

% Driftmap
param.driftmap.ylabel = {'tempAmpl'}; % ampl label
param.driftmap.xlabel = {'time (min)'}; % time label
param.driftmap.ytick  = 'auto'; % ampl ticks
param.driftmap.xtick  = 'auto'; % time ticks
param.driftmap.yticklabels = {'auto'}; % ampl label
param.driftmap.xticklabels = {'auto'}; % time label

%% Get details
% Figure size
if strcmpi('adaptive', param.size)
    screen.size = get(0, 'ScreenSize');  
    screen.width =  screen.size(3);
    screen.height = screen.size(4);
else
    screen.width =  param.size(1); 
    screen.height = param.size(2);
end

%% Initialize levels
if param.treatment == true
    param.levels = 1; 
    param.trial_change = 1;
    % check existence of defined treatments
    if isfield(opt,'trEvents')
        % how many
        ntreatments = numel(opt.trEvents);
        % have they start/end or just a single application time?
        for i = 1:ntreatments
            % add levels accordingly, e.g. basal/treatment_present/post (+2) or basal/post (+1)
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

%% Whole trial raster and PSH, cluster by cluster
% Use ItiOn alignment
a = 1;
% for each cluster
for c = 1:length(neurons.(opt.alignto{a}))
    % Figure
    param.raster.title = 'Whole trial';
    param.raster.subtitle = ['cluster: ', spike.label{c}];
    
    % Initialize figure
    fig = figure('visible', param.visible); % switch visibility
    set(fig, 'Position', [0, 0, round(screen.width), round(screen.height)]); % Set fig size as screen
    
    % Event-aligned Spike Raster
    subplot(3,2,[1,2])
        trialCounter = 1;
        for lvl=1:param.levels
            % raster plot
            trialCounter = plotRaster(neurons.(opt.alignto{a}){c}(param.trial_change(lvl):param.trial_change(lvl+1)), ... % spikes
                                    trialCounter,               ... % trialCounter
                                   'plotcol',    param.plotcol(lvl,:), ... % color per align (for now)
                                   'spkwidth',   param.spkWidth,       ...
                                   'linelength', param.lineLength,     ...
                                   'plotstyle',  param.plotStyle);
        end
        % Prettify
        prettify(param.raster, a);
        xlim([0 8000])

    % Event-aligned PSH
    subplot(3,2,[3,4])
        for lvl=1:param.levels
            upperY = plotPSTH(neurons.(opt.alignto{a}){c}(param.trial_change(lvl):param.trial_change(lvl+1)), ... % spikes
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
        if upperY <= 5, ylim([0 5]); end

    % Plot a single 'driftmap'
    subplot(3,2,[5,6])
        scatter(spike.timestamp{c}/60, spike.templampl{c}, 5, "black", "filled");
        % Prettify
        prettify(param.driftmap);
            ylim([min(spike.templampl{c})*0.9 max(spike.templampl{c})*1.1]);
            xlim([0 spike.timestamp{c}(end)/60+.2]);

        % treatment window 
        if param.treatment == true
            for i = 1:ntreatments
                if isfield(events,(opt.trEvents{i})) % If there is a na3 treatment field
                    xy = [events.(opt.trEvents{i}).time{1}(1)/60, 0];
                    w = (events.(opt.trEvents{i}).time{1}(2)-events.(opt.trEvents{i}).time{1}(1))/60; 
                    h = (max(spike.templampl{c})*1.1-min(spike.templampl{c})*0.9);
                    
                    rectangle('Position', [xy(1), xy(2), w, h], ...             % [x, y, w, h]
                              'FaceColor',[.5 .5 .5 .35], 'LineStyle', 'none'); % [r g b alpha], no line.
                end
            end
        end

    % Save figure per alignment&cluster
    exportgraphics(fig,fullfile(opt.analysis,'genplots',['fulltrial_',spike.label{c},'.png']),'Resolution',600);
    close all
end

end