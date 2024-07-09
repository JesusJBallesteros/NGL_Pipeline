function plot_trialclusters(neurons, events, spike, opt, param)
% Will take trial-long data and plot a series of basic rasters,
% to inspect long dynamics in relation to task events.

%% Default options.
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = [1000 600];   end
if ~isfield(param,'treatment'),     param.treatment      = false;        end
if ~isfield(param,'plotcol'),       param.plotcol        = [ 0  0  0;
                                                            .6 .6 .6;
                                                            .3 .3 .3]; 
end
% Rasters
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'lines';      end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 1;            end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
                                    param.timelim        = [0 8000];    % Hard coded, TODO
if ~isfield(param,'plotevent'),     param.plotevent      = [1 3 7];      end
% PSH
if ~isfield(param,'binSize'),       param.binSize        = 200;          end
if ~isfield(param,'stepSz'),        param.stepSz         = 10;           end
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end

%% Default figure attributes. 
% Use ItiOn alignment only
a = 1;

% Raster plot
param.raster.ylabel = {'Trial #'}; % trial label
param.raster.ytick = 0:50:1000; % trial ticks
param.raster.yticklabels = {mat2cell(param.raster.ytick,1)}; % ticks label
param.raster.xlabel = {'s'};   % time label
param.raster.xtick = param.timelim(1):1000:param.timelim(2); % time ticks
param.raster.xticklabels = {mat2cell(param.raster.xtick/1000,1)}; % ticks label

% PSH
param.psh.ylabel = {'spikes/s'}; % rate label
param.psh.ytick = 0:5:60; % fire rate ticks
param.psh.yticklabels = {mat2cell(param.psh.ytick,1)}; % rate labels
param.psh.xlabel = {'s'}; % time label
param.psh.xtick = (param.timelim(1):1000:param.timelim(2))/param.stepSz; % time ticks
param.psh.xticklabels = {mat2cell(param.psh.xtick*param.stepSz/1000,1)}; % time labels

% Driftmap
param.driftmap.ylabel = {'tempAmpl'}; % ampl label
param.driftmap.xlabel = {'min'}; % time label
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
% for each cluster
for c = 1:length(neurons.(opt.alignto{a}))
    % Titles
    param.raster.title = 'Whole trial';
    param.raster.subtitle = ['cluster: ', spike.label{c}];
    
    % Initialize
    fig = figure('visible', param.visible); % switch visibility
    set(fig, 'Position', [0, 0, round(screen.width), round(screen.height)]); % Set fig size as screen
    
    % Trial-long Spike Raster
    subplot(3,1,1)
        trialCounter = 1;
        for lvl = 1:param.levels
            % prepare range of trials to plot
            trialrange = param.trial_change(lvl):param.trial_change(lvl+1);
                if lvl > 1, trialrange(1) = []; end % remove repeated trial at beginning, when treatments exist
                if length(trialrange) < 2, continue, end

            % raster plot
            trialCounter = plotRaster(neurons.(opt.alignto{a}){c}(trialrange), ... % spikes
                                    trialCounter,               ... % trialCounter
                                   'plotcol',    param.plotcol(lvl,:), ... % color per align (for now)
                                   'spkwidth',   param.spkWidth,       ...
                                   'linelength', param.lineLength,     ...
                                   'plotstyle',  param.plotStyle);
        end

            % Plot requested events (param.plotevent)
            if ~isempty(param.plotevent)
                for trial = 1:size(events.(opt.alignto{a}).code,1)
                    for ev = param.plotevent
                        % Check event presence in trial
                        evidx = find(events.(opt.alignto{a}).code{trial,1} == ev);
                        % If any, plot
                        if ~isempty(evidx)
                            color = [];
                            if ev == 1, color = 'b'; end % stimOn1
                            if ev == 3, color = 'r'; end % bhv
                            if ev == 7, color = 'g'; end % rwd
                            line(1000*[events.(opt.alignto{a}).time{trial,1}(evidx) events.(opt.alignto{a}).time{trial,1}(evidx)], ...
                                 [trial trial+1], 'Color', color, 'LineWidth', 2)
                        end
                    end
                end
            end

        % Prettify
        prettify(param.raster);
            xlim(param.timelim)

    % Trial-long PSH
    subplot(3,2,[3,4])
        for lvl=1:param.levels
            % prepare range of trials to plot
            trialrange = param.trial_change(lvl):param.trial_change(lvl+1);
                if lvl > 1, trialrange(1) = []; end % remove repeated trial at beginning, when treatments exist
                if length(trialrange) < 2, continue, end
            % Plot
            upperY = plotPSTH(neurons.(opt.alignto{a}){c}(param.trial_change(lvl):param.trial_change(lvl+1)), ... % spikes
                            param.stepSz,   ... % stepSz
                            param.binSize,  ...  % binSize
                            param.timelim,  ... % interval
                            param.smpRate,  ...  % samples per second in the feeded data
                            'plotcol',      param.plotcol(lvl,:),...
                            'meanline',     '-',...
                            'smoothplot',   true);
        end
        % Prettify
        prettify(param.psh)
        if upperY <= 5, ylim([0 5]); end

    % Session-long 'driftmap'
    subplot(3,2,[5,6])
        scatter(spike.timestamp{c}/60, spike.templampl{c}, 5, "black", "filled");
        % Prettify
        prettify(param.driftmap);
            ylim([min(spike.templampl{c})*0.9 max(spike.templampl{c})*1.1]);
            xlim([0 spike.timestamp{c}(end)/60+.2]);

        % treatment window 
        if param.treatment
            for i = 1:ntreatments
                if isfield(events,(opt.trEvents{i})) % If there is a na3 treatment field
                    xy = [events.(opt.trEvents{i}).time{1}(1)/60, 0];
                    w = (events.(opt.trEvents{i}).time{1}(2)-events.(opt.trEvents{i}).time{1}(1))/60; 
                    h = 1000; % hard coded
                    
                    rectangle('Position', [xy(1), xy(2), w, h], ...             % [x, y, w, h]
                              'FaceColor',[.6 .6 .6 .35], 'LineStyle', 'none'); % [r g b alpha], no line.
                end
            end
        end

    % Save figure per alignment&cluster
    exportgraphics(fig,fullfile(opt.analysis,'genplots',['fulltrial_',spike.label{c},'.png']),'Resolution',600);
    close all
end

end