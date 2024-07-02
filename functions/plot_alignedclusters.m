function plot_alignedclusters(neurons, events, spike, opt, param)
% Will take neuron-trial data and plot a series of basic rasters,
% histograms and other statistics to inspect clusters in relation to task
% events. A variable number of options can be given to modify plots without
% need to go to low level functions.

%% Default options.
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = [600 1000];   end
if ~isfield(param,'treatment'),     param.treatment      = false;        end
if ~isfield(param,'plotcol'),       param.plotcol        = [0.4 0.4 0.4; 0.6 0.1 0.2; 0.0 0.0 0.0]; end
% Rasters
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'square';     end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 1;            end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
if ~isfield(param,'timelim'),       param.timelim        = [ -500 1500]; end
if ~isfield(param,'timelimItiOn'),  param.timelimItiOn   = [    0 2000]; end % Fringe case (no time before 0)
if ~isfield(param,'timelimRwd'),    param.timelimRwd     = [-2000    0]; end % Fringe case (no time after 0)
% PSH
if ~isfield(param,'binSize'),       param.binSize        = 100;          end
if ~isfield(param,'stepSz'),        param.stepSz         = 10;           end
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end
if ~isfield(param,'interval'),      param.interval       = [0 2000];     end

%% Default figure attributes. 
% Raster plot
param.raster.ylabel = {'Trial #'}; % trial label
param.raster.xlabel = {'ms'};   % time label
param.raster.ytick = 0:50:1000; % trial ticks
param.raster.xtick = param.timelim(1):500:param.timelim(2); % time ticks
param.raster.yticklabels = {mat2cell(param.raster.ytick,1)}; % ticks label
param.raster.xticklabels = {mat2cell(param.raster.xtick,1)}; % ticks label

% Waveform
param.wf.ylabel = {'voltage'}; % volt label
param.wf.xlabel = {'time (ms)'}; % time label
param.wf.ytick = -100:50:100; % volt ticks
param.wf.xtick = 0:32:96; % time ticks
param.wf.yticklabels = {''}; % volt labels
param.wf.xticklabels = {mat2cell(floor((param.wf.xtick-32)/32),1)}; % time labels

% PSH
param.psh.ylabel = {'spikes/s'}; % rate label
param.psh.xlabel = {'ms'}; % time label
param.psh.ytick = 0:5:60; % fire rate ticks
param.psh.xtick = (param.interval(1):500:param.interval(2))/param.stepSz; % time ticks
param.psh.yticklabels = {mat2cell(param.psh.ytick,1)}; % rate labels
param.psh.xticklabels = {mat2cell(param.psh.xtick*param.stepSz-500,1)}; % time labels

% ISI Hist
param.isihist.ylabel = {'Rel. prob.'}; % probability label
param.isihist.xlabel = {'ISI (ms)'};   % ISI time label
param.isihist.ytick  = 'auto'; % prob. ticks
param.isihist.xtick  = 0:20:400; % time ticks
param.isihist.yticklabels = {'auto'}; % prob label
param.isihist.xticklabels = {mat2cell(param.isihist.xtick/2,1)}; % time label

% Driftmap
param.driftmap.ylabel = {'tempAmpl'}; % ampl label
param.driftmap.xlabel = {'time (min)'}; % time label
param.driftmap.ytick  = 'auto'; % ampl ticks
param.driftmap.xtick  = 'auto'; % time ticks
param.driftmap.yticklabels = {'auto'}; % ampl label
param.driftmap.xticklabels = {'auto'}; % time label

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

% Individual stats, cluster by cluster
% For each alignment
for a = 1:nalign
    % for each cluster
    for c = 1:length(neurons.(opt.alignto{a}))
        %% Figure
        param.raster.title = ['Aligned to: ', opt.alignto{a}];
        param.raster.subtitle = ['cluster: ', spike.label{c}];
        
        % Special cases for itiON and Rwd, as fringe events
        if strcmp(opt.alignto{a},'itiOn')     
            tlim = param.timelimItiOn; 
            param.psh.xticklabels = {mat2cell(param.psh.xtick*param.stepSz,1)}; % time labels
        elseif strcmp(opt.alignto{a},'rwd')  
            tlim = param.timelimRwd;
        else, tlim = param.timelim;
        end

        % Initialize figure
        fig = figure('visible', param.visible); % switch visibility
        set(fig, 'Position', [0, 0, round(screen.width), round(screen.height)]); % Set fig size as screen
        
        % Event-aligned Spike Raster
        subplot(3,2,1)
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
            prettify(param.raster, a);

        % Plot waveform
        subplot(3,2,2) 
            plot(mean(spike.waveform{c}, 2, "omitnan"), ...
                'LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
            prettify(param.wf)

        % Event-aligned PSH
        subplot(3,2,3)
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
            if upperY < 5, ylim([0 5]); end

        % Plot ISI
        subplot(3,2,4)
            histogram('BinEdges',0:400,'BinCounts',spike.isihist{c}, 'EdgeColor', 'none', 'FaceColor', 'k', ...
                    'Normalization','probability');
            % Prettify
            prettify(param.isihist); 
                xlim([-5 200])

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
        exportgraphics(fig,fullfile(opt.analysis,'genplots',['raster_',opt.alignto{a},'_',spike.label{c},'.png']),'Resolution',600);
        close all
    end
end

end