function plot_alignedclusters(neurons, events, spike, opt, param)
% Will take neuron-trial data and plot a series of basic rasters,
% histograms and other statistics to inspect clusters in relation to task
% events. A variable number of options can be given to modify plots without
% need to go to low level functions.

%% Default options.
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = [600 1000];   end
if ~isfield(param,'treatment'),     param.treatment      = false;        end
if ~isfield(param,'plotcol'),       param.plotcol        = [ 0  0  0;
                                                            .6 .6 .6;
                                                            .3 .3 .3];
end
if ~isfield(param,'baseline'),      param.baseline       = 500;          end
if ~isfield(param,'post'),          param.post           = 2500;         end
% Rasters
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'lines';      end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 1;            end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
if ~isfield(param,'timelim'),       param.timelim        = [-param.baseline param.post]; end
% PSH
if ~isfield(param,'binSize'),       param.binSize        = 100;          end
if ~isfield(param,'stepSz'),        param.stepSz         = 10;           end
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end
if ~isfield(param,'interval'),      param.interval       = [0 diff(param.timelim)]; end

%% Default figure attributes. 
% Raster plot
param.raster.ylabel = {'Trial #'}; % trial label
param.raster.xlabel = {'ms'};   % time label
param.raster.ytick = 0:50:1000; % trial ticks
param.raster.xtick = param.timelim(1):param.baseline:param.timelim(2); % time ticks
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
param.psh.xtick = (param.interval(1):param.baseline:param.interval(2))/param.stepSz; % time ticks
param.psh.yticklabels = {mat2cell(param.psh.ytick,1)}; % rate labels
param.psh.xticklabels = {mat2cell(param.psh.xtick*param.stepSz-param.baseline,1)}; % time labels

% ISI Hist
param.isihist.ylabel = {'Rel. prob.'}; % probability label
param.isihist.xlabel = {'ISI (ms)'};   % ISI time label
param.isihist.ytick  = 'auto'; % prob. ticks
param.isihist.xtick  = 0:20:400; % time ticks
param.isihist.yticklabels = {'auto'}; % prob label
param.isihist.xticklabels = {mat2cell(param.isihist.xtick/2,1)}; % time label

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
% For each alignment, skipping itiOn
for a = 2:nalign
    % for each cluster
    for c = 1:length(neurons.(opt.alignto{a}))
        %% Figure
        param.raster.title = ['Aligned to: ', opt.alignto{a}];
        param.raster.subtitle = ['cluster: ', spike.label{c}];

        % Initialize figure
        fig = figure('visible', param.visible); % switch visibility
        set(fig, 'Position', [0, 0, round(screen.width), round(screen.height)]); % Set fig size as screen
        
        % Event-aligned Spike Raster
        subplot(2,2,1)
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
                                       'plotstyle',  param.plotStyle,      ...
                                       'timelim',    param.timelim);
            end

            % Plot requested events (param.plotevent)
            if ~isempty(param.plotevent)
                for trial = 1:size(events.(opt.alignto{a}).code,1)
                    for ev = param.plotevent
                        % Check event presence in trial
                        evidx = find(events.(opt.alignto{a}).code{trial,1} == ev);
                        % If any, plot
                        if ~isempty(evidx)
                            color = [1 1 1]; 
                            if ev == 1 %&& strcmpi(opt.alignto{a},'stimOn1')
                                color = 'b'; end % stimOn1
                            if ev == 3 %&& strcmpi(opt.alignto{a},'rwd')
                                color = 'r'; end % bhv
                            if ev == 7 %&& strcmpi(opt.alignto{a},'rwd')
                                color = 'g'; end % rwd
                            line(1000*[events.(opt.alignto{a}).time{trial,1}(evidx) events.(opt.alignto{a}).time{trial,1}(evidx)], ...
                                 [trial trial+1], 'Color', color, 'LineWidth', 1)
                        end
                    end
                end
            end

            % Prettify
            prettify(param.raster);
                ylim([0 trialCounter+1])

        % Plot waveform
        subplot(2,2,2) 
            plot(mean(spike.waveform{c}, 2, "omitnan"), ...
                'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
            prettify(param.wf)

        % Event-aligned PSH
        subplot(2,2,3)
            for lvl=1:param.levels
                % prepare range of trials to plot
                trialrange = param.trial_change(lvl):param.trial_change(lvl+1);
                    if lvl > 1, trialrange(1) = []; end % remove repeated trial at beginning, when treatments exist
                    if length(trialrange) < 2, continue, end
                % plot
                upperY = plotPSTH(neurons.(opt.alignto{a}){c}(trialrange), ... % spikes
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
                xline(param.psh.xtick(find(param.psh.xticklabels{1}{1} == 0)),'--k');
                if upperY <= 5, ylim([0 5]); end

        % Plot ISI
        subplot(2,2,4)
            histogram('BinEdges', 0:400, 'BinCounts', spike.isihist{c}, ...
                      'EdgeColor', 'none', 'FaceColor', 'k', ...
                      'Normalization','probability');
            % Prettify
            prettify(param.isihist); 
                xlim([-5 200])

        % Save figure per alignment&cluster
        exportgraphics(fig,fullfile(opt.analysis,'genplots',['raster_',opt.alignto{a},'_',spike.label{c},'.png']),'Resolution',600);
        close all
    end
end

end