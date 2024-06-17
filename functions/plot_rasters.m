function plot_rasters(neurons, events, spike, opt, param)
% Will take neuron-trial data and plot a series of basic rasters,
% histograms and other statistics to inspect clusters in relation to task
% events. A variable number of options can be given to modify plots without
% need to go to low level functions.

%% Default options.
if ~isfield(param,'res'),           param.res            = false;        end
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = [600 1000];   end
if ~isfield(param,'treatment'),     param.treatment      = false;        end
if ~isfield(param,'genstats'),      param.genstats       = true;         end
if ~isfield(param,'plotcol'),       param.plotcol        = zeros(1,3);   end
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'square';     end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 3;            end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
if ~isfield(param,'timelim'),       param.timelim        = [-500 2000];  end
if ~isfield(param,'timelimItiOn'),  param.timelimItiOn   = [0 2500];     end
if ~isfield(param,'timelimRwd'),    param.timelimRwd     = [-2000 500];  end
if ~isfield(param,'binSize'),       param.binSize        = 100;          end
if ~isfield(param,'stepSz'),        param.stepSz         = 10;           end
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end
if ~isfield(param,'interval'),      param.interval       = [0 (param.timelim(2)-param.timelim(1))/10];  end

%% Default figure attributes. 
% Raster plots
param.raster.ylabel = {'Trial #'}; % trial label
param.raster.xlabel = {'ms'};   % time label
param.raster.ytick = 0:50:1000; % trial ticks
param.raster.xtick = param.timelim(1):250:param.timelim(2); % time ticks
param.raster.yticklabels = {mat2cell(param.raster.ytick,1)}; % ticks label
param.raster.xticklabels = {mat2cell(param.raster.xtick,1)}; % ticks label

% Waveform
param.wf.ylabel = {'voltage (uv)'}; % trial label
param.wf.xlabel = {'time (ms)'};   % time label
param.wf.ytick = -100:10:100; % trial ticks
param.wf.xtick = 0:32:96; % time ticks
param.wf.yticklabels = {mat2cell(param.wf.ytick,1)}; % ticks label
param.wf.xticklabels = {mat2cell(floor(param.wf.xtick/32),1)}; % ticks label

% PSH
param.psh.ylabel = {'spikes/s'}; % trial label
param.psh.xlabel = {'ms'};   % time label
param.psh.ytick = 0:5:60; % trial ticks
param.psh.xtick = 0:param.stepSz*5:param.interval(2)/param.stepSz; % time ticks
param.psh.yticklabels = {mat2cell(param.psh.ytick,1)}; % ticks label
param.psh.xticklabels = {mat2cell(param.psh.xtick*param.stepSz,1)}; % ticks label

% ISIHist
param.isihist.ylabel = {'P (pdf est.)'}; % trial label
param.isihist.xlabel = {'ISI (ms)'};   % time label
param.isihist.ytick  = 'auto'; % time ticks
param.isihist.xtick  = 0:10:400; % time ticks
param.isihist.xticklabels = {mat2cell(param.isihist.xtick/2,1)}; % ticks label
param.isihist.yticklabels = {'auto'}; % ticks label

% Driftmap
param.driftmap.ylabel = {'tempAmpl'}; % trial label
param.driftmap.xlabel = {'time (min)'};   % time label
param.driftmap.ytick  = 'auto'; % time ticks
param.driftmap.xtick  = 'auto'; % time ticks
param.driftmap.yticklabels = {'auto'}; % ticks label
param.driftmap.xticklabels = {'auto'}; % ticks label

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
% For each alignment
for a = 1:nalign
    % for each cluster
    for c = 1:length(neurons.(opt.alignto{a}))
         if param.genstats
            % Prepare figure
            param.raster.title = ['Aligned to: ', opt.alignto{a}];
            param.raster.subtitle = ['cluster: ', spike.label{c}];
            
            % Special cases for itiON and Rwd, as fringe events
            if strcmp(opt.alignto{a},'itiOn'),      tlim = param.timelimItiOn;
            elseif strcmp(opt.alignto{a},'rwd'),    tlim = param.timelimRwd;
            else, tlim = param.timelim;
            end

            % Initialize figure
            fig = figure('visible', param.visible); % switch visibility
            set(fig, 'Position', [0, 0, round(screen.width), round(screen.height)]); % Set fig size as screen
            
            %% Event-aligned Spike Raster
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

            %% Plot waveform
            subplot(3,2,2) 
                plot(mean(spike.waveform{c}, 2, "omitnan"),'LineWidth',2,'Color','b','LineStyle','-');
            
            prettify(param.wf)
                ylim('auto')

            %% Event-aligned PSH
            subplot(3,2,3)
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

            %% Plot ISI
            subplot(3,2,4)
                histogram(spike.isihist{c}, 400,'EdgeColor','none','Normalization','pdf');
            
            % Prettify
            prettify(param.isihist); 
                xlim([-5 100])

            %% Plot a single 'driftmap'
            subplot(3,2,[5,6])
                scatter(spike.timestamp{c}/60000, spike.ampl{c}, 5, "black", "filled");

                if isfield(events,'na3') % If there is a na3 treatment field
                    rectangle('Position', [events.na3.time{1}(1)/60, min(spike.ampl{c})*0.9, ... % [x, y ...
                                          (events.na3.time{1}(2)-events.na3.time{1}(1))/60, (max(spike.ampl{c})*1.1-min(spike.ampl{c})*0.9)], ... % ... w, h]
                                          'FaceColor',[.5 .5 .5 .35], 'LineStyle', 'none'); % [r g b alpha], no line.
                end

            % Prettify
            prettify(param.driftmap);
                ylim([min(spike.ampl{c})*0.9 max(spike.ampl{c})*1.1]);
                xlim([-.2 spike.timestamp{c}(end)/60000+.2]);

            %% Save figure per alignment&cluster
            exportgraphics(fig,fullfile(opt.analysis,'genplots',['raster_',opt.alignto{a},'_',spike.label{c},'.png']),'Resolution',600);
            close all

         end
    end
end