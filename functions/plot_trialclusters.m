function plot_trialclusters(neurons, events, spike, opt, param)
% Will take trial-long data and plot a series of basic rasters,
% to inspect long dynamics in relation to task events.

%% Default options.
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = [1000 600];   end
if ~isfield(param,'Resolution'),    param.Resolution     = 300;          end
% if ~isfield(param,'treatment'),     param.treatment      = false;        end
if ~isfield(param,'plotcol'),       param.plotcol        = [ .482  .125  .302;
                                                             .220  .161  .420;
                                                             .435  .588  .196];
end
% Rasters
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'lines';      end
if ~isfield(param,'spkWidth'),      param.spkWidth       = 1;            end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
if ~isfield(param,'plotevent'),     param.plotevent      = [1 3 7];      end
% PSH
if ~isfield(param,'binSize'),       param.binSize        = 200;          end
if ~isfield(param,'stepSz'),        param.stepSz         = 10;           end
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end

param.timelim        = [0 6000]; % Hard coded, TODO
param.treatment = opt.treatment; % overrule para.treatment, to deprecate (TODO)

%% Default figure attributes. 
param.raster.title = 'Whole trial';
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
toalignto = opt.alignto;

% How many alignments
nalign = numel(toalignto); 

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
param.levels = 1; % force to integer
param.trial_change = 1;

if param.treatment
    % conds = fieldnames(conditions);
    conds = {'AllTrials'};
    relevant = {1};

    % On top, check existence of defined treatments
    if isfield(opt,'trEvents')
        % how many
        ntreatments = numel(opt.trEvents);

        % have they start/end or just a single application time?
        for i = 1:ntreatments
            % add levels accordingly, e.g. basal/treatment_present/post (+2) or basal/post (+1)
            param.levels = param.levels + numel(events.(opt.trEvents{i}).time{1});
            param.trial_change = [param.trial_change events.(opt.trEvents{i}).trial{1}(2:end)];
        end
    end

% Social yes/no assessment
elseif opt.plotSocial && ~param.treatment
    % Use Social assessment
    conds = fieldnames(events.social);
    % Set really relevant social conditions for an ItiOn alignment
    relevant = {[1:7,10,12,14,15], [], []};
    % in the social assessment, treatments are 1/0. Add levels accordingly,
    param.levels = 2;
    param.trial_change = ones(1,param.levels);

% or conditions (TODO)
elseif opt.useConditions && ~param.treatment
    conds = fieldnames(events.social);
    % Set really relevant social conditions for an ItiOn alignment
    relevant = {[1:7,10,12,14,15], [], []};
    % in the social assessment, treatments are 1/0. Add levels accordingly,
    param.levels = 2;
    param.trial_change = ones(1,param.levels);

end

switch param.levels
    case 1
    param.plotcol      = [0  0  0];
    case 2
    param.plotcol      = [0  0  0;
                         .3 .3 .3];
    case 3
    param.plotcol      = [ .482  .125  .302;
                           .220  .161  .420;
                           .435  .588  .196];     
    otherwise
        for i = 4:3:param.levels
            param.plotcol(i:i+2,:) = param.plotcol(1:3,:);
        end
end

%% Whole trial raster and PSH, cluster by cluster
% for each cluster
for a = 1:nalign
    if strcmp(toalignto{a},'itiOn')
        % add last trial to the level-limits vector
        param.trial_change = [param.trial_change length(neurons.(toalignto{a}){1})];

        % For each indexing condition
        for cc = relevant{a}    
            
            % For each cluster
            for c  = 1:length(neurons.(toalignto{a}))
                jump = 0; % reset zero-trial-index switch

                % Initialize
                param.raster.subtitle = ['cluster: ', spike.label{c}, '. ', conds{cc}];
                fig = figure('visible', param.visible); % switch visibility
                set(fig, 'Position', [0, 0, round(screen.width), round(screen.height)]); % Set fig size as screen
                
                % Trial-long Spike Raster
                subplot(3,1,1)
                trialCounter = 1;
                    for lvl = 1:param.levels-1
                        % prepare range of trials to plot
                        trialrange = param.trial_change(lvl):param.trial_change(lvl+1);
                        condsrange = trialrange;
                        if lvl > 1, trialrange(1) = []; end % remove repeated trial at beginning, when treatments exist
                        if size(conds,1) > 1
                            trialrange = param.trial_change(1):param.trial_change(end);
                            cndidx = events.social.(conds{cc});
                            % restrict trials to those indexed as 
                            if lvl == 1,    condsrange = cndidx; % true
                            elseif lvl==2,  condsrange = ~cndidx; % false
                            end
                            trialCounter = 1; % do not re-start from last trial
                        end
                        
                        if length(trialrange) < 2, continue, end
                        if sum(condsrange) < 1, jump = jump + 1; continue, end % if no trials for this specific conditions, add to off-switch
                        
                        % Now, take spikes to plot
                        toplot = neurons.(toalignto{a}){c}(trialrange);
    
                        % empty trials not indexed by condition but keep trial ordinal
                        toplot(~condsrange) = {[]}; 

                        % raster plot
                        trialCounter = plotRaster(toplot, ... % spikes
                                                trialCounter,               ... % trialCounter
                                               'plotcol',    param.plotcol(lvl,:), ... % color per align (for now)
                                               'spkwidth',   param.spkWidth,       ...
                                               'linelength', param.lineLength,     ...
                                               'plotstyle',  param.plotStyle);
                    end
            
                % Plot requested events (param.plotevent)
                if ~isempty(param.plotevent)
                    color = 'k';
                    for ev = param.plotevent
                        if ev == 1 || 2, color = [.365 .200 .043]; end % stimOn1 || stimOn2 'b'
                        if ev == 3, color = [.345 .067 .231]; end % bhv 'r'
                        if ev == 7, color = [.071 .067 .231]; end % rwd 'g'
                        
                        evidx = cellfun(@(x) x==ev, events.(toalignto{a}).code, 'UniformOutput', 0);
%                         evidx = find([events.(toalignto{a}).code{:}] == ev);
%                         for trial = param.trial_change(1):param.trial_change(end)                            
%                             % Check event presence in trial
%                             index = find([C{:}] == 5);
%                             evidx = cellfun(@find(X==ev), events.(toalignto{a}).code{trial,1}));
%                             % If any, plot
%                             if ~isempty(evidx)
%                                 line(1000*[events.(toalignto{a}).time{trial,1}(evidx) events.(toalignto{a}).time{trial,1}(evidx)], ...
%                                      [trial trial+1], 'Color', color, 'LineWidth', 2)
%                             end
%                         end
                        
                        for trial = param.trial_change(1):param.trial_change(end)  
                            if any(evidx{trial})
                                line(1000*[max(events.(toalignto{a}).time{trial,1}(evidx{trial})) ...
                                           max(events.(toalignto{a}).time{trial,1}(evidx{trial}))], ...
                                     [trial trial+1], 'Color', color, 'LineWidth', 2)
                            end
                        end
                    end
                end
            
                if jump == 3, continue, end % if no trials at any level, cancel figure

                % Prettify
                prettify(param.raster);
                    xlim(param.timelim)
                    ylim([0 param.trial_change(lvl+1)])

                % Trial-long PSH
                jump = 0;
                subplot(3,2,[3,4])
                    for lvl=1:param.levels-1
                        % prepare range of trials to plot
                        trialrange = param.trial_change(lvl):param.trial_change(lvl+1);
                        if lvl > 1, trialrange(1) = []; end % remove repeated trial at beginning, when treatments exist
                        if size(conds,1) > 1
                            trialrange = param.trial_change(1):param.trial_change(end);
                            cndidx = events.social.(conds{cc});
                            % restrict trials to those indexed as 
                            if lvl == 1,    condsrange = cndidx; % true
                            elseif lvl==2,  condsrange = ~cndidx; % false
                            end
                        end
                        
                        if length(trialrange) < 2, continue, end
                        % if no trials for this specific conditions, add to off-switch
                        if sum(condsrange) < 1, jump = jump + 1; continue, end 
    
                        % Now, take spikes to plot
                        toplot = neurons.(toalignto{a}){c}(trialrange);
    
                        % empty trials not indexed by condition but keep trial ordinal
                        toplot(~condsrange) = {[]}; 

                        % Plot
                        upperY(lvl) = plotPSTH(toplot, ... % spikes
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
                    maxY = max(upperY); if maxY <= 5, maxY = 6; end % force a minimum y-axis scale
                    ylim([0 maxY*1.1]);
            
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
                exportgraphics(fig,fullfile(opt.analysis,'genplots', ...
                    ['fulltrial_', spike.label{c}, '_', conds{cc}, '.png']), ...
                    'Resolution', param.Resolution);
                close all
            end
        end
    end
end
end