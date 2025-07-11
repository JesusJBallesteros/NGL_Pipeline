function plot_alignedclusters(neurons, events, spike, conditions, opt, param)
% Will take neuron-trial data and plot a series of basic rasters,
% histograms and other statistics to inspect clusters in relation to task
% events. A variable number of options can be given to modify plots without
% need to go to low level functions.

%% Default options.
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = [600 1000];   end
if ~isfield(param,'Resolution'),    param.Resolution     = 300;          end
if ~isfield(param,'treatment'),     param.treatment      = false;        end
if ~isfield(param,'plotcol'),       param.plotcol        = [ .482  .125  .302;
                                                             .220  .161  .420;
                                                             .435  .588  .196];
end
if ~isfield(param,'baseline'),      param.baseline       = 500;          end
if ~isfield(param,'post'),          param.post           = 2500;         end
if ~isfield(param,'plotevent'),     param.plotevent      = 1;            end
% Rasters
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'lines';      end
if ~isfield(param,'spkWidth'),      param.spkWidth       = .5;            end
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
param.psh.xticklabels = {mat2cell((param.psh.xtick*param.stepSz)-param.baseline,1)}; % time labels

% ISI Hist
param.isihist.ylabel = {'Rel. prob.'}; % probability label
param.isihist.xlabel = {'ISI (ms)'};   % ISI time label
param.isihist.ytick  = 'auto'; % prob. ticks
param.isihist.xtick  = 0:20:400; % time ticks
param.isihist.yticklabels = {'auto'}; % prob label
param.isihist.xticklabels = {mat2cell(param.isihist.xtick/2,1)}; % time label

% Figure size
if strcmpi('adaptive', param.size)
    screen.size = get(0, 'ScreenSize');  
    screen.width =  screen.size(3);
    screen.height = screen.size(4);
else
    screen.width =  param.size(1); 
    screen.height = param.size(2);
end

%% Initialize stuff
toalignto = opt.alignto;

% How many alignments
nalign = numel(toalignto); 

param.levels{1}         = 1; % force cell integer
param.trial_change{1}   = 1; % force cell integer
ntreatments             = 1;
conds                   = {};
bhv = 0;

%% Prepare treatments
if param.treatment
    % conds = fieldnames(conditions);
    conds = {'AllTrials'};
    relevant = {1};

    % On top, check existence of defined treatments
    if isfield(opt,'trEvents')
        % how many
        ntreatments = numel(opt.trEvents);

        % for each
        for i = 1:ntreatments
            % is it an inter-trial appearance (trial subfield)?
            if isfield(events.(opt.trEvents{i}), 'trial')
                % add levels accordingly, e.g. basal/treatment_present/post (+2) or basal/post (+1)
                param.levels{i} = param.levels{1}(1) + numel(events.(opt.trEvents{i}).time{1});
                param.trial_change{i} = [param.trial_change{1}(1) events.(opt.trEvents{i}).trial{1}(2:end)];
            end

            % Small tweak for 'tr2' case in Extintion Arena, here it
            % signifies that Novel Stimulus was shown in the trial,
            % rather than a block/level change in the paradigm.
            if param.levels{i} > 9, param.levels{i} = 2; end
        end
    end
end 

% Social/Conditions assessment
if opt.plotSocial
    conds = [conds; fieldnames(events.social)];
    % Set really relevant social conditions for an ItiOn alignment
    relevant = {[1:7,10,12,14,15], [], []};
    % in the social assessment, treatments are 1/0. Add levels accordingly,
    param.levels{3} = 2;
    param.trial_change{3} = ones(1,param.levels{i});

elseif opt.useConditions
    conds = [conds; fieldnames(conditions)];
    % Set relevant conditions for alignments (leave alltrials out)
    relevant = {[], 2:length(conds), 2:length(conds)};
end

% if param.treatment
%     % conds = fieldnames(conditions);
%     conds = {'AllTrials'};
%     relevant = {1};
%     % On top, check existence of defined treatments
%     if isfield(opt,'trEvents')
%         % how many
%         ntreatments = numel(opt.trEvents);
%         relevant = {1, 1, 1, 1, 1};
%         % have they start/end or just a single application time?
%         for i = 1:ntreatments
%             % add levels accordingly, e.g. basal/treatment_present/post (+2) or basal/post (+1)
%             param.levels{i} = param.levels{i} + numel(events.(opt.trEvents{i}).time{1});
%             param.trial_change = [param.trial_change events.(opt.trEvents{i}).trial{1}(2:end)];
%         end
%     end
% % Social yes/no assessment
% elseif opt.plotSocial && ~param.treatment
%     % Use Social assessment
%     conds = fieldnames(events.social);
%     % Set really relevant social conditions for StimOn, rwd, alignment
%     relevant = {[], [1:7,10,12,14,15], [1:7,10,12,14,15]};
%     % in the social assessment, treatments are 1/0. Add levels accordingly,
%     param.levels{i} = 2;
%     param.trial_change = ones(1,param.levels{i});
%     param.plotcol      = [0  0  0;
%                           1  0  0];
% % or conditions (TODO)
% elseif opt.useConditions && ~param.treatment
%     conds = fieldnames(events.social);
%     % Set really relevant social conditions for an ItiOn alignment
%     relevant = {[], [], []}; % (TODO)
%     % in the social assessment, treatments are 1/0. Add levels accordingly,
%     param.levels{i} = 2;
%     param.trial_change = ones(1,param.levels{i});
% end

switch param.levels{i}
    case 1
    param.plotcol      = [0  0  0];
    case 2
    param.plotcol      = [0  0  0;
                         .3 .3 .3];
    case 3
    param.plotcol      = [ 0  0  0;
                          .6 .6 .6;
                          .3 .3 .3];
    otherwise
        for i = 4:3:param.levels{i}
            param.plotcol(i:i+2,:) = param.plotcol(1:3,:);
        end
end

% add last trial to the level-limits vector
param.trial_change = [param.trial_change length(neurons.(toalignto{1}){1})];

%% Figure
% For each alignment.
for a = 1:nalign
    % Skip itiOn aligment
    if strcmp(toalignto{a},'itiOn'), continue, end 
    if strcmp(toalignto{a},'bhv'), bhv = bhv + 1; end 
    
    % Alignment title
    param.raster.title = ['Aligned to: ', toalignto{a}];
        if bhv == 2, param.raster.title = ['Aligned to: bhv2']; end
    
    % For each indexing condition
    for cc = relevant{a}

        % For each cluster
        for c = 1:length(neurons.(toalignto{a}))
            jump = 0; % reset zero-trial-index switch
            
            % Initialize figure
            param.raster.subtitle = ['cluster: ', spike.label{c}, '. ', conds{cc}];
            fig = figure('visible', param.visible); % switch visibility
            set(fig, 'Position', [0, 0, round(screen.width), round(screen.height)]); % Set fig size as screen
            
            % Event-aligned Spike Raster
            subplot(2,2,1)
            trialCounter = 1;
                for lvl = 1:param.levels{i}

                    % For no condition, skip last level
                    if param.treatment && lvl==param.levels{i}(end), continue, end
                    
                    % prepare range of trials to plot
                    trialrange = param.trial_change{a}(lvl):param.trial_change{a}(lvl+1);
                    condsrange = trialrange;
                    
                    % remove repeated trial at beginning, when treatments exist
                    if lvl > 1, trialrange(1) = []; end 
                    
                    % When plotting trials with condition true/false
                    if size(conds,1) > 1
                        trialrange = param.trial_change{a}(1):param.trial_change{a}(end);
                                                
                        if opt.plotSocial || opt.useConditions % Here, we plot trials true for a given event
                            cndidx = conditions.(conds{cc});
                           % cndidx = events.xxxxx.(conds{cc});  TODO
                           
                           % restrict trials to those indexed above 
                           if     lvl == 1,  condsrange = cndidx;  % use true       
                           elseif lvl == 2,  condsrange = ~cndidx; % use false
                           end
                        end

                        trialCounter = 1; % always plot from first trial
                    end
                    
                    if sum(condsrange) < 1, jump = jump + 1; continue, end % if no trials for this specific conditions, add to off-switch
                    if length(trialrange) < 2, continue, end
                       
                    % Now, take spikes to plot
                    toplot = neurons.(toalignto{a}){c}(trialrange);

                    % empty trials not indexed by condition but keep trial ordinal
                    if opt.plotSocial || opt.useConditions
                        toplot(~condsrange) = {[]}; 
                    end

                    % plot
                    trialCounter = plotRaster(toplot,       ... % spikes
                                trialCounter,               ... % trialCounter
                               'plotcol',    param.plotcol(lvl,:), ... % color per align (for now)
                               'spkwidth',   param.spkWidth,       ...
                               'linelength', param.lineLength,     ...
                               'plotstyle',  param.plotStyle,      ...
                               'timelim',    param.timelim);
                end

            % Plot requested events (param.plotevent)
            if ~isempty(param.plotevent)
                color = 'k';
                for ev = param.plotevent
                    if ev == 1 || 2, color = 'b'; end % stimOn1 || stimOn2
                    if ev == 3, color = 'r'; end % bhv
                    if ev == 7, color = 'g'; end % rwd
                    
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
                    
                    for trial = param.trial_change{a}(1):param.trial_change{a}(end)  
                        if any(evidx{trial})
                            if bhv == 1
                                line(1000*[min(events.(toalignto{a}).time{trial,1}(evidx{trial})) ...
                                           min(events.(toalignto{a}).time{trial,1}(evidx{trial}))], ...
                                     [trial trial+1], 'Color', color, 'LineWidth', 2)
                            elseif bhv == 2
                                line(1000*[max(events.(toalignto{a}).time{trial,1}(evidx{trial})) ...
                                           max(events.(toalignto{a}).time{trial,1}(evidx{trial}))], ...
                                     [trial trial+1], 'Color', color, 'LineWidth', 2)
                            end
                        end
                    end
                end
            end
    
            if jump == 3, continue, end % if no trials at any level, cancel figure
                
            % Prettify
            prettify(param.raster);
                if param.treatment, ylimmax = param.trial_change(lvl);
                else,  ylimmax = param.trial_change(lvl+1); end
                ylim([0 ylimmax])
            
            % Plot waveform
            subplot(2,2,2) 
            if isfield(spike,'waveform')
                plot(mean(spike.waveform{c}, 2, "omitnan"), ...
                    'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
                prettify(param.wf)
            end
        
            % Event-aligned PSH
            jump = 0;
            subplot(2,2,3)
                for lvl = 1:param.levels{i}-1
                    % prepare range of trials to plot
                    trialrange = param.trial_change{a}(lvl):param.trial_change{a}(lvl+1);
                    if lvl > 1, trialrange(1) = []; end % remove repeated trial at beginning, when treatments exist
                    if size(conds,1) > 1
                        trialrange = param.trial_change{a}(1):param.trial_change{a}(end);
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
    
                    % plot
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
    
            % Plot ISI
            subplot(2,2,4)
                histogram('BinEdges', 0:400, 'BinCounts', spike.isihist{c}, ...
                          'EdgeColor', 'none', 'FaceColor', 'k', ...
                          'Normalization','probability');
            % Prettify
            prettify(param.isihist); 
                xlim([-5 200])
            
            % Save figure per alignment&cluster
            if bhv == 1, algmnt = toalignto{a};
            elseif bhv == 2, algmnt = 'bhv2'; end

            exportgraphics(fig,fullfile(opt.analysis,'genplots', ...
                            ['raster_', algmnt, '_', spike.label{c}, '_', conds{cc}, '.png']), ...
                            'Resolution', param.Resolution);
            close all
        end
    end
end
end