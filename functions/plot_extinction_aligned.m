function plot_extinction_aligned(neurons, events, spike, conditions, opt, param)
% Will take trial-long data and plot a series of basic rasters,
% to inspect long dynamics in relation to task events.

%% Default options.
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = 'adaptive';   end
if ~isfield(param,'Resolution'),    param.Resolution     = 300;          end
if ~isfield(param,'treatment'),     param.treatment      = true;         end
if ~isfield(param,'nBlocks'),       param.nBlocks        = 8;            end
if ~isfield(param,'binSize'),       param.binSize        = 500;          end
if ~isfield(param,'stepSz'),        param.stepSz         = 50;           end
if ~isfield(param,'timelim'),       param.timelim        = [-2000 3000]; end
if ~isfield(param,'trial2plot'),    param.trial2plot     = 'correct';    end
if ~isfield(param,'FS2plot'),       param.FS2plot        = false;        end

% Rasters
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'lines';      end
if ~isfield(param,'spkWidth'),      param.spkWidth       = .5;           end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end
if ~isfield(param,'plotevent'),     param.plotevent      = [1 3 7];      end

% PSH
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end

%% Default figure attributes. 
% Raster plot
param.raster.ylabel = {'Trial #'}; % trial label
param.raster.ytick = 0:50:1000; % trial ticks
param.raster.yticklabels = {mat2cell(param.raster.ytick,1)}; % ticks label
param.raster.xlabel = {'time from alignment (s)'};   % time label
param.raster.xtick = param.timelim(1):1000:param.timelim(2); % time ticks
param.raster.xticklabels = {mat2cell(param.raster.xtick/1000,1)}; % ticks label

% PSH
param.psh.ylabel = {'spikes/s'}; % rate label
param.psh.ytick = 0:5:60; % fire rate ticks
param.psh.yticklabels = {mat2cell(param.psh.ytick,1)}; % rate labels
param.psh.xlabel = {'time from alignment (s)'}; % time label
param.psh.xtick = (0:1000:diff([param.timelim(1) param.timelim(2)]))/param.stepSz; % time ticks
param.psh.xticklabels = {mat2cell((param.psh.xtick-40)*param.stepSz/1000,1)}; % time labels

% Figure size
if strcmpi('adaptive', param.size)
    param.screen.size = get(0, 'ScreenSize');  
    param.screen.width =  param.screen.size(3);
    param.screen.height = param.screen.size(4);
else
    param.screen.width =  param.size(1); 
    param.screen.height = param.size(2);
end

%% Initialize
toalignto = opt.alignto;
param.levels{1}         = 1; % force cell integer
param.trial_change{1}   = 1; % force cell integer
param.cond = 'tr1';
param.plotcol = [0 0 1; 1 0 0; 1 0 1];
blks = cell(2, param.nBlocks);
blocks.Pos = {'A' 'B' 'A' 'B' 'A' 'C' 'A' 'C'};
blocks.Amb = {'A' 'B' 'A' 'B' 'C' 'B' 'A' 'C'};
blocks.Neg = {'A' 'B' 'A' 'B' 'C' 'B' 'C' 'A'};

%% Prepare treatments
% add levels accordingly, e.g. basal/treatment_present/post (+2) or basal/post (+1)
param.levels{1} = numel(events.(opt.trEvents{1}).trial{1});

% For blocks, add the first and last trials to complete the 8 blocks
if events.(opt.trEvents{1}).trial{1}(1) == 0
    param.trial_change{1} = [1 events.(opt.trEvents{1}).trial{1}(2:end) length(neurons.itiOn{1})];
else
    param.trial_change{1} = [1 events.(opt.trEvents{1}).trial{1}(1:end) length(neurons.itiOn{1})];
end

%% Aligned raster and PSH, cluster by cluster
if strcmpi(param.trial2plot, 'allInitiated'), param.maxalign = 2;
else, param.maxalign = length(toalignto)-1; end

for j = 1:param.maxalign %(e.g. Stim1 (ini), Stim2 (choice), bhv (bhv2) ...)
    % Grab neuron data to plot
    neuronSet{j} = neurons.(toalignto{j+1});
    
    % Small modification to fit the field name
    if strcmp(toalignto{j+1},'bhv2'), toalignto{j+1} = 'bhv'; end

    % Now we collect the events
    eventsSet{j} = events.(toalignto{j+1}); 
    
    % Let's empty those of no interest, and collect the trial range
    if strcmp(param.trial2plot, 'correct') 
        eventsSet{j}.code(conditions.aborted | conditions.incorrect | conditions.omission) = {[]};
    elseif strcmp(param.trial2plot, 'incorrect')
        eventsSet{j}.code(conditions.aborted | conditions.correct | conditions.omission) = {[]};
    elseif strcmp(param.trial2plot, 'omission')
        eventsSet{j}.code(conditions.aborted | conditions.correct | conditions.incorrect) = {[]};
    elseif strcmp(param.trial2plot, 'allInitiated')
        eventsSet{j}.code(logical(conditions.aborted)) = {[]};
    end
    trialrange = 1:length(conditions.correct);
end

% For each cluster
for c = 1:length(neuronSet{1})
    % Empty toplot cell array
    toplot = cell(param.maxalign, param.nBlocks);

    % cluster label and waveforms
    param.label = spike.label{c};
    spikewf = mean(spike.waveform{1,c},2);
    param.raster.title = [param.label, ' ', param.trial2plot];
    if param.FS2plot, param.raster.title = [param.raster.title, ' FS&NS'];
    else,             param.raster.title = [param.raster.title, ' NS']; end

    % For each Block
    for i = 1:param.nBlocks
        % Find out block type
        blks{1,i} = blocks.(conditions.type){i};     
        
        % Get the range trial borders
        blks{2,i} = [param.trial_change{1}(i) param.trial_change{1}(i+1)];

        % Here, we find trials for the current block range
        cndidx = ismember(trialrange, param.trial_change{1}(i):param.trial_change{1}(i+1));
        
        % And we keep ONLY those for Novel Stimuli (tr2)
        if ~param.FS2plot, FSidx = ismember(trialrange, events.tr2.trial{1});
        else, FSidx = ones(1,length(trialrange)); end
        cndidx = trialrange(FSidx & cndidx);
        
        % We find the complementary set of trials (non-block (and FS, if so))
        notcndidx = trialrange(setdiff(1:end,cndidx));

        % And we select trials from cluster
        for j = 1:param.maxalign       
            % Take all trials first
            toplot{j,i} = neuronSet{j}{c};

            % Remove set of trials of non-interest  
            toplot{j,i}(notcndidx) = {[]}; 

            % Also, empty trials with conditions of no interest
            if strcmp(param.trial2plot, 'correct')
                toplot{j,i}(conditions.aborted | conditions.incorrect | conditions.omission) = {[]};
            elseif strcmp(param.trial2plot, 'incorrect')
                toplot{j,i}(conditions.aborted | conditions.correct | conditions.omission) = {[]};
            elseif strcmp(param.trial2plot, 'omission')
                toplot{j,i}(conditions.aborted | conditions.correct | conditions.incorrect) = {[]};
            elseif strcmp(param.trial2plot, 'allInitiated')
                toplot{j,i}(logical(conditions.aborted)) = {[]};
            end
        end
    end

    %% FIGURE
    plot_blocks(toplot, spikewf, blks, eventsSet, param, opt)
    
    close all hidden
    close all force
end
end % Main function

function plot_blocks(toplot,spikewf, blks, eventsSet, param, opt)
    color = cell(1,param.nBlocks);
    upperY = nan(1,param.nBlocks);

    fig = figure('visible', param.visible); % switch visibility
    set(fig, 'Position', [0, 0, round(param.screen.width), round(param.screen.height)]); % Set fig size as screen

    % Trial-long Spike Rasters
    for j = 1:param.maxalign
        sbpl = subplot(2,param.maxalign,j);
        trialCounter = 1; % always plot from first trial
    
        for i = 1:param.nBlocks
            if strcmp(blks(1,i),'A'), color{i} = param.plotcol(1,:); end
            if strcmp(blks(1,i),'B'), color{i} = param.plotcol(2,:); end
            if strcmp(blks(1,i),'C'), color{i} = param.plotcol(3,:); end
    
            % raster
            plotRaster(toplot{j,i}, ... % spikes
                        trialCounter,               ... % trialCounter
                       'plotcol',    color{i}, ... % color per align (for now)
                       'spkwidth',   1,       ...
                       'linelength', param.lineLength,     ...
                       'plotstyle',  param.plotStyle);
            
            yline(blks{2,i}(2)+1, '-', blks{1,i}, 'Color', color{i}, 'LabelVerticalAlignment', 'bottom',  'LineWidth', 1);
        end
    
        % Plot events
        plot_events(param.plotevent, eventsSet{j});
    
        % Prettify
        prettify(param.raster);
            xlim([-1000 param.timelim(2)])
            ylim([0 blks{2,i}(2)+1])
            if j == 1
                xline(0,'--k')
                subtitle('Aligned to Ini Stimulus')
            elseif j == 2       
                sbpl.YAxis.Visible = 'off';    
                subtitle('Aligned to Choice Stimulus')
            elseif j == 3
                sbpl.YAxis.Visible = 'off'; 
                subtitle('Aligned to Choice behavior')
            end

    % PSH plots
        subplot(2,param.maxalign,j+param.maxalign)
        for i = 1:param.nBlocks
            upperY(1,i) = plotPSTH(toplot{j,i}, ... % spikes
                        param.stepSz,   ... % stepSz
                        param.binSize,  ...  % binSize
                        param.timelim,  ... % interval
                        param.smpRate,  ...  % samples per second in the feeded data
                        'plotcol',      color{i},...
                        'meanline',     '-',...
                        'smoothplot',   true, ...
                        'erralpha',     0.2);
        end
    
        % Prettify
        prettify(param.psh)
            xline(param.psh.xtick(find(param.psh.xticklabels{1}{1} == 0)),'--k');
            maxY = max(upperY); if maxY <= 5, maxY = 5; end % force a minimum y-axis scale
            ylim([0 maxY*1.1]);
            xlim([param.psh.xtick(2) param.psh.xtick(end)])
        if j > 1,  sbpl.YAxis.Visible = 'off'; end
    end

    % inset waveform 
    axes('Position',[.85 .34 .05 .1])
        plot(spikewf, 'Color', 'k', 'LineWidth', 2), box off
        xlabel('ms'), ylabel('Ampl. (a.u.)')
        xlim([-2 96])
        xticks([0 32 64 96]); xticklabels({'-1' '0' '1' '2'});
        yticks([]); yticklabels({});

    % Save figure per alignment&cluster    
    sgtitle(param.raster.title)
    exportgraphics(fig, fullfile(opt.analysis,'newplots', ...
                        ['Aligned_', param.raster.title, '.png']), ...
                        'Resolution', param.Resolution);
end

function plot_events(evidx, events)
% Plot requested events (param.plotevent)
for ev = evidx
    if ev == 2, color = 'k'; end % StimOn2 
    if ev == 3, color = [0 0.6 0]; end % bhv 'g'
    if ev == 7, color = 'k'; end % rwd 'g'
    idx = find(~cellfun(@isempty, events.code))';

    for i = idx
        evnt = find(events.code{i,1}==ev);
        if ~isempty(evnt)
            if ev == 3 && length(evnt)>1
                evtime1 = min(events.time{i,1}(evnt))*1000;
                evtime2 = max(events.time{i,1}(evnt))*1000;
    
                line([evtime1 evtime1], [i i+1], 'Color', color, 'LineWidth', 1)
    
                if evtime2~=0
                    line([evtime2 evtime2]-1000, [i i+1], 'Color', color, 'LineWidth', 1) % Fixing behavior artificial delay of 1 second (to Stim2!)
                else
                    line([evtime2 evtime2], [i i+1], 'Color', color, 'LineWidth', 1) % Fixing behavior artificial delay of 1 second (to Stim2!)
                end
            else
                evtime = events.time{i,1}(evnt)*1000;
                line([evtime evtime], [i i+1], 'Color', color, 'LineWidth', 1)
            end
        end
    end
end
end