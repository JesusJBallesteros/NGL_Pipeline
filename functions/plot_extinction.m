function plot_extinction(neurons, events, spike, conditions, opt, param)
% Will take trial-long data and plot a series of basic rasters,
% to inspect long dynamics in relation to task events.

%% Default options.
if ~isfield(param,'visible'),       param.visible        = 'off';        end
if ~isfield(param,'size'),          param.size           = 'adaptive';   end
if ~isfield(param,'Resolution'),    param.Resolution     = 300;          end
if ~isfield(param,'treatment'),     param.treatment      = true;         end
if ~isfield(param,'trial2plot'),    param.trial2plot     = 'correct';    end
if ~isfield(param,'nBlocks'),       param.nBlocks        = 8;            end
if ~isfield(param,'binSize'),       param.binSize        = 500;          end
if ~isfield(param,'stepSz'),        param.stepSz         = 50;           end
if ~isfield(param,'timelim'),       param.timelim        = [-2000 10000]; end
if ~isfield(param,'plotevent'),     param.plotevent      = [1 3 7];      end
if ~isfield(param,'FS2plot'),       param.FS2plot        = false;        end

% Rasters
if ~isfield(param,'plotStyle'),     param.plotStyle      = 'lines';      end
if ~isfield(param,'spkWidth'),      param.spkWidth       = .5;           end
if ~isfield(param,'lineLength'),    param.lineLength     = 1;            end

% PSH
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end

%% Default figure attributes. 
% Raster plot
param.raster.ylabel = {'Trial #'}; % trial label
param.raster.ytick = 0:50:1000; % trial ticks
param.raster.yticklabels = {mat2cell(param.raster.ytick,1)}; % ticks label
param.raster.xlabel = {'time from trial start (s)'};   % time label
param.raster.xtick = param.timelim(1):1000:param.timelim(2); % time ticks
param.raster.xticklabels = {mat2cell(param.raster.xtick/1000,1)}; % ticks label

% PSH
param.psh.ylabel = {'spikes/s'}; % rate label
param.psh.ytick = 0:5:60; % fire rate ticks
param.psh.yticklabels = {mat2cell(param.psh.ytick,1)}; % rate labels
param.psh.xlabel = {'time from trial start (s)'}; % time label
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
param.trial_change{1}   = 1; % force cell integer
blks = cell(2, param.nBlocks);
blocks.Pos = {'A' 'B' 'A' 'B' 'A' 'C' 'A' 'C'};
blocks.Amb = {'A' 'B' 'A' 'B' 'C' 'B' 'A' 'C'};
blocks.Neg = {'A' 'B' 'A' 'B' 'C' 'B' 'C' 'A'};

%% Prepare treatments
% add levels accordingly, e.g. basal/treatment_present/post (+2) or basal/post (+1)
param.levels{1} = numel(events.(opt.trEvents{1}).trial{1});
param.levels{2} = 2;

% For blocks, add the first and last trials to complete the 8 blocks
if events.(opt.trEvents{1}).trial{1}(1) == 0
    param.trial_change{1} = [1 events.(opt.trEvents{1}).trial{1}(2:end) length(neurons.itiOn{1})];
else
    param.trial_change{1} = [1 events.(opt.trEvents{1}).trial{1}(1:end) length(neurons.itiOn{1})];
end
param.trial_change{2} = events.(opt.trEvents{2}).trial{1}(1:end);

%% Whole trial raster and PSH, cluster by cluster
neuronSet = neurons.(toalignto{1}); % (HERE ONLY itiOn)
eventsSet = events.(toalignto{1});

if strcmp(param.trial2plot, 'correct') 
    eventsSet.code(conditions.aborted | conditions.incorrect | conditions.omission) = {[]};
elseif strcmp(param.trial2plot, 'incorrect')
    eventsSet.code(conditions.aborted | conditions.correct | conditions.omission) = {[]};
elseif strcmp(param.trial2plot, 'omission')
    eventsSet.code(conditions.aborted | conditions.incorrect | conditions.correct) = {[]};
elseif strcmp(param.trial2plot, 'allInitiated')
    eventsSet.code(logical(conditions.aborted)) = {[]};
end
trialrange = 1:length(conditions.correct);

% For each cluster
for c = 1:length(neuronSet)
    % cluster label and waveform
    param.label = spike.label{c};
    spikewf = mean(spike.waveform{1,c},2);

    % Do this for each requested situation (#levels)
    for p = 1:length(param.levels)
        if param.levels{p} > 2
            param.cond = 'tr1';
            param.plotcol = [0 0 1; 1 0 0; 1 0 1];
            param.raster.title = [param.label, ' ', param.trial2plot];    
            if param.FS2plot, param.raster.title = [param.raster.title, ' FS&NS'];
            else,             param.raster.title = [param.raster.title, ' NS']; end

            % Empy plot array for this cluster
            toplot = cell(1, param.nBlocks);

            %% For each Block
            for i = 1:param.nBlocks
                % Take all trials for both treatments
                toplot{i} = neuronSet{c};
    
                % Find out block type
                blks{1,i} = blocks.(conditions.type){i};

                % Get the range trial borders
                blks{2,i} = [param.trial_change{p}(i) param.trial_change{p}(i+1)];
    
                % Here, we find trials for the current block range
                cndidx = ismember(trialrange, param.trial_change{1}(i):param.trial_change{1}(i+1));
                
                % And we keep ONLY those for Novel Stimuli (tr2)
                if ~param.FS2plot, FSidx = ismember(trialrange, events.tr2.trial{1});
                else, FSidx = ones(1,length(trialrange)); end
                cndidx = trialrange(FSidx & cndidx);
                
                % We find the complementary set of trials (non-block (and FS, if so))
                notcndidx = trialrange(setdiff(1:end,cndidx));
                
                % Remove the complementary to NS
                toplot{i}(notcndidx) = {[]};
                
                % Empty trials of no interest
                if strcmp(param.trial2plot, 'correct')
                    toplot{i}(conditions.aborted | conditions.omission | conditions.incorrect) = {[]};
                    % param.ext = '_Block_NS_Correct.png';
                    param.subt{1} = 'CORRECT trials';
                elseif strcmp(param.trial2plot, 'incorrect')
                    toplot{i}(conditions.aborted | conditions.omission | conditions.correct) = {[]};
                    % param.ext = '_Block_NS_Incorrect.png';
                    param.subt{1} = 'INCORRECT trials';
                elseif strcmp(param.trial2plot, 'omission')
                    toplot{i}(conditions.aborted | conditions.correct | conditions.incorrect) = {[]};
                    % param.ext = '_Block_NS_Omission.png';
                    param.subt{1} = 'Omission trials';
                elseif strcmp(param.trial2plot, 'allInitiated')
                    toplot{i}(logical(conditions.aborted)) = {[]};
                    % param.ext = '_Block_NS_Omission.png';
                    param.subt{1} = 'All valid trials';
                end
            end

            %% FIGURE
            plot_blocks(toplot, spikewf, blks, eventsSet, param, opt)

        elseif param.levels{p} == 2
            param.cond = 'tr2';
            param.plotcol = [0 0 1; 1 0 0];
            param.raster.title = [param.label, ' ', 'NS vs FS'];
            
            toplot = cell(1,2);

            %% Now, take all trials for both treatments
            toplot{1} = neuronSet{c};
            toplot{2} = neuronSet{c};

            % Here, we index the tr2 code (13) for NS
            cndidx = events.(param.cond).trial{1};
            
            % find the complementary set of trials, the FS
            notcndidx = trialrange(setdiff(1:end,cndidx));

            if ~isempty(cndidx)
                % Empty FS trials
                toplot{1}(notcndidx) = {[]};

                % Empty NS trials 
                toplot{2}(cndidx) = {[]};
                
                % Also, empty Aborted and omited trials
                if strcmp(param.trial2plot, 'correct')
                    toplot{1}(conditions.aborted | conditions.omission | conditions.incorrect) = {[]};
                    toplot{2}(conditions.aborted | conditions.omission | conditions.incorrect) = {[]};
                    param.subt{1} = 'CORRECT NS trials';
                    param.subt{2} = 'CORRECT FS trials';
                    % param.ext = '_Correct.png';
                elseif strcmp(param.trial2plot, 'incorrect')
                    toplot{1}(conditions.aborted | conditions.omission | conditions.correct) = {[]};
                    toplot{2}(conditions.aborted | conditions.omission | conditions.correct) = {[]};
                    param.subt{1} = 'INCORRECT NS trials';
                    param.subt{2} = 'INCORRECT FS trials';
                    % param.ext = '_Incorrect.png';
                elseif strcmp(param.trial2plot, 'omission')
                    toplot{1}(conditions.aborted | conditions.correct | conditions.incorrect) = {[]};
                    toplot{2}(conditions.aborted | conditions.correct | conditions.incorrect) = {[]};
                    param.subt{1} = 'Omission NS trials';
                    param.subt{2} = 'Omission FS trials';
                    % param.ext = '_Omission.png'; 
                elseif strcmp(param.trial2plot, 'allInitiated')
                    toplot{1}(logical(conditions.aborted)) = {[]};
                    toplot{2}(logical(conditions.aborted)) = {[]};
                    param.subt{1} = 'All valid NS trials';
                    param.subt{2} = 'All valid FS trials';
                    % param.ext = '_All valid.png'; 
                end
            end

            %% FIGURE
            plot_treatments(toplot, blks, eventsSet, param, opt)

            close all hidden
            close all force
        end
    end
end
end % Main function

function plot_blocks(toplot, spikewf, blks, eventsSet, param, opt)
    color = cell(1,param.nBlocks);
    upperY = nan(1,param.nBlocks);

    fig = figure('visible', param.visible); % switch visibility
    set(fig, 'Position', [0, 0, round(param.screen.height), round(param.screen.width)]); % Set fig size as screen

    % Trial-long Spike Raster
    subplot(2,1,1) 
    trialCounter = 1; % always plot from first trial

    for i = 1:param.nBlocks
        if strcmp(blks(1,i),'A'), color{i} = param.plotcol(1,:); end
        if strcmp(blks(1,i),'B'), color{i} = param.plotcol(2,:); end
        if strcmp(blks(1,i),'C'), color{i} = param.plotcol(3,:); end

        % raster
        plotRaster(toplot{i}, ... % spikes
                    trialCounter, ... % trialCounter
                   'plotcol',    color{i}, ... % color per align (for now)
                   'spkwidth',   1,       ...
                   'linelength', param.lineLength, ...
                   'plotstyle',  param.plotStyle);
        
        yline(blks{2,i}(2)+1, '-', blks{1,i}, 'Color', color{i}, 'LabelVerticalAlignment', 'bottom',  'LineWidth', 1);
    end

    % Plot events
    plot_events(param.plotevent, eventsSet);

    % Prettify
    prettify(param.raster);
        xlim([-1500 param.timelim(2)])
        ylim([0 blks{2,i}(2)+1])
    xline(0,'--k')
    subtitle(param.subt{1})

    % PSH
    subplot(2,1,2)
    for i = 1:param.nBlocks
        upperY(1,i) = plotPSTH(toplot{i}, ... % spikes
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
        xlim([10 param.psh.xtick(end)])

    % inset waveform 
    axes('Position',[.79 .34 .1 .1])
        plot(spikewf, 'Color', 'k', 'LineWidth', 2), box off
        xlabel('ms'), ylabel('Amplitude (a.u.)')
        xlim([-2 96])
        xticks([0 32 64 96]); xticklabels({'-1' '0' '1' '2'});
        yticks([]); yticklabels({});

    % Save figure per alignment&cluster    
    sgtitle(param.raster.title)
    exportgraphics(fig, fullfile(opt.analysis,'newplots', ...
                        ['fulltrial_', param.raster.title, '.png']), ...
                        'Resolution', param.Resolution);
end

function plot_treatments(toplot, blks, eventsSet, param, opt)
    fig = figure('visible', param.visible); % switch visibility
    set(fig, 'Position', [0, 0, round(param.screen.width), round(param.screen.height)]); % Set fig size as screen

    % Trial-long Spike Raster
    for j = 1:2
        subplot(2,2,j) % NS
        trialCounter = 1; % always plot from first trial
        % raster
        plotRaster(toplot{j}, ... % spikes
                    trialCounter,               ... % trialCounter
                   'plotcol',    param.plotcol(j,:), ... % color per align (for now)
                   'spkwidth',   param.spkWidth,       ...
                   'linelength', param.lineLength,     ...
                   'plotstyle',  param.plotStyle);
    
        lasttr = cellfun(@isempty, toplot{j});
        lasttr = find(~lasttr,1,"last");
    
        % Plot events
        plot_events(param.plotevent, eventsSet);
    
        % Prettify
        prettify(param.raster);
            xlim([-1500 param.timelim(2)])
            ylim([0 lasttr+1])
        subtitle(param.subt{j})
        hold on

        for i = 1:param.nBlocks
            rectangle('position',[-2000 blks{2,i}(1)+1 ...
                                  12000 blks{2,i}(2)-blks{2,i}(1)+1], ...
                      "LineWidth", 1, ...
                      "EdgeColor", [.5 .5 .5]);
        end
    end

    % Trial-long PSH
    subplot(2,2,[3 4])
    for j = 1:2
        upperY(j) = plotPSTH(toplot{j}, ... % spikes
                param.stepSz,   ... % stepSz
                param.binSize,  ...  % binSize
                param.timelim,  ... % interval
                param.smpRate,  ...  % samples per second in the feeded data
                'plotcol',      param.plotcol(j,:),...
                'meanline',     '-',...
                'smoothplot',   true);  
        hold on
    end

    % Prettify
    prettify(param.psh)
        xline(param.psh.xtick(find(param.psh.xticklabels{1}{1} == 0)),'--k');
        maxY = max(upperY); if maxY <= 5, maxY = 5; end % force a minimum y-axis scale
        ylim([0 maxY*1.1]);
        xlim([10 param.psh.xtick(end)])
    legend({'', 'Novel', '', 'Familiar'})

    % Save figure per alignment&cluster
    sgtitle(param.raster.title)
    exportgraphics(fig,fullfile(opt.analysis,'newplots', ...
                        ['fulltrial_', param.raster.title, ' ', param.trial2plot, '.png']), ...
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
    
                line([evtime2 evtime2]-1000, [i i+1], 'Color', color, 'LineWidth', 1) % Fixing behavior artificial delay of 1 second (to Stim2!)
            else
                evtime = events.time{i,1}(evnt)*1000;
                line([evtime evtime], [i i+1], 'Color', color, 'LineWidth', 1)
            end
        end
    end
end
end