%% NGL03_plotting
% To run after all standard data have been sorted, curated and saved. In
% principle, different tyoes of plots could be selected to be done or not,
% and new plots could be added for personalization. Parameters are taken
% from main script or from an additional one.
%
% Jesus 04.07.2024

%% 00. Check current inputs.
% Check if input variable exist already. Parse values.
if ~exist("input","var")
    input = struct( 'datadrive' , datadrive , ...   % force char array
                    'studyName' , studyname , ...   % force char array
                    'toolbox'   , toolbox   , ...   % force char array
                    'subjects'  , [], ...           % do NOT force char array
                    'dates'     , []        );      % do NOT force char array
    input.dates     = dates;    % place as it comes
    input.subjects  = subjects; % place as it comes
else
    disp('Using INPUTS from NGL01_MAIN.')
end

% Set default inputs and dependencies. In case NGL01 did not before.
input = set_default(input, opt);
% skip =1 ;

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);
dosave = 1;

for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        input.run = [x y]; % Current run, to pass to functions.
        
        %% 02. Proceed only once
        if exist(fullfile(input.analysis, "data_all.mat"),"file")
            load(fullfile(input.analysis, "data_all.mat"));
            dosave = 0;
            break
        end
        
        %% 03. Prepare to proceed with a single session.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);           
        if ~exist(fullfile(opt.analysis,'newplots'),"dir")
            mkdir(fullfile(opt.analysis,'newplots'))
        end

        %% 04. Recover all Project data if not collected yet
        % Recover trial definitions created after event extraction and processing. 
        % Will have as many variations as requested at that time. Needs to be ran 
        % again to create new alignments.
        if ~exist('events','var'),    load(fullfile(opt.analysis, "events.mat")),      end
        if ~exist('neurons','var'),   load(fullfile(opt.analysis, "neurons.mat")),     end
        if ~exist('condition','var'), load(fullfile(opt.analysis, "condition.mat")),   end 
        if ~exist('spike','var'),     load(fullfile(opt.spikeSorted, "spike.mat")),    end
        if ~exist('blob','var')     && opt.plot_SocLear
            load(fullfile(opt.analysis, "blob.mat")),
        end
        
        % Collect all data into single all variables: events, neurons, conditions, blob
        allneurons{x,y}     = neurons;
        allevents{x,y}      = events;
        allconditions{x,y}  = conditions;
        allspike{x,y}       = spike;
        if opt.plot_SocLear
           allblobs{x,y}    = blob;
        end
    
        clear events neurons spike blob spike
    end
end

%% 02. Save it for next iterations
if dosave && opt.plot_SocLear
    save(fullfile(input.analysis, "data_all.mat"), "allblobs","allspike","allconditions","allevents","allneurons")
elseif dosave && ~opt.plot_SocLear 
    save(fullfile(input.analysis, "data_all.mat"), "allspike","allconditions","allevents","allneurons")
end

%% 03. Plot session by session
for x = 1:input.nsubjects % Subjects.   
    for y = 1:input.sessions(x).nsessions % Sessions.
        input.run = [x y]; % Current run, to pass to functions.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);           

        %% 3.1 All clusters piled, ...
        % aligned to requested events, for all clusters
        if opt.pooledstats
            plot_pooledstats(allneurons{x,y}, opt, param, input)
        end
                   
        %% 3.2 Raster, trace, PSH, ISI
        % whole trial
        if opt.plot_trial
            if ~exist('spike','var'),     load(fullfile(opt.spikeSorted, "spike.mat")),    end
            plot_trialclusters(allneurons{x,y}, allevents{x,y}, spike, allconditions{x,y}, opt, param)
        end
        
        % aligned to requested events, per cluster
        if opt.plot_align && numel(opt.alignto) > 1
            plot_alignedclusters(allneurons{x,y}, allevents{x,y}, spike, allconditions{x,y}, opt, param)
        end
        
        %% 3.3 SocialLearning-specific plots
        if opt.plot_SocLear
            plot_SocLear(allneurons{x,y}, allevents{x,y}, allspike{x,y}, allconditions{x,y}, allblobs{x,y}, opt, input, param)
        end
    end
end

%% 04. Across sessions plots
if opt.plot_sessions
    %% 4.1 Distances and velocities
    if ~skip
    % Title
    param.title = 'Subject-Feeder and Subject-Tutor distance';
    for x = 1:input.nsubjects % Subjects.
        param.subtitle = ['Subject:', input.subjects(x).name];
        fig = figure('visible', param.visible); % switch visibility
        set(fig, 'Position', [0, 0, 1080, 1920]); % Set fig size as screen

        for y = 1:input.sessions(x).nsessions % Sessions.
            input.run = [x y]; % Current run
            
            path = fullfile(input.analysis, input.subjects(input.run(1)).name, input.sessions(input.run(1)).list{input.run(2)});
            if exist(fullfile(path, 'Dist2Feed_Vel.mat'),"file")
                load(fullfile(path, 'Dist2Feed_Vel.mat'));
            else, continue
            end
    
            % Get data
            eucl     = allblobs{x,y}.euclDist(1:30:end);
            merges   = allblobs{x,y}.Merges;
            mindist     = Dist2feeder(1:30:end,:);
            vel      = AgentVel; % sampled by seconds
            % and times
            sesstime = allblobs{x,y}.Time(1:30:end); 
            veltime  = allblobs{x,y}.Time(1:60:end);
    
            % Check data sizes
            if length(veltime)>length(vel), veltime = veltime(1:length(vel)); 
            else,                                vel = vel(1:length(veltime)); end
    
            if length(eucl)> length(sesstime), eucl = eucl(1:length(sesstime));
            else,                              sesstime = sesstime(1:length(eucl)); end
            
            if length(mindist)> length(sesstime), mindist = mindist(1:length(sesstime),:);
            else,                              sesstime = sesstime(1:length(mindist)); end
      
            % index session solo time
            sesssoloidx = sesstime<5*60 | sesstime>16*60;
            mindist(~sesssoloidx,:) = nan; 
    
            % Plot time series of euclidean distance between blobs and between blob and feeders
            subp = subplot(input.sessions(x).nsessions,1, y);
                pos = subp.Position;
                
                % Distance between blobs(at tutor time)
                plot(subp, sesstime, eucl, 'LineStyle','-','Color',[.25 .25 .25],'LineWidth', 1)
                    box off, hold on
    
                    % Add merge periods as black transparent patches
                    for r = 1:size(merges,1)
                        rectangle('Position', [merges(r,1) 0 merges(r,2)-merges(r,1) 500], ...
                            'FaceColor', [0, 0, 0, 0.2], 'EdgeColor', 'none');
                    end
                    xline([5*60,16*60], 'LineStyle', ':','Color', 'k', 'LineWidth', 1)
                    
                    % Remove solo times, i.e show only 250 < t < 1000. Add labels.
                    xlim([1 sesstime(end)]); xlabel('min');
                    xticks(0:120:sesstime(end)); xticklabels(string(xticks/60))
        
                    ylim([-1 500]); ylabel('Eucl. dist. (px)');  
                    yticks(0:100:400);
    
                % Distance from subject to Feeders
                plot(subp, sesstime, mindist, 'LineStyle','-', 'LineWidth', 1)
            
                % Add the agent velocity to the top, axis on the left at proper scale
                ax = axes('Position', [pos(1), pos(2)+0.84*pos(4), pos(3), 0.025]); 
                    plot(ax, veltime, AgentVel,  'LineStyle','-','Color', 'r', 'LineWidth', 0.5)
                        ax.XAxis.Visible = 'off';
                        ax.YAxis.Color = 'r';
                        set(ax, 'YAxisLocation', 'left');
                        xlim([0 veltime(end-1)]);
                        ylim([0 max(AgentVel)*1.1]); ylabel('px/s');  
                hold off
    
        end

        % Set title and save
        sgtitle([{param.title}, {param.subtitle}]);
        savepath = fullfile(input.analysis, input.subjects(input.run(1)).name);
        if ~exist(fullfile(savepath, 'plots'),"dir"), mkdir(fullfile(savepath, 'plots')); end
        exportgraphics(fig, fullfile(savepath, 'plots', ...
                        ['DistVel_', input.subjects(input.run(1)).name, '.png']), ...
                        'Resolution', param.Resolution);
        close all
    end
    %% 4.2 Distance to Feeder aligned to StimOn 
    % Title
    param.title = 'Distance to Feeder, aligned';
    plotcol = lines(4); % create colormap with as 4 levels
    for x = 1:input.nsubjects % Subjects.
        param.subtitle = ['Subject:', input.subjects(x).name];
        fig = figure('visible', param.visible); % switch visibility
        set(fig, 'Position', [0, 0, 1080, 1920]); % Set fig size as screen
        tlength = 6; % seconds to extract around stimOn
        sub = 0;
       
        for y = 1:input.sessions(x).nsessions % Sessions.
            input.run = [x y]; % Current run
            
            path = fullfile(input.analysis, input.subjects(input.run(1)).name, input.sessions(input.run(1)).list{input.run(2)});
            if exist(fullfile(path, 'Dist2Feed_Vel.mat'), "file")
                load(fullfile(path, 'Dist2Feed_Vel.mat'));
            else, continue
            end

            path = fullfile(input.trialSorted, input.subjects(input.run(1)).name, input.sessions(input.run(1)).list{input.run(2)});
            if exist(fullfile(path, 'trialdef.mat'), "file")
                load(fullfile(path, 'trialdef.mat'));
            else, continue
            end

            % and times
            sesstime = allblobs{x,y}.Time; % sec
            trialdef = trialdef{2,2}/1000; % sec

            % Get data
            ntrials = size(trialdef,1);
            dist = Dist2feeder;
            for i = 1:ntrials
                for ii = 1:4
                    distTrial{i}(:,ii) = dist(sesstime > trialdef(i,3)-1 & sesstime < trialdef(i,3)+(tlength-1),ii);
                    if length(distTrial{i}(:,ii)) > tlength*60
                        distTrial{i} = distTrial{i}(1:tlength*60,ii);
                    end
                end
            end
            clear dist

            for ii = 1:4
                sub = sub+1;
                subplot(input.sessions(x).nsessions,4,sub); hold on
                for t = 1:ntrials
                    if t > 2 
                        if length(distTrial{t}(:,ii)) < length(dist(:,t-1))
                            distTrial{t}(length(dist(:,t-1)),:) = nan;
                        end
                    end
                    dist(:,t) = distTrial{t}(:,ii);
                    plot(1:tlength*60, dist(:,t), 'Color', plotcol(ii,:), 'LineStyle', '-', 'LineWidth', 0.2);
                end
                plot(1:tlength*60, mean(dist(:,t), 2, "omitnan"), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
                ylabel('Eucl. dist (px)');  xlabel('time (s)'); 
                ylim([0 125]); yticks([0:25:125]);
                xlim([0 tlength*60]); xticks(0:60:tlength*60); xticklabels({'-1' '0' '1' '2' '3' '4' '5'});
                xline(60, LineStyle="--", LineWidth=2);
                hold off
            end
        clear distTrial

        end
        
        % Set title and save
        sgtitle([{param.title}, {param.subtitle}]);
        savepath = fullfile(input.analysis, input.subjects(input.run(1)).name);
        if ~exist(fullfile(savepath, 'plots'),"dir"), mkdir(fullfile(savepath, 'plots')); end
        exportgraphics(fig, fullfile(savepath, 'plots', ...
                        ['DistVel_', input.subjects(input.run(1)).name, '.png']), ...
                        'Resolution', param.Resolution);
        close all
    end
    end
    %% 4.2. Position heatmaps 
     for x = 1:input.nsubjects % Subjects.
        for y = 1:input.sessions(x).nsessions % Sessions.
            input.run = [x y]; % Current run
            
            C = allblobs{x,y}.Pos;
            fun = @(v) v(1,:);
            posDatall = {};

            % Get whole session time for position
            sesstime = allblobs{x,y}.Time'; 

            % define times
            preidx = sesstime < 5*60;
            postidx = sesstime > 16*60;
            sesssoloidx = sesstime < 5*60 | sesstime > 16*60;

            % define non-empty time bins
            notemptyidx = ~cellfun('isempty',C);

            % Get blob detection indexed by non-empty and 'solo' time
            posData_pre = cellfun(fun, C(preidx & notemptyidx), 'UniformOutput',false);
            posData_post = cellfun(fun, C(postidx & notemptyidx), 'UniformOutput',false);
            posDatall = cellfun(fun, C(sesssoloidx & notemptyidx), 'UniformOutput',false);

            % Convert to continuous matrix
            posData_pre = cell2mat([posData_pre(:)]);
            posData_post = cell2mat([posData_post(:)]);
            posDatall = cell2mat([posDatall(:)]);
 
            % Plot Parameters
            heatmappar.title    = 'Position density';
            heatmappar.subtitle = [input.subjects(x).name, ' ', input.sessions(x).list(y)];
            heatmappar.NumBins  = [76 76];
            heatmappar.BinWidth = [10 10];
            heatmappar.visible  = 'off';
            heatmappar.AXvisible= 'off';
            heatmappar.shape    =  [473 438; 533 433; 616 431; 626 347; 640 292; 729 293; 825 304;  % Vertex of the polygon fitting the 
                                    832 361; 830 436; 719 547; 592 645; 517 636; 462 622; 462 526] ...
                                    - [410 240]; % In this experiment, we removed one half
            heatmappar.flipY    = 1;
            heatmappar.norm     = 'probability';
            heatmappar.colmap   = 'hot';
            heatmappar.savepath = fullfile(input.analysis, input.subjects(input.run(1)).name, 'plots');

            % Plot and save
            PlotBinSeq([], posDatall, heatmappar);

            % Plot and save
            heatmappar.title = 'Position density (PreTut)';
            PlotBinSeq([], posData_pre, heatmappar);

            % Plot and save
            heatmappar.title = 'Position density (PostTut)';
            PlotBinSeq([], posData_post, heatmappar);
        end
     end

    %% 4.3. Position/Spike heatmaps
    for x = 1:input.nsubjects % Subjects.
        for y = 1:input.sessions(x).nsessions % Sessions.
            input.run = [x y]; % Current run
            
            C = allblobs{x,y}.Pos;
            fun = @(v) v(1,:);
            posDatall = {};

            % Get whole session time for position
            sesstime = allblobs{x,y}.Time'; 

            % define times
            preidx = sesstime < 5*60;
            postidx = sesstime > 16*60;
            sesssoloidx = sesstime < 5*60 | sesstime > 16*60;

            % define non-empty time bins
            notemptyidx = ~cellfun('isempty',C);

            % Get blob detection indexed by non-empty and 'solo' time
            posData_pre = cellfun(fun, C(preidx & notemptyidx), 'UniformOutput',false);
            posData_post = cellfun(fun, C(postidx & notemptyidx), 'UniformOutput',false);
            posDatall = cellfun(fun, C(sesssoloidx & notemptyidx), 'UniformOutput',false);

            % Convert to continuous matrix
            posDat{1} = cell2mat([posData_pre(:)]);
            posDat{2} = cell2mat([posData_post(:)]);
            posDat{3} = cell2mat([posDatall(:)]);

            % reduce session time to same idexing, for later
            ts_t{1} = sesstime(preidx & notemptyidx); 
            ts_t{2} = sesstime(postidx & notemptyidx); 
            ts_t{3} = sesstime(sesssoloidx & notemptyidx); 

           % Plot parameters
           frparam.title = 'FR_position_';  %
           frparam.subt{1} = '_bsl';
           frparam.subt{2} = '_post';
           frparam.nbins   = 50;        % Define grid resolution (higher = finer)
           frparam.sr      = 60;     % Sampling rate in Hz (samples per second)
           frparam.sig     = 1;     % Gaussian kernel width (second)
           frparam.shape   = [473 438; 533 433; 616 431; 626 347; 640 292; 729 293; 825 304;  % Vertex of the polygon fitting the 
                              832 361; 830 436; 719 547; 592 645; 517 636; 462 622; 462 526] ...
                              - [410 240];
           frparam.visible   = 'off';      % see the figure
           frparam.AXvisible = 'off';      % hide axes
           frparam.colmap    = 'hot';      % colormap
           frparam.flipY     = 1;          % Reverse Y axis
           frparam.barunit   = 'sp/s';
           frparam.savepath  = fullfile(input.analysis, 'plots', input.subjects(input.run(1)).name, input.sessions(x).list{y});

            % go per cluster
            for c = 1:size(allspike{x,y}.timestamp,2)
               tstamps{1} = [];
               tstamps{2} = [];
               tstamps{3} = [];

               % From all spikes, use only those during solo time
               D = allspike{x,y}.timestamp{1,c};
               
               preidx = D < 5*60;
               postidx = D > 16*60;
               sesssoloidx = D < 5*60 | D > 16*60;
               
               tstamps{1} = D(preidx);
               tstamps{2} = D(postidx);
               tstamps{3} = D(sesssoloidx);
               clear D

               for i = 1:2
                  % Plot parameters
                  frparam.subtitle  = ['c', allspike{x,y}.label{1,c}, frparam.subt{i}];
                  
                  % Plot
                  plot_fr_position(tstamps{i}, ts_t{i}, posDat{i}, frparam)
               end
                
            end
        end
    end

end