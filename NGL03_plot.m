%% NGL03_plot Script 
% To create some preliminary figures with the info extracted from NGL02.
% In the works. Provisional until we have the events extracted with the
% neural data.
% Do not clear the variables 'input', 'results', 'opt' from previous step.

%% Options
% For now, testing in Pilot_SocialLearning: load '*par.mat'
input.bhvfolder = fullfile(input.datadrive, input.studyName, '\data\behavior\');

% To plot or not to plot, that's ...
opt.plotrasters = true;

% Set EventCode meaning
itiOn       = 1;
stimOn      = 2;
rwd         = 5;
tutor       = 4;

% Set Event to align rasters t0 to.
alignto = rwd;

%% 01. From all subjects and sessions, compile some data
% Prepare fields
dat.nclusters   = [];
dat.clusID      = cell(size(results));
dat.clusCAT     = cell(size(results));
dat.amps        = cell(size(results));
dat.trials      = cell(size(results));
dat.spikes      = cell(size(results));
dat.spikesID    = cell(size(results));

% Extract data from 'results' to 'dat'
for a = 1:input.nsubjects % animals
    for ss = 1:size(results,1) % sessions
        if ~isempty(results{ss,a})
            dat.info(1,a)        = {input.subjects(a).name};        % animal code
            dat.info(2,a)        = {sessions(a).list'};             % session dates
            dat.nclusters(ss,a)  = length(results{ss,a}.spike.cgs); % #clusters per animal/session
            dat.clusID(ss,a)     = {results{ss,a}.spike.cids};      % Cluster's IDs
            dat.clusCAT(ss,a)    = {results{ss,a}.spike.cgs};       % Clsuter's category (good, mua, noise)
            dat.spikes(ss,a)     = {results{ss,a}.spike.st};        % Spike times
            dat.spikesID(ss,a)   = {results{ss,a}.spike.clu};       % Spike's cluster IDs 
            dat.amps(ss,a)       = {results{ss,a}.spike.Amps};      % Spike's amplitude

        end
    end
end

clear a ss clus

%% 02. Plot General stuff
if opt.plotgeneral
    % Clusters in total
    X = categorical({input.subjects.name});
    figure,
    b = bar(X,[sum(dat.nclusters(:,1)); sum(dat.nclusters(:,2))]);
        ylabel('total clusters')
        ax = gca; ax.FontSize = 14;
        box("off")
        sdf('bold_600')
    
    % Clusters per sessions
    figure,
    plot(dat.nclusters);
        xlim([0 size(dat.nclusters,1)+1]), ylim([0 max(max(dat.nclusters))+1]);
        xlabel('Session'), ylabel('n clusters')
        ax = gca; ax.FontSize = 14;
        hold on
        line([1 size(dat.nclusters,1)],[nanmean(dat.nclusters(:,1)) nanmean(dat.nclusters(:,1))],'LineStyle',':','Color','b')
        line([1 size(dat.nclusters,1)],[nanmean(dat.nclusters(:,2)) nanmean(dat.nclusters(:,2))],'LineStyle','-.','Color','r')
        legend({input.subjects.name}); legend(Box="off")
        hold off, box("off")
        sdf('bold_600')
    
    clear ax b X
end

%% 03. Trial parse.
% Loop over subjects and session to collect and organize timestamps and spikes
for s = 1:input.nsubjects
    for ss = 1:sessions(s).nsessions

       % 01. Collect eventcodes.
       file = ls([fullfile(input.bhvfolder,input.subjects(s).name,sessions(s).list{ss}), '\*par.mat']); % find session
       evnt = load(fullfile(input.bhvfolder,input.subjects(s).name,sessions(s).list{ss}, file),...
                    'SaveEvnts'); % collect events in temp variable
       events{ss,s} = evnt.SaveEvnts; % place in final variable
       clear evnt file

       % Convert spike times to seconds.
       dat.spikes{ss,s}    = seconds(dat.spikes{ss,s});

       % Create a variable for trial start times from the behavioral results
       % initialize trial matrix with length as spikes and 2 columns (trial, tutor)
       dat.trials{ss,s}  = zeros(size(dat.spikes{ss,s}));
       dat.tutor{ss,s}   = zeros(2,1);
       dat.tutorID{ss,s} = zeros(size(dat.spikes{ss,s}));

       % Set events to codes
       dat.events.itiOn{ss,s} = events{ss,s}(events{ss,s}(:,2)==itiOn); % Get all timestamps for t zero
       dat.events.tutor{ss,s} = events{ss,s}(events{ss,s}(:,2)==tutor); % Get timestamps for tutor in & out
       dat.events.stimOn{ss,s} = events{ss,s}(events{ss,s}(:,2)==stimOn); % Get timestamps for stimOn
       dat.events.rwd{ss,s} = events{ss,s}(events{ss,s}(:,2)==rwd); % Get timestamps for stimOn

       % Convert times relative to session start and to seconds
       dat.events.itiOn{ss,s} = dat.events.itiOn{ss,s}-events{ss,s}(1,1); % use first event as session start time
       dat.events.itiOn{ss,s} = seconds(dat.events.itiOn{ss,s});          % Convert values to seconds

       dat.events.tutor{ss,s} = dat.events.tutor{ss,s}-events{ss,s}(1,1); % use first event as session start time
       dat.events.tutor{ss,s} = seconds(dat.events.tutor{ss,s});          % Convert values to seconds

       dat.events.stimOn{ss,s} = dat.events.stimOn{ss,s}-events{ss,s}(1,1); % use first event as session start time
       dat.events.stimOn{ss,s} = seconds(dat.events.stimOn{ss,s});          % Convert values to seconds

       dat.events.rwd{ss,s} = dat.events.rwd{ss,s}-events{ss,s}(1,1); % use first event as session start time
       dat.events.rwd{ss,s} = seconds(dat.events.rwd{ss,s});          % Convert values to seconds

       % Identify which spike times belong to which trial.
       tr = 1; % Initialize at trial 1
       tutorin = 0; % Initially Tutor is NOT present
       
       % Evaluate trial for every spike time respect to a certain time
       % when we align zero.
%        if alignto == 1 %bc some sessions don't have stimOn
           dat.t2align{ss,s} = dat.events.itiOn{ss,s};
%        elseif alignto == 2 % then we use itiOn
%            dat.t2align{ss,s} = dat.events.stimOn{ss,s};
%        elseif alignto == 5 
%            dat.t2align{ss,s} = dat.events.rwd{ss,s};
%        end

       for j = 1:length(dat.spikes{ss,s}) 
            % Spike times before first trial, tag as 0
            if tr == 1 && dat.spikes{ss,s}(j) < dat.t2align{ss,s}(tr) 
                dat.trials{ss,s}(j) = 0;

            % Spike times between current and next trial starts, tag as current trial
            elseif (dat.spikes{ss,s}(j) >= dat.t2align{ss,s}(tr)) && (dat.spikes{ss,s}(j) < dat.t2align{ss,s}(tr+1))
                dat.trials{ss,s}(j) = tr;

                % Take trial of tutor event and set flag
                if ~isempty(dat.events.tutor{ss,s}) % If tutor events exist at all
                   if (dat.spikes{ss,s}(j) >= dat.events.tutor{ss,s}(1)) && (dat.spikes{ss,s}(j) < dat.events.tutor{ss,s}(2))
                      dat.tutorID{ss,s}(j) = 1; % tag spike as 'with tutor'
                      
                      if dat.tutor{ss,s}(1)==0
                         dat.tutor{ss,s}(1) = tr; % Set entry trial
                      end
                   end
                end

            % Once spike time is equal or greater than next start trial, 
            % increment trial and tag this first spike time as belonging to it
            % Next spike time should fall in the above condition again.
            else 
                tr = tr+1;
                dat.trials{ss,s}(j) = tr;

                % Evaluate if tutor is still in at tr+1
                if ~isempty(dat.events.tutor{ss,s}) % If tutor events exist at all
                   if (dat.spikes{ss,s}(j)>dat.events.tutor{ss,s}(2)) && (dat.tutor{ss,s}(1)>0)
                       dat.tutor{ss,s}(2) = tr; % Trial when tutor is gone
                       % dat.tutor{ss,s}(2) = tr-1; % Or last trial tutor was in?
                   end
                end

            end
        
            % Once we reach the last trial, we can not evaluate to the next one.
            % All remaining spike times are tagged as last trial. (not optimal)
            if tr == length(dat.t2align{ss,s})-1
                dat.trials{ss,s}(j:end, 1) = length(dat.t2align{ss,s});
                break % and we break the for loop.
                      % until fixed
            end
       end
       clear j
       
%        % Remove all assigned to trial 0. (Do not belong to the session)
%        dat.spikes{ss,s}(dat.trials{ss,s}==0)   = [];
%        dat.spikesID{ss,s}(dat.trials{ss,s}==0) = [];
%        dat.tutorID{ss,s}(dat.trials{ss,s}==0)  = [];
% 
%        % And remove themselves
%        dat.trials{ss,s}(dat.trials{ss,s}==0)   = [];
    end
end
    
%% 04. Raster plots.
% Will plot a raster plot per cluster, or a figure with all clusters.
if opt.plotrasters
   % add rasters toolbox
   addpath('C:\Code\spike-raster-plot')

   % Loop subjects, sessions and clusterIDs to plot their rasters
    for s = 1:3%input.nsubjects
        for ss = 1:sessions(s).nsessions
           % prepare output file
           savefig = [fullfile(input.sorted,input.subjects(s).name,sessions(s).list{ss}), '\spikeID_'];
           
           % To plot clusters separatedly, we loop the plotting function.
           i = 1; % Initilize to plot handle s.plot(1)
           for clus = dat.clusID{ss,s}
                figure;
                p.plot{i} = spikeRasterPlot(dat.spikes{ss,s}(dat.spikesID{ss,s}==clus), ...
                                            dat.trials{ss,s}(dat.spikesID{ss,s}==clus));

                    p.plot{i}.AlignmentTimes = dat.events.itiOn{ss,s} ;% dat.t2align{ss,s}; % Align all trials to the stablished zero
                    p.plot{i}.GroupData      = dat.tutorID{ss,s}(dat.spikesID{ss,s}==clus);
                    p.plot{i}.LegendTitle    = 'Tutor';
                    p.plot{i}.TitleText      = ['Cluster ',int2str(clus)];
                    p.plot{i}.YLabelText     = 'Trial';
                    if alignto == 1
                        p.plot{i}.XLimits        = seconds([-1 12]);    % Show from -1 to +10 seconds
                        p.plot{i}.XLabelText     = 'Time since itiOn (s)';
                    elseif alignto == 2
                        p.plot{i}.XLimits        = seconds([-1 10]);    % Show from -1 to +35 seconds
                        p.plot{i}.XLabelText     = 'Time since stimOn (s)';
                    elseif alignto == 5
                        p.plot{i}.XLimits        = seconds([-7 8]);    % Show from -1 to +35 seconds
                        p.plot{i}.XLabelText     = 'Time since rwd (s)';
                    end

                % Saving
                saveas(gcf,[savefig,int2str(clus),'.png']);
    
                close(gcf)
                i = i+1;
           end 
    
           % All clusters can be plotted together, with different colors by
           % running the commented piece below.
    %             p = spikeRasterPlot(dat.spikes{ss,s}, dat.trials{ss,s});
    %                 p.AlignmentTimes = dat.events.itiOn{ss,s};
    %                 p.GroupData      = dat.spikesID{ss,s};  
    %                 p.XLimits        = seconds([-1 25]);
        end
    end
end

%% Single cluster firing rates TODO
% if opt.plotfrs
% 
% end