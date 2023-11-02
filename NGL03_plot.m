%% NGL03_plot Script 
% to create some preliminary figures with the info extracted from NGL02

% Do not clear the variables 'input' and 'results'

% For rasters set here which animal and session you want.
% General plots do not need this.
animal  = 1;
session = 6;

%% 01. Open behavior file with collected eventcodes.
% Needed for rasters and any single cluster plot
% For now, testing in Pilot_SocialLearning: load '*res.mat'
[file,path] = uigetfile('*.mat');

if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
   behavior = load(fullfile(path,file), 'res');
   behavior = behavior.res;
end
clear path file

%% 02. From all subjects and sessions, compile some data
dat.nclusters = zeros(size(results));
dat.clusID = cell(size(results));
dat.clusCAT = cell(size(results));
dat.times = cell(size(results));
dat.timesID = cell(size(results));
dat.amps = cell(size(results));

for a = 1:size(results,2) % animals
    for t = 1:size(results,1) % sessions
        if ~isempty(results{t,a}.spike.spikeTemplates)
            dat.nclusters(t,a)  = length(results{t,a}.spike.cgs);
            dat.clusID(t,a)     = {results{t,a}.spike.cids};
            dat.clusCAT(t,a)    = {results{t,a}.spike.cgs};
            dat.times(t,a)      = {results{t,a}.spike.st};
            dat.timesID(t,a)    = {results{t,a}.spike.clu};
            dat.amps(t,a)       = {results{t,a}.spike.Amps};
        end
    end
end
clear a t

%% Plot General stuff
% Clusters in total
    X = categorical({input.subjects.name});
    figure,
    b = bar(X,[sum(dat.nclusters(:,1)); sum(dat.nclusters(:,2))]);
    ylabel('total clusters')
    ax = gca; ax.FontSize = 14;

% Clusters per sessions
    figure,
    plot(dat.nclusters);
    xlim([0 size(dat.nclusters,1)+1]), ylim([0 max(max(dat.nclusters))+1]);
    xlabel('Session'), ylabel('n clusters')
    ax = gca; ax.FontSize = 14;
    hold on
    line([1 size(dat.nclusters,1)],[mean(dat.nclusters(:,1)) mean(dat.nclusters(:,1))],'LineStyle',':','Color','b')
    line([1 size(dat.nclusters,1)],[mean(dat.nclusters(:,2)) mean(dat.nclusters(:,2))],'LineStyle','-.','Color','r')
    legend({input.subjects.name}); legend(Box="off")
    hold off

% How many 'good' and 'mua'
    for t = 1:length(dat.clusCAT)
        good(t,:) = [sum(dat.clusCAT{t,1}==1) sum(dat.clusCAT{t,2}==1)];
        mua(t,:)  = [sum(dat.clusCAT{t,1}==2) sum(dat.clusCAT{t,2}==2)];
    end
    
    figure,
    subplot(121)
    plot(good);
    xlim([0 size(dat.nclusters,1)+1]), ylim([0 max(max(dat.nclusters))+1]);
    xlabel('Session'), ylabel('n GOOD clusters')
    ax = gca; ax.FontSize = 14;
    
    subplot(122)
    plot(mua);
    xlim([0 size(dat.nclusters,1)+1]), ylim([0 max(max(dat.nclusters))+1]);
    xlabel('Session'), ylabel('n MUA clusters')
    ax = gca; ax.FontSize = 14;
    legend({input.subjects.name}); legend(Box="off")

    clear t ax b X s ss

%% Plot rasters
% add rasters toolbox
addpath('C:\Code\spike-raster-plot')

    % Extract cluster ID, spike times and spike IDs from one session
    % and make easy to access
    clusIDset = dat.clusID{session, animal};
    spikes    = seconds(dat.times{session, animal});
    spikesID  = dat.timesID{session, animal};
    
    % Create a variable for trial start times from the behavioral results
    trials      = zeros(size(spikes));
    trialstarts = cell2mat({behavior(2:end).Timelapse}.');
    trialstarts = minutes(trialstarts);
    
    % Identify which spike times belong to which trial.
    t = 1; % Initialize trial 1
    for j = 1:length(spikes) % evaluate every spike time
        if t == 1 && spikes(j)<trialstarts(t) 
            % spikes before first trial, tag as 0
            trials(j) = 0; 
        elseif (spikes(j)>=trialstarts(t)) && (spikes(j)<trialstarts(t+1))
            % spikes between current and next trial starts, tag as current trial
            trials(j) = t; 
        else 
            % Once spike time is equal or greater than next start trial, 
            % increment trial and tag the first spike time as belonging to it
            % Next spike time should fall in the above condition again.
            t = t+1; 
            trials(j) = t;
        end
    
        % This is until we reach the last trial, where we can not evaluate
        % to the unexistent next one. In taht case, we tag all remaining
        % spike times as last trial. (Which in not totally right)
        if t == length(trialstarts)-1
            trials(j:end, 1) = length(trialstarts);
            break % and we break the for loop bc we do not need to go one by one.
                  % until fixed
        end
    end
    clear t j
    
    % Remove all assigned to trial 0. It messes up the plotting.
    spikes(trials==0)   = [];
    spikesID(trials==0) = [];
    trials(trials==0)   = [];
    
    % And proceed to plot. To plot clusters separatedly, we loop over the
    % plotting function.
    % All clusters can be plotted together, with different colors by
    % running the commented piece below.
    i = 1; % Initilize to plot handle s.plot(1)
    for t = clusIDset
        figure,
        s.plot{i} = spikeRasterPlot(spikes(spikesID==t), trials(spikesID==t));
            s.plot{i}.AlignmentTimes = trialstarts;
            s.plot{i}.XLimits        = seconds([-1 25]);
        i = i+1;
    end

%     % All clusters together (Uncomment if wanted)
%     s = spikeRasterPlot(spikes, trials);
%         s.AlignmentTimes = trialstarts;
%         s.GroupData      = spikesID;  
%         s.XLimits        = seconds([-1 25]);


%% Single cluster firing rates
% TODO
    