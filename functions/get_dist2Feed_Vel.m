function [Dist2feeder, AgentVel] = get_dist2Feed_Vel(blob, Fs, tutorTime, input)
% Distance to feeder and agent velocity calculations
%
% INPUT:    blob; struct containing all data from blob tracking process
%           Fs; integer with video sampling rate, to adjust velocity to px per second
%           tutorTime; 2x1 integer vector with time when tutor comes in and gets out (in minutes)

% Get blob config data for this project, mostly bc of feeder positions
feeder = [];
blobConfig; % 'blobConfig.m' file should exist under '.../analysisCode'
path = fullfile(input.analysis, input.subjects(input.run(1)).name, input.sessions(input.run(1)).list{input.run(2)});

% Calculate times in/out
timeOut = [tutorTime(1)*60*Fs tutorTime(2)*60*Fs]; % pre/post tutor time (t<5 & t>15 min, at 60Hz)

% Re-find last ID
lastID = size(blob.Tracking,3);

% allocate variable
Dist2feeder = nan(length(blob.Time),lastID,4);

%% Calculate distance to feeders.
% Calculate Eucl. dist from blob position to each fixed feeder position, in
% a frame by frame manner. Use the process to check for 'too-static' and
% 'teleportation' artifacts. The former usually happens when another blob
% is detected by mistake and stays fixed in a position. The later happens
% because the blobs cross each other and the tracker exchanges them,
% normally.

% if already calculated
if exist(fullfile(path, 'Dist2Feed.mat'),"file")
    load(fullfile(path, 'Dist2Feed.mat'));
else
    % Use a 'static' measure to determine when a blob is too still
    static = 0;
    % per feeder
    for ff = 1:size(feeder,1)
        % per frame (pre and post tutor) and blob ID
        for ID = 1:lastID
            for f = [2:timeOut(1) , timeOut(2):length(blob.Time)]
                
                % Euclidean distance calculation
                Dist2feeder(f,ID,ff) = pdist([squeeze(blob.Tracking(f,:,ID)); feeder(ff,:)]); 
    
                % static detection by comparing adjacent samples
                if Dist2feeder(f,ID,ff) == Dist2feeder(f-1,ID,ff)
                    static = static+1; % add up staticity
                    if static > 4 % limit reached
                        Dist2feeder(f-5:f-1,ID,ff) = nan; % remove static samples
                        static = 0; % reset
                    end
                else
                    static = 0; % reset
                end
    
                % teleportation checker, if distance suddenly increases over a limit
                if abs(Dist2feeder(f,ID,ff) - Dist2feeder(f-1,ID,ff)) > 20
                    Dist2feeder(f,ID,ff) = nan; % remove teleportation
                end
    
            end
    
            % Label the ID if it is empty after the calculations and first
            % cleaning. Later, will reduce use of memory.
            if ff == 1
                isnanidx(ID) = all(isnan(Dist2feeder(:,ID,ff)));
            end
        end
    end

    % % All data feeder 1;
    % subplot(5,1,1); plot(Dist2feeder(:,:,1));
    % xlim([1 length(blob.Time)]);  ylim([ 0 400]);
    
    % Destroy all empty IDs for all feeders
    if any(isnanidx), Dist2feeder(:,isnanidx,:) =  []; end
    clear isnanidx static
    
    save(fullfile(path, 'Dist2Feed.mat'), "Dist2feeder")
end

%% Detect double tracking at a given times
% 

% % Re-find last ID
% lastID = size(Dist2feeder,2);
t = 0; % initialize a double track counter
ff = 1; % Use only one feeder

% For each frame (pre and post tutor)
for f = [2:timeOut(1) , timeOut(2):length(blob.Time)]
    % Find times when more than one ID has distance value
    % and it did not before (beggining of a double tracking)
    if sum(~isnan(Dist2feeder(f,:,ff))) > 1 && sum(~isnan(Dist2feeder(f-1,:,ff))) < 2
        t = t+1;        % double track ordinal
        idx(1,t) = f;   % first frame of double tracking
    end

    % Find times when only one ID has a distance value 
    % and there were two before (end of a double tracking)
    if sum(~isnan(Dist2feeder(f,:,ff))) < 2 && sum(~isnan(Dist2feeder(f-1,:,ff))) > 1
        % uses same t ordinal
        idx(2,t) = f-1; % get the last frame when there were two (previous)
    end
end

% In (one, so far) weird cases, the double tracking begins right before timeOut (or end), getting
% lost inside tutor time and therefore its ending not being found. Force it.
if size(idx,1) < 2 % t has begins but no ends
    idx(2,:) = [timeOut(1) length(blob.Time)]; % set tutor entry time and session end
end

% Check the double tracking dynamics. We keep the most variant one, as it
% is most likely to be an active agent, and not an static artifact. 
% *This applies to all feeders, even though we are evaluating only one.
% DOES NOT WORK WELL FOR ALL POSSIBLE CASES.
dyn = [];
if ~isempty(idx) % there are double track periods
    for i = 1:t
        % Find IDs with distance value
        idsIdx = find(~isnan(Dist2feeder(idx(1,i),:,1))); 
        % Use IDs and detection periods to extract the trajectories of each
        % blob during this time.
        dyn(:,idsIdx) = Dist2feeder(idx(1,i):idx(2,i),idsIdx,1);

            % By comparing their range of values, we can estimate if one is
            % 'static' vs a dynamic one. Make the static invalid.
            % IN RARE CASES, THE DYNAMIC ONE IS A SPURIOUS TRAJECTORY
            % DIVERGING FROM THE ACTUAL BLOB, OR A TELEPORTATION. IT WILL FAIL.
            if range(dyn(:,idsIdx(1))) > range(dyn(:,idsIdx(2)))
                Dist2feeder(idx(1,i):idx(2,i),idsIdx(2),:) = nan;
            elseif range(dyn(:,idsIdx(1))) < range(dyn(:,idsIdx(2)))
                Dist2feeder(idx(1,i):idx(2,i),idsIdx(1),:) = nan;
            end

        clear dyn idsIdx
    end
end
clear idx t

% % All data feeder 1;
% subplot(5,1,2); plot(Dist2feeder(:,:,1));
% xlim([1 length(blob.Time)]);  ylim([ 0 400]);

%% Collapse all IDs
% Now, we have hopefully greatly reduced the overlap od tracking with 
% different IDs for a same time. Some samples may remain, tho. Ideally, the
% average will take the value of one ID per time, but if we have done it
% well enought before, the noise after averaging IDs should remain small.

% allocate
DistFeed = nan(length(blob.Time), size(feeder,1));

% We do this for all feeders and all valid samples
for ff = 1:size(feeder,1)
    for f = [1:timeOut(1) , timeOut(2):length(blob.Time)]
        % Only if there is one or more values
        notnan = ~isnan(Dist2feeder(f,:,ff));

        if sum(notnan) == 1 % Only one value found, use it
            DistFeed(f,ff) = Dist2feeder(f,notnan,ff);
        elseif sum(notnan) > 1 % More than one, average (PROB NOT BEST SOLUTION)
            DistFeed(f,ff) = mean(Dist2feeder(f,notnan,ff),2,'omitnan');
        end
    end
end
clear Dist2feeder

%% Fill gaps using autoregressive modeling
Dist2feeder([1:timeOut(1),timeOut(2):length(blob.Time)],:) = fillgaps(DistFeed([1:timeOut(1),timeOut(2):length(blob.Time)],:));
%clear DistFeed

%% Remove improbable zeroes
Dist2feeder(Dist2feeder==0) = nan;

%% Compute velocity.
% Sample time/pos each 60 points (i.e each second)
% Do so for each feeder and then average all computations
AgentVel = mean(abs(diff(Dist2feeder(1:Fs:end,:)) ./ diff(blob.Time(1:Fs:end))'), 2);

savepath = fullfile(input.analysis, input.subjects(input.run(1)).name, input.sessions(input.run(1)).list{input.run(2)});
save(fullfile(savepath, 'Dist2Feed_Vel.mat'), "Dist2feeder", "AgentVel")

%% Plotting
% Optional Plots for distance and velocity
fig = figure('visible', 'off'); % switch visibility
set(fig, 'Position', [100, 100, 800, 600]); % Set fig size as screen

subplot(2,1,1); 
    plot(Dist2feeder);
    xlim([1 length(blob.Time)]); 
    ylim([ 0 405]);
    xticks(1:60*Fs*5:length(blob.Time))
    xticklabels({0:5:ceil(length(blob.Time)/60*Fs)});
    legend('West/Left', 'West/Right', 'North/Left', 'North/Rigth', Location='best');
    legend('boxoff');
    ylabel('distance to feeder (px)');
    xlabel('time (min)');
    grid("on");

subplot(2,1,2);
    plot(AgentVel);
    xlim([1 length(blob.Time)/Fs]); 
    ylim([ 0 max(AgentVel)*1.1]);
    xticks([1:5*Fs:length(blob.Time)/Fs])
    xticklabels({0:5:ceil(length(blob.Time)/Fs)});
    ylabel('agent velocity (px/s)');
    xlabel('time (min)');
    grid("on");

% Set title and save
sgtitle('Blob distance to feeder(s) and velocity'); 
path = fullfile(input.bhvfolder, input.subjects(input.run(1)).name, input.sessions(input.run(1)).list{input.run(2)});
exportgraphics(fig, fullfile(path, 'Dist2Feed_Velocity.png'), 'Resolution', 300);
close all
end