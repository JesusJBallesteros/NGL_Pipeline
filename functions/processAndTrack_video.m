function [blob] = processAndTrack_video(opt)
% Extract euclidean distances between two subjects and tag interaction times.
% This function allows to treat and reduce videos from central cenital
% camera, to run a frame-by-frame blob tracking analysis and to extract the
% euclidean distance between a set of detected objects. If objects merge,
% social interaccion is assumed and tagged as such. Once they split, we
% resume the distance measure.
%
% Subject to improvement and modification.
% % TODO work (01) initial processing of video chunks into single video. 
%        work on individual distances to feeders (only tutee? single ID not possible yet)
%        

% Jesus 20.08.2024

if ~isfield(opt,'visualize') || isempty(opt.visualize),     opt.visualize = false; end
loaddata = false;

%% 00. Default/Custom Parameters
% Inputs
vidName = fullfile(opt.behavFiles, [opt.SavFileName '.mp4']);

% Explicit main variables here...
    mask = []; % series of vertices to hide all except for area of interest in arena
    trim = []; % values to trim down the video to contain only the action
    
    % Output structure. Will contain the input parameters as well.
    blob = struct('IDs', [], 'Pos', [], 'Time', [], ... % Basic variables from main blob analysis
                  'Par', [], ... % Default/custom parameters
                  'Tracking', [], 'euclDist', [], 'Merges', []); % Processed data of main interest

% ... or fill up main variables with a pre-made script
    blobConfig; % 'blobConfig.m' file should exist under '.../analysisCode'

%% 01. TODO: concatenate original 5-minute chunk videos into full session.
% Find pieces (sessionname_001 *_002 ...)
% Rotate 6 degrees (CW/CCW?) as Juan manually does that.
% Concatenate (sessionname_001 *_002 ...) into vidName: Create videoWriter 
%   object, create videoReader and sequentially read chinks, counter of frames 
%   adds up the frames from individual shorts into a longer count for full video.
% Crop video to include only from 'All LED flash to All LED flash', as session 
%   start and ends. Should be detectable. ** Perhaps in next step 'TreatVideo' **
% Save 'vidName'.

%% 02. Treat original video to extreme constrast and reduce size.
% Translate the original video arena masking values to the final ones.
arena(:,1) = mask(:,1)-trim(1);
arena(:,2) = -mask(:,2)+trim(2);

% Proceed
if ~isfile([vidName(1:end-4) '_t.mp4'])
    disp('Treating video...')
    vidName = TreatVideo(vidName, mask, trim); % updates video fine name
else
    disp('Treated video already exist.')
    vidName = [vidName(1:end-4) '_t.mp4'];
end

%% 03. Blob Analysis
% Get parameters out of blobconfig.m
blob.Par.visualize = opt.visualize;

% Proceed with processed video
if ~isfile([vidName(1:end-6) '_t_ID.mp4'])
    disp('Processing blob detection...')
    [blob.Pos, blob.IDs, blob.Time] = MotionBasedObjectTracking(vidName, blob.Par);
else
    disp('Blob tracked video already exists. Delete it to reprocess.')
end

% Proceed with processed video and parameters 
if isfile(fullfile(opt.behavFiles, "blob.mat"))
    disp('Blob tracked data already exists. Delete it to reprocess.')
    loaddata = true; % To reload data for debugging. Comment next line
    return % if we reach this point, nothing else happens
end

%% 04. Aggregate data into a 3d matrix with dim (frame, [x,y], ID)
% Find out frames with sucessful detection and assignment (ie, with ID)
if ~loaddata
    fidx = ~cellfun(@isempty, blob.IDs);
    
    % Make empty frames ID == 0
    blob.IDs(~fidx) = {0};
    
    % Last frame with assigned ID
    lastID = find(fidx, 1, "last");
    lastID = max(blob.IDs{lastID});
    
    % Allocate maximum 3-dim matrix
    blob.Tracking = nan(size(fidx,1), 2, lastID);
    
    % Re-order data for later plotting tracks in space
    for f = 1:size(blob.Tracking,1) % per frame
        if fidx(f) && ~isempty(blob.Pos{f}) % ID > 0 and valid centroid
            if size(blob.Pos{f}, 1) < numel(blob.IDs{f}) % more IDs than centroids
                for i = 0:size(blob.Pos{f}, 1)-1
                    blob.Tracking(f,:,blob.IDs{f}(end-i)) = blob.Pos{f}(i+1,:);
                end
            else % as many IDs as centroids (or more)
                for i = 0:numel(blob.IDs{f})-1
                    blob.Tracking(f,:,blob.IDs{f}(i+1)) = blob.Pos{f}(end-i,:);
                end
            end
        end
    end
else
    % already processed data
    load(fullfile(opt.behavFiles, "blob.mat"), 'blob');
    
    % Re-find last ID
    fidx = ~cellfun(@isempty, blob.IDs);
    % Last frame with assigned ID
    lastID = find(fidx, 1, "last");
    lastID = max(blob.IDs{lastID});
end
clear fidx f i

%% 05. Analyze tracks and IDs to extract track merging state
if ~isfile(fullfile(opt.behavFiles, "blob.mat"))
% Create vectors: 
%   IDs = numIDs, number of objects in frame (int)
%   euclDist = euclidean distance between objects centroids (double) [0:max]
%   mergeState = merge state, 0|1
% only those 0-distances with two objects are considered merge (interaction) 
goodID        = zeros(lastID,1);
IDs           = zeros(length(blob.Time),1);
blob.euclDist = IDs;
mergeState    = IDs;

% Find good IDs 
for ID = 1:lastID
    % frames with this ID
    lidx = cellfun(@(x) ismember(ID,x), blob.IDs); 
    % if it spans longer than 180 frames (3s)
    if any(lidx) && sum(lidx) > blob.Par.minVisibleCount*3
          goodID(ID) = 1; % add good ID
    else, goodID(ID) = 0;
    end
end

% Evaluate each frame (f(1) always == [0 0 0])
for f = 2:length(blob.Time) 
    % If we are in the first 5 minutes, any second object is an artifact
    if f < 60*300 && sum(blob.IDs{f} > 1)
        IDs(f) = 1; % force to 1
        % blob.euclDist(f) = 0; % already 0
        % mergeState(f) = 0; % already 0 
    else
        IDs(f) = sum(blob.IDs{f}~=0); % get number of objects (not counting 0s)

        % If there is no objects, var(f) remains as the previous frame
        if IDs(f) == 0 % Avoiding some jitter of blob/centroid detection
            blob.euclDist(f) = blob.euclDist(f-1);
            mergeState(f) = mergeState(f-1);

        % If there is ONLY one object
        elseif IDs(f) == 1
            % to consider a merge, prior distance small and either previous frame 
            % were 2 objects or it was a merge
            if IDs(f-1) == 2 || mergeState(f-1) == 1
                mergeState(f) = 1; % start/continue the merge state;

                if mergeState(f-1) ~= 1 && any(blob.euclDist(f-31:f-1) < 80)
                    mergeState(f) = 1; % allow start a merge state if distance was low;
                end
            end
            % if not, it must be an actual single object, or one object got
            % lost by mislabeling or other reasons
    
        % if there are TWO objects and TWO centroids
        elseif IDs(f) == 2 && size(blob.Pos{f},1) == 2 
            % Get Euclidean distance btw centroids
            blob.euclDist(f) = pdist(blob.Pos{f}); 
            % Force merge state to 0, as clearly separated objects.
            mergeState(f) = 0; 
    
        % if TWO objects BUT mismatching number of centroids
        elseif IDs(f) == 2 && size(blob.Pos{f},1) ~= 2 
            % Keep merge state as previous. 
            mergeState(f) = mergeState(f-1); 
            % Sometimes the centroid prediction fails, or take a number of
            % frames to predict. My guess is that we could skip the centroid
            % perdiction and work with a direct calculation from the area (or
            % box) already detected, but not sure if really necessary. 
            blob.euclDist(f) = blob.euclDist(f-1);

        % there are MORE than two objects.
        elseif IDs(f) > 2
            % At times in these sort of videos there are new elements appearing
            % briefly (human, doors, shadows...) or artifactual detections. The
            % parameters should try to avoid or reduce these, but it may be
            % inevitable to a point. In this case just assume no merge and hope
            % for a brief period mislabeled.
            mergeState(f) = 0; % Force merge state to 0
            blob.euclDist(f) = nan;

        end
    
        % After tutor time is up, Any possible artifact resulting in merge
        % is removed. This can happen due to different effects after 
        % having various blobs and then the tutor leaving, existing only
        % one blob until the end. Also removes the last merge if it was 
        % started within 30 seconds before this time.
        if f == 60*960 && mergeState(f)==1
            mergeState(f) = 0;
            solotime = find(IDs(f-3600:f-1)~=1,1,'last');
            mergeState(f-3600+solotime:f-1) = 0;
        end

        % Last check at last frame. If we have been carrying a merge state
        % by mistake, we will remove it since it started.
        if f == length(mergeState) && mergeState(f)==1
            solotime = find(mergeState~=1,1,'last');
            mergeState(solotime:f) = 0;
        end

    end
end

% Identify merge start/end times
fidx = NaN(2,size(blob.Time(mergeState==1),2));
fidx(1,:) = blob.Time(mergeState==1);       % frames idx as merge
fidx(2,1) = 1;                              % first fidx can only be a merge start
fidx(2,2:end) = diff(fidx(1,:)) > 0.0167;   % non-continuous merge states belong to different merges (over 1 sec)
fidx(2,end) = 2;                            % last fidx can only be a merge end
endmerge = fidx(2,2:end)==1;                % check all merge starts and locate ends right before
fidx(2,endmerge) = 2;                       % tag those positions as merge end

% now get the times associated to 1s and 2s as merge start/ends
blob.Merges = [fidx(1,fidx(2,:)==1); fidx(1,fidx(2,:)==2)]';

clear fidx endmerge f

%% 06. Figures
% IDs lifespan and merge states
fig = figure('visible', 'off'); % switch visibility
set(fig, 'Position', [0, 0, 1600, 900]); % Set fig size as screen
label = cell(numel(goodID),1);

h = subplot(2,2,1);
    for ID = 1:lastID
        % frames with this ID
        lidx = cellfun(@(x) ismember(ID,x), blob.IDs); 
        if goodID(ID)
            label{sum(goodID(1:ID))} = int2str(ID); % take it's ordinal
            scatter(h, blob.Time, int32(lidx)*sum(goodID(1:ID)), 5); % plot its presence
            hold on
        end
    end 
    
    % Add merge periods as black transparent patches
    for r = 1:size(blob.Merges,1)
        rectangle('Position', [blob.Merges(r,1) 0 blob.Merges(r,2)-blob.Merges(r,1) lastID+1], ...
            'FaceColor', [0, 0, 0, 0.3], 'EdgeColor', 'none');
    end
    
    % Adjust
    ylim([0.5 sum(goodID)+0.5]); ylabel('Good IDs'); % hide IDs == 0 
    yticks(1:1:sum(goodID)); yticklabels(label);
    xlim([blob.Time(1) blob.Time(end)]); xlabel ('min');
    xticks(0:300:blob.Time(end)); xticklabels(string(xticks/60))
    cmap = h.ColorOrder; % Get colors and order used

% Centroids trajectory
j = subplot(2,2,2); 
    j.ColorOrder = cmap;
    for ID = 1:lastID
        if goodID(ID)
            trackID = squeeze(blob.Tracking(:,:,ID));    
            plot(j, trackID(:,1),-trackID(:,2));
            hold on
        end
    end
    % Add arena contour. Close it by concatenating first element at the end.
    plot(j, [arena(:,1); arena(1:1)] , [arena(:,2); arena(1,2)], 'LineWidth',2, 'Color', 'k');
    xlim([0 550]), ylim([-550 0]);
    set(j,'XAxisLocation','top');
    xlabel('width (px)'); ylabel('hight (px)');

k = subplot(2,2,[3 4]); 
    % Plot simple time series of euclidean distance between blobs
    plot(blob.Time, blob.euclDist, 'LineStyle','-','Color','k','LineWidth',1)
    box off
    hold on
    % Add merge periods as black transparent patches
    for r = 1:size(blob.Merges,1)
        rectangle('Position', [blob.Merges(r,1) 0 blob.Merges(r,2)-blob.Merges(r,1) max(blob.euclDist)*1.1], ...
            'FaceColor', [0, 0, 0, 0.3], 'EdgeColor', 'none');
    
    end
    % Remove solo times, i.e show only 250 < t < 1000. Add labels.
    xlim([250 1000]); xlabel('min');
    xticks(0:60:blob.Time(end));  xticklabels(string(xticks/60))
    ylim([-2 max(blob.euclDist)*1.1]); ylabel('Eucl. dist. (px)');  

exportgraphics(fig,fullfile(opt.behavFiles,'blobTracking.png'),'Resolution',600);
close all

end