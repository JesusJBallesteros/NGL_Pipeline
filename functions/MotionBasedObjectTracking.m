function [blobPos, blobIDs, blobTime] = MotionBasedObjectTracking(vidName, par)
%% Motion-Based Multiple Object Tracking
% This example shows how to perform automatic detection and motion-based
% tracking of moving objects in a video from a stationary camera.
%   Copyright 2014 The MathWorks, Inc.
% Detection of moving objects and motion-based tracking are important 
% components of many computer vision applications, including activity
% recognition, traffic monitoring, and automotive safety.  The problem 
% of motion-based object tracking can be divided into two parts:
%
% # Detecting moving objects in each frame 
% # Associating the detections corresponding to the same object over time
%
% The detection of moving objects uses a background subtraction algorithm
% based on Gaussian mixture models. Morphological operations are applied to
% the resulting foreground mask to eliminate noise. Finally, blob analysis
% detects groups of connected pixels, which are likely to correspond to
% moving objects. 
%
% The association of detections to the same object is based solely on
% motion. The motion of each track is estimated by a Kalman filter. The
% filter is used to predict the track's location in each frame, and
% determine the likelihood of each detection being assigned to each 
% track.
%
% Track maintenance becomes an important aspect of this example. In any
% given frame, some detections may be assigned to tracks, while other
% detections and tracks may remain unassigned. The assigned tracks are
% updated using the corresponding detections. The unassigned tracks are 
% marked invisible. An unassigned detection begins a new track. 
%
% Each track keeps count of the number of consecutive frames, where it
% remained unassigned. If the count exceeds a specified threshold, the
% example assumes that the object left the field of view and it deletes 
% the track.  
%
% % Jesus 17.07.2024

%% Create System objects used for reading video, detecting moving objects,
% and displaying the results.
par.newVidName = [vidName(1:end-4), '_ID.mp4'];
obj = setupSystemObjects(vidName, par);

% Create an empty array of tracks and initialize tracks ID
tracks = initializeTracks();
nextId = 1;

% Frames counter and output variables
Nfs = obj.reader.NumFrames;
f = 1;
blobIDs  = cell(Nfs,1);
blobPos  = cell(Nfs,1);
blobTime = zeros(1,Nfs);

%% Main process. Detect moving objects, and track them across video frames.
while hasFrame(obj.reader)
    frame = readFrame(obj.reader);  % Get frame

    [centroids, bboxes, mask] = detectObjects(frame); % Process frame to locate blobs

    predictNewLocationsOfTracks();  % If blobs go invisible, predict its location

    [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment(centroids, par);

    updateAssignedTracks(centroids);

    updateUnassignedTracks();

    deleteLostTracks(par);

    createNewTracks(centroids, par);

    ids = displayTrackingResults(centroids, par);

    writeVideo(obj.writeVid, frame);

    % Collect frame data
    blobIDs{f} = ids;
    blobPos{f} = centroids;
    if isempty(centroids) && ~isempty(bboxes)
        blobPos{f} = centroids;
    end
    blobTime(1,f) = obj.reader.CurrentTime-(1/60);
    f = f + 1;
end

%% Initialize functions 
% Create System Objects
function obj = setupSystemObjects(vid_file, par)
% Initialize Video I/O. Create objects for reading a video from a file,
% drawing the tracked objects in each frame, and playing the video.
    
    % Video Reader.
    obj.reader = VideoReader(vid_file);
    
    % Video Players.
    if par.visualize
        obj.maskPlayer  = vision.VideoPlayer('Position', [2029 -248 654 683]); % specific for a vertical second screen
        obj.videoPlayer = vision.VideoPlayer('Position', [2029 428 678 676]); % specific for a vertical second screen
    end

    % The foreground detector is used to segment moving objects from the background. 
    % It outputs a binary mask, where the pixel value of 1 corresponds to the 
    % foreground and the value of 0 corresponds to the background. 
    obj.detector = vision.ForegroundDetector('NumGaussians', par.NumGaussians, 'NumTrainingFrames', par.NumTrainingFrames, ...
        'MinimumBackgroundRatio', par.MinimumBackgroundRatio, 'LearningRate', par.LearningRate);
    
    % Connected groups of foreground pixels are likely to correspond to moving objects.  
    % The blob analysis System object finds such groups ('blobs'), and compute their 
    % characteristics, such as area, centroid, and the bounding box.
    obj.blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort', par.BoundingBoxOutputPort, 'AreaOutputPort', par.AreaOutputPort, ...
        'CentroidOutputPort', par.CentroidOutputPort, 'MinimumBlobArea', par.MinimumBlobArea, 'MaximumBlobArea', par.MaximumBlobArea);

    % Object for writing the ID'd video
    obj.writeVid = VideoWriter(par.newVidName, 'MPEG-4');
        obj.writeVid.FrameRate = 60;
        open(obj.writeVid);

end

% Initialize Tracks
function tracks = initializeTracks()
% Creates an array of tracks, each is a structure representing a moving object
% in the video. The purpose is to maintain the state of a tracked object, consisting 
% of information used for detection, track assignment, track termination and display. 
%
% The structure contains the following fields:
% id:                  the integer ID of the track
% bbox:                the current bounding box of the object; used for display
% kalmanFilter:        a Kalman filter object used for motion-based tracking
% age:                 the number of frames since the track was first detected
% totalVisibleCount:   the total number of frames in which the track was detected (visible)
% consecutiveInvisibleCount: the number of consecutive frames for which the track was not detected (invisible).
%
% Noisy detections tend to result in short-lived tracks. This happens when 
% 'totalVisibleCount' exceeds a specified threshold.    
%
% When no detections are associated with a track for several consecutive frames, 
% is assumed that the object has left the field of view and deletes the track. 
% This happens when 'consecutiveInvisibleCount' exceeds a specified threshold. 
% A track may also get deleted as noise if it was tracked for a short time, 
% and marked invisible for most of the frames.        
tracks = struct(...
                'id',                           {}, ...
                'bbox',                         {}, ...
                'kalmanFilter',                 {}, ...
                'age',                          {}, ...
                'totalVisibleCount',            {}, ...
                'consecutiveInvisibleCount',    {});
end

%% Tracking functions
% Detect Objects
function [centroids, bboxes, mask] = detectObjects(frame)
% Returns the centroids, the bounding boxes and the binary mask, of the detected
% objects. Pixels with a value of 1 correspond to the foreground, and pixels with
% a value of 0 correspond to the background. The function performs motion segmentation
% using the foreground detector. It then performs morphological operations on the
% resulting binary mask to remove noisy pixels and to fill the holes in the remaining blobs.  
    
    % Detect foreground.
    mask = obj.detector.step(frame);

    % Apply morphological operations.
    mask = imopen(mask, strel('rectangle', [3, 3]));
    mask = imclose(mask, strel('rectangle', [8, 8])); % [15, 15]
    mask = imfill(mask, 'holes');
    
    % Blob analysis of determined masks
    [~, centroids, bboxes] = obj.blobAnalyser.step(mask);
end

% Predict New Locations of Existing Tracks
function predictNewLocationsOfTracks()
% Use the Kalman filter to predict the centroid of each track in the
% current frame, and update its bounding box accordingly.

    for i = 1:length(tracks)
        bbox = tracks(i).bbox;
        
        % Predict the current location of the track.
        predictedCentroid = predict(tracks(i).kalmanFilter);
        
        % Shift the bounding box so that its center is at 
        % the predicted location.
        predictedCentroid = int32(predictedCentroid) - bbox(3:4) / 2;
        tracks(i).bbox = [predictedCentroid, bbox(3:4)];
    end
end

% Assign Detections to Tracks
function [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment(centroids, par)
% Assigning object detections in the current frame to existing tracks is done by minimizing cost.
% The cost is defined as the negative log-likelihood of a detection corresponding to a track.  
%
% Step 1: Compute the cost of assigning every detection to each track using the 'distance'
% method. The cost takes into account the Euclidean distance between the predicted centroid
% of the track and the centroid of the detection. It also includes the confidence of the 
% prediction, which is maintained by the Kalman filter. The results are an MxN matrix, where 
% M is the number of tracks, and N is the number of detections.   
%
% Step 2: Solve the assignment problem represented by the cost matrix using the 
% 'assignDetectionsToTracks' function. The function takes the cost matrix and the cost of
% not assigning any detections to a track.  
%
% The value for the cost of not assigning a detection to a track depends on the range of values
% returned by the 'distance' method of the  'vision.KalmanFilter'. This value must be tuned 
% experimentally. Setting it too low increases the likelihood of creating a new track, and may
% result in track fragmentation. Setting it too high may result in a single track corresponding 
% to a series of separate moving objects.   
%
% It returns an M x 2 matrix containing the corresponding indices of assigned tracks and
% detections in its two columns. It also returns the indices of tracks and detections
% that remained unassigned.         
    nTracks = length(tracks);
    nDetections = size(centroids, 1);
    
    % Compute the cost of assigning each detection to each track.
    cost = zeros(nTracks, nDetections);
    for i = 1:nTracks
        cost(i, :) = distance(tracks(i).kalmanFilter, centroids);
    end
    
    % Solve the assignment problem.
    [assignments, unassignedTracks, unassignedDetections] = assignDetectionsToTracks(cost, par.costOfNonAssignment);
end

% Update Assigned Tracks
function updateAssignedTracks(centroids)
% The |updateAssignedTracks| function updates each assigned track with the
% corresponding detection. It calls the |correct| method of
% |vision.KalmanFilter| to correct the location estimate. Next, it stores
% the new bounding box, and increases the age of the track and the total
% visible count by 1. Finally, the function sets the invisible count to 0. 

    numAssignedTracks = size(assignments, 1);
    for i = 1:numAssignedTracks
        trackIdx = assignments(i, 1);
        detectionIdx = assignments(i, 2);
        centroid = centroids(detectionIdx, :);
        bbox = bboxes(detectionIdx, :);
        
        % Correct the estimate of the object's location using the new detection.
        correct(tracks(trackIdx).kalmanFilter, centroid);
        
        % Replace predicted bounding box with detected bounding box.
        tracks(trackIdx).bbox = bbox;
        
        % Update track's age.
        tracks(trackIdx).age = tracks(trackIdx).age + 1;
        
        % Update visibility.
        tracks(trackIdx).totalVisibleCount = tracks(trackIdx).totalVisibleCount + 1;
        tracks(trackIdx).consecutiveInvisibleCount = 0;
    end
end

% Update Unassigned Tracks
function updateUnassignedTracks()
% Mark each unassigned track as invisible, and increase its age by 1.
    for i = 1:length(unassignedTracks)
        ind = unassignedTracks(i);
        tracks(ind).age = tracks(ind).age + 1;
        tracks(ind).consecutiveInvisibleCount = tracks(ind).consecutiveInvisibleCount + 1;
    end
end

% Delete Lost Tracks
function deleteLostTracks(par)
% The |deleteLostTracks| function deletes tracks that have been invisible
% for too many consecutive frames. It also deletes recently created tracks
% that have been invisible for too many frames overall. 

    if isempty(tracks)
        return;
    end
        
    % Compute the fraction of the track's age for which it was visible.
    ages = [tracks(:).age];
    totalVisibleCounts = [tracks(:).totalVisibleCount];
    visibility = totalVisibleCounts ./ ages;
    
    % Find the indices of 'lost' tracks.
    lostInds =  (ages < par.ageThreshold & visibility < 0.7) | ...
                [tracks(:).consecutiveInvisibleCount] >= par.invisibleForTooLong;
    
    % Delete lost tracks.
    tracks = tracks(~lostInds);
end

% % Delete Static Tracks
% function deleteLostTracks(par)
% % The deleteStaticTracks function deletes tracks that have been static
% % for too many consecutive frames.
% 
%     if isempty(tracks)
%         return;
%     end
%         
%     % Compute the fraction of the track's age for which it was static.
%     ages = [tracks(:).age];
%     totalVisibleCounts = [tracks(:).totalVisibleCount];
%     visibility = totalVisibleCounts ./ ages;
%     
%     % Find the indices of 'lost' tracks.
%     lostInds =  (ages < par.ageThreshold & visibility < 0.7) | ...
%                 [tracks(:).consecutiveInvisibleCount] >= par.invisibleForTooLong;
%     
%     % Delete lost tracks.
%     tracks = tracks(~lostInds);
% end

% Create New Tracks
function createNewTracks(centroids, par)
% Create new tracks from unassigned detections. Assume that any unassigned
% detection is a start of a new track. In practice, you can use other cues
% to eliminate noisy detections, such as size, location, or appearance.

    centroids = centroids(unassignedDetections, :);
    bboxes = bboxes(unassignedDetections, :);
    
    for i = 1:size(centroids, 1)
        centroid = centroids(i,:);
        bbox = bboxes(i, :);
        
        % Create a Kalman filter object.
        kalmanFilter = configureKalmanFilter('ConstantVelocity', centroid, ...
                    par.InitEstErr, par.MotionNoise, par.MeasurementNoise);
        
        % Create a new track.
        newTrack = struct(...
            'id', nextId, ...
            'bbox', bbox, ...
            'kalmanFilter', kalmanFilter, ... % 
            'age', 1, ...
            'totalVisibleCount', 1, ...
            'consecutiveInvisibleCount', 0);
        
        % Add it to the array of tracks.
        tracks(end + 1) = newTrack;
        
        % Increment the next id.
        nextId = nextId + 1;
    end
end

%% Display Function
function ids = displayTrackingResults(centroids, par)
% The |displayTrackingResults| function draws a bounding box and label ID 
% for each track on the video frame and the foreground mask. It then 
% displays the frame and the mask in their respective video players. 
    % Convert the frame and the mask to uint8 RGB.
    frame = im2uint8(frame);
    mask = uint8(repmat(mask, [1, 1, 3])) .* 255;
    ids = [];

    if ~isempty(tracks)
          
        % Noisy detections tend to result in short-lived tracks.
        % Only display tracks that have been visible for more than 
        % a minimum number of frames.
        reliableTrackInds = [tracks(:).totalVisibleCount] > par.minVisibleCount;
        reliableTracks = tracks(reliableTrackInds);
        
        % Display the objects. If an object has not been detected
        % in this frame, display its predicted bounding box.
        if ~isempty(reliableTracks)
            % Get bounding boxes.
            bboxes = cat(1, reliableTracks.bbox);
            
            % Get ids.
            ids = int32([reliableTracks(:).id]);
            
            % Create labels for objects indicating the ones for 
            % which we display the predicted rather than the actual 
            % location.
            labels = cellstr(int2str(ids'));
            
            r = ones(size(centroids,1),1)*2;
            % Draw the objects on the frame and mask
            frame = insertObjectAnnotation(frame, 'rectangle', bboxes, labels);
            if length(labels) == size(centroids,1)
                frame = insertObjectAnnotation(frame, 'circle', [centroids r], labels);
            end
            mask = insertObjectAnnotation(mask, 'rectangle', bboxes, labels);
        end
    end
    
    % Display the mask and the frame.
    if par.visualize
        obj.maskPlayer.step(mask);        
        obj.videoPlayer.step(frame);
    end
end

end
