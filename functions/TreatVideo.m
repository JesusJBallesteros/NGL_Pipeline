function newVidName = TreatVideo(vidName, mask, trim)
% Go through video, create and apply mask, trim frames and create new video file.
% The mask is a polygon adaptable to any possible arena.
% The trim is a rectangle that contains the mask, fitted to a reduced size.
%
% Jesus 17.07.2024

%% Video to Process
readVid = VideoReader(vidName);
    numberOfFrames = readVid.NumFrames;
    vidHeight = readVid.Height;
    vidWidth = readVid.Width;
       
%% Processed video
newVidName= [vidName(1:end-4), '_t.mp4'];
writeVid = VideoWriter(newVidName, 'MPEG-4');
    writeVid.FrameRate = 60;
    open(writeVid);
    
%% Go through frames 
for frame = 1:numberOfFrames
    thisFrame = read(readVid, frame); % Extract frame
    thisFrame = rgb2gray(thisFrame);  % Conver to gray scale
    thisFrame = imadjust(thisFrame, [0.45 0.55]); % adjust intensity levels

    % Create mask only on first frame
    if frame == 1
        figure('Visible','off');
        imshow(thisFrame)
        % Create a polygonal ROI
        h = images.roi.Polygon(gca,'Position', mask);
        % Create a ROI based mask
        BW = createMask(h, vidHeight, vidWidth);
        close all
    end

    thisFrame = thisFrame.*uint8(BW);    % Apply mask to frame
    framecut = imcrop(thisFrame, trim);  % Crop frame to reduce size
    writeVideo(writeVid, framecut);      % Write frame
end

% Finish writing new file
close(writeVid)
end