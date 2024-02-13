function Deuteron_PlotMotionSensors(data, timestamps, opt, ver, varargin)
% Description is progress
%
%
%
% Jesus

%% Check input variables
if nargin < 5
    visual = 0; 
    record = 0; 
elseif nargin == 5
    visual = varargin{1}; % ver dynamic figure to visualize the rotation.
    record = 0; 
elseif nargin == 6
    visual = varargin{1};
    record = varargin{2}; % If 'visual'=1, save it as avi.
end

% defaults
fs = 1000; % sampling rate.
% framerate = 1/50; % As 1/Hz of pause for next measurement. Default to 50Hz

%% Sensor readings
if ver == 1
    f = figure;
    for i=1:3
        if      i==1, MSData = data.acc; tit = 'Accelerometer'; units = 'm/s^2';
        elseif  i==2, MSData = data.gyr; tit = 'Gyroscope';     units = 'deg/s';
        elseif  i==3, MSData = data.mag; tit = 'Magnetometer';  units = 'Tesla';
        end
    
        if ~isempty(MSData)
            subplot(3,1,i)
            % General ver
            plot(timestamps, MSData.X); hold on
            plot(timestamps, MSData.Y); hold on
            plot(timestamps, MSData.Z);
            xlim([timestamps(1) timestamps(end)]);
            ylabel(units);
            if i<3, ylim([-MSData.max*1.1 MSData.max*1.1]);  
            else,   ylim([-1e-4 1e-4]); xlabel('ms');
            end
            legend({'X' 'Y' 'Z'}, 'Box', 'off');
            title(tit);
        end
    end
    f.Position(3:4) = [1600 600];
    % Save
    exportgraphics(gcf, fullfile(opt.FolderProcDataMat, strcat('MotionRaw.png')), 'Resolution', 300)
    close gcf

elseif ver == 2
    f = figure;
    for i=1:3
        if      i==1, MSData = data.acc; tit = 'Accelerometer'; units = 'm/s^2';
        elseif  i==2, MSData = data.gyr; tit = 'Gyroscope';     units = 'deg/s';
        elseif  i==3, MSData = data.magcorr; tit = 'Mag Corrected';  units = 'Tesla';
        end
    
        if ~isempty(MSData)
            subplot(3,1,i)
            % General ver
            plot(timestamps, MSData.X); hold on
            plot(timestamps, MSData.Y); hold on
            plot(timestamps, MSData.Z);
            xlim([timestamps(1) timestamps(end)]);
            ylabel(units);
            if i<3, ylim([-MSData.max*1.1 MSData.max*1.1]);  
            else,   ylim([-1e-4 1e-4]); xlabel('ms');
            end
            legend({'X' 'Y' 'Z'}, 'Box', 'off');
            title(tit);
        end
    end
    f.Position(3:4) = [1600 600];

    % Save
    exportgraphics(gcf, fullfile(opt.FolderProcDataMat, strcat('MotionRaw_corrected.png')),'Resolution',300)
    close gcf

elseif ver == 3
    %% Plot the orientation in Euler angles in degrees over time.
    orientation_v = eulerd(data, 'ZYX', 'frame');
    timeVector = (0:length(timestamps)-1).'/fs;

    f = figure;
    plot(timeVector, orientation_v);
        xlim([timeVector(1) timeVector(end)]);
        xlabel('Time (s)')
        ylabel('Rotation (degrees)')
    legend({'X, Roll', 'Y, Yaw', 'Z, Pitch'}, 'Box','off');
    title('Orientation Estimate');
    f.Position(3:4) = [1600 600];

    exportgraphics(gcf, fullfile(opt.FolderProcDataMat, strcat('ecomp_OrEstim.png')), 'Resolution', 300)
    close gcf

    %% Plot the helper viewer example from MATLAB
    % Initialize objects and set timer.
    if visual
        stopTimer = timeVector(end); % seconds to run simulation
        
        % Creates a very specific figure object provided by Matlab.
        pp = poseplot;
            title("Pose, NED");
            xlabel('North')
            ylabel('East')
            zlabel('Down')

        % initialize the VideoWriter object.
        if record
            writerObj = VideoWriter(fullfile(opt.FolderProcDataMat, strcat('ecomp_PoseStim.avi')),'Motion JPEG AVI'); 
            open(writerObj); % Opens the file.
        end
    
        % Timer
        ts = tic; % start timer
        pause(0.0009) % Let clock tic to a first milisecond
        % Run until elapsed time reaches set 'stopTimer' (-.5 sec to avoid breaks)
        while(toc(ts) < stopTimer-0.5) 
            t = round(toc(ts)*1000); % takes the approximated msec of the run.
            % plot it in the dynamic figure
            set(pp, "Orientation", data(t))
            drawnow limitrate
            % Get frame and write it to video.
            if record
                F = getframe(gcf);           % Capture the frame
                writeVideo(writerObj, F) % add the frame to the movie
            end
        end
    
        % Close video file.
        if record
            close(writerObj);
        end
    end

%     %% plot Dynamic figure where the heading vector moves as the rotation happens
%     % Assuming you have an array of quaternions "rotators" with dimensions (n, 4)
%     % where n is the number of time points. Convert the quaternions to rotation matrices
% %     rotMat = quat2rotm(data);
% 
%     % Define an initial vector.
%     curr_pos = [1, 1, 0]; % do not use [0, 0, 0]
%     
%     if visual
%         stopTimer = timeVector(end); % seconds to run simulation
% 
%         if record
%             % Initialize the VideoWriter object.
%             writerObj = VideoWriter(fullfile(opt.FolderProcDataMat, strcat('ecomp_HeadingStim.avi')),'Motion JPEG AVI'); 
%             open(writerObj); % Opens the file.
%         end
%         % Create a figure with initial, non-visible vector.
%         c = quiver3(0,0,0,0,0,0,'off');
%             c.Color = 'r'; % Arrow color
%             c.MaxHeadSize = 2; % Arroy head size
%             c.LineWidth = 2; % Arroy line width
%             xlabel('X');  ylabel('Y');  zlabel('Z'); % Label axes
%             xlim([-2,2]); ylim([-2,2]); zlim([-2,2]); % Fix axes scale
%             daspect([1 1 1]); % Set the aspect ratio to be equal.
%             title('Estimated Heading');
%         
%         % Timer
%         ts = tic; % start timer
%         pause(0.0009) % Let clock tic to a first milisecond
%         % Run until elapsed time reaches set 'stopTimer' (-5 msec to avoid breaks)
%         while(toc(ts) < stopTimer-0.5) 
%             t = round(toc(ts)*1000); % takes the approximated msec of the run.
%             % Get the rotation matrix at the current time point and rotate the previous position by the rotation matrix.
%             txt = ['Time: ', num2str(t/1000), ' sec'];
% %             curr_pos = rotMat(:,:,t) * curr_pos; % Get current position and rotate according to 'rotMat' step
%             curr_pos = rotatepoint(data(t), curr_pos); % Get current position and rotate according to 'quaternion' step
%             % Collect new datapoints 
%             c.UData = curr_pos(1);
%             c.VData = curr_pos(2);
%             c.WData = curr_pos(3);
%             
%             % Update figure.
%             drawnow;
%                 title(txt);
% 
%             % Get the frame and write it to video.
%             if record
%                 F = getframe;           % Capture the frame
%                 writeVideo(writerObj,F) % add the frame to the movie
%             end
%             pause(framerate) % pause the run to an approx. framerate.
%         end
%     
%         if record
%             % Close video file.
%             close(writerObj);
%         end
%     end
end

end
