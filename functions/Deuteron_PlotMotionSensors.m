function Deuteron_PlotMotionSensors(data, opt, timestamps, stream, varargin)
% Description is progress
%
%
%
% Jesus

%% Check input variables
if nargin < 5
    visual = 0; % Plot dynamic figure to visualize the rotation.
    record = 0; % If plotting the figure, create a video from it. If 'visual' = 0, this has no effect.
elseif nargin == 5
    visual = varargin{1};
    record = 0; 
elseif nargin == 6
    visual = varargin{1};
    record = varargin{2}; 
end

%% Sensor readings
if stream == 1
    figure,
    for i=1:3
        if      i==1, MSData = data.acc; tit = 'Accelerometer'; units = 'm/s^2';
        elseif  i==2, MSData = data.gyr; tit = 'Gyroscope';     units = 'deg/s';
        elseif  i==3, MSData = data.mag; tit = 'Magnetometer';  units = 'Tesla';
        end
    
        if ~isempty(MSData)
            subplot(3,1,i)
            % General plot
            plot(timestamps, MSData.X); hold on
            plot(timestamps, MSData.Y); hold on
            plot(timestamps, MSData.Z);
            xlim([timestamps(1) timestamps(end)]);
            ylabel(units);
            if i<3, ylim([-MSData.max*1.1 MSData.max*1.1]);  
            else,   ylim([-1e-4 1e-4]); xlabel('ms');
                    legend({'X' 'Y' 'Z'}, 'Box', 'off');
            end
            title(tit);
        end
    end

elseif stream == 2
    %% Plot the orientation in Euler angles in degrees over time.
<<<<<<< HEAD
    orientation_v = eulerd(data, 'ZYX', 'frame');
    timeVector = (0:length(timestamps)-1).'/1000;

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
=======
    time = timestamps-(timestamps(1));
    orientation_deg = eulerd(data, 'ZYX', 'frame');
>>>>>>> d38898a4857c88bbc90db712ff88b029452ef85b

    figure,
    plot(time/1000,orientation_deg);
        xlim([time(1)/1000 time(end)/1000]);
        xlabel('Time (sec)');
        ylabel('Rot. around axis (deg)');
        legend({'Z, Pitch', 'Y, Yaw', 'X, Roll'}, 'Box','off');
        title('Orientation Estimate');
    
    %% Plot the helper viewer example from MATLAB
    % Initialize objects and set timer.
    if visual
        stopTimer = (timestamps(end)-timestamps(1))/1000; % seconds to run simulation
        framerate = 1/50; % As 1/Hz of pause for next frame. Default to 50Hz.
    
        % Creates a very specific figure object provided by Matlab.
<<<<<<< HEAD
        pp = poseplot(data(1));
            title("Pose, NED");
            xlabel('North')
            ylabel('East')
            zlabel('Down')
=======
        viewer = HelperOrientationViewer('Title',{'AHRS Filter'});
>>>>>>> d38898a4857c88bbc90db712ff88b029452ef85b

        if record
            % initialize the VideoWriter object.
            writerObj = VideoWriter(fullfile(opt.FolderProcDataMat, strcat('Magnetometer_rotations.avi')),'Motion JPEG AVI'); 
            open(writerObj); % Opens the file.
        end
    
        % Timer
        ts = tic; % start timer
        pause(0.0009) % Let clock tic to a first milisecond
        % Run until elapsed time reaches set 'stopTimer' (-1 msec to avoid breaks)
        while(toc(ts) < stopTimer-0.001) 
            t = round(toc(ts)*1000); % takes the approximated msec of the run.
<<<<<<< HEAD
            % plot it in the dynamic figure
            set(pp, "Orientation", data(t))
            subtitle(sprintf('%0.3f s', t/1000))
            drawnow limitrate
            % Get frame and write it to video.
            if record
                F = getframe(gcf); % Capture the frame
                writeVideo(writerObj, F) % add the frame to the movie
            end
            % pause(1/50);
=======
            % Plot it in the dynamic figure
            viewer(data(t));
            % If desired, get the frame and write it to video.
            if record
                F = getframe;           % Capture the frame
                writeVideo(writerObj, F) % add the frame to the movie
            end
            pause(framerate) % pause the run for 1/Hz msec, to an approx framerate.
>>>>>>> d38898a4857c88bbc90db712ff88b029452ef85b
        end
    
        if record
            % Close video file.
            close(writerObj);
        end
    end
        
    %% Plot Dynamic figure where the heading vector moves as the rotation happens
    % Assuming you have an array of quaternions "rotators" with dimensions (n, 4)
    % where n is the number of time points. Convert the quaternions to rotation matrices
%     rotMat = quat2rotm(data);

    % Define an initial vector.
    curr_pos = [1, 1, 0]; % do not use [0, 0, 0]
    
    if visual
        stopTimer = (timestamps(end)-timestamps(1))/1000; % seconds to run simulation
        framerate = 1/50; % As 1/Hz of pause for next measurement. Default to 50Hz
        if record
            % Initialize the VideoWriter object.
            writerObj = VideoWriter(fullfile(opt.FolderProcDataMat, strcat('Agent_estim_heading.avi')),'Motion JPEG AVI'); 
            open(writerObj); % Opens the file.
        end
        % Create a figure with initial, non-visible vector.
        c = quiver3(0,0,0,0,0,0,'off');
            c.Color = 'r'; % Arrow color
            c.MaxHeadSize = 2; % Arroy head size
            c.LineWidth = 2; % Arroy line width
            xlabel('X');  ylabel('Y');  zlabel('Z'); % Label axes
            xlim([-2,2]); ylim([-2,2]); zlim([-2,2]); % Fix axes scale
            daspect([1 1 1]); % Set the aspect ratio to be equal.
            title('Estimated Heading');
        
        % Timer
        ts = tic; % start timer
        pause(0.0009) % Let clock tic to a first milisecond
        % Run until elapsed time reaches set 'stopTimer' (-5 msec to avoid breaks)
        while(toc(ts) < stopTimer-0.005) 
            t = round(toc(ts)*1000); % takes the approximated msec of the run.
            % Get the rotation matrix at the current time point and rotate the previous position by the rotation matrix.
            txt = ['Time: ', num2str(t/1000), ' sec'];
%             curr_pos = rotMat(:,:,t) * curr_pos; % Get current position and rotate according to 'rotMat' step
            curr_pos = rotatepoint(data(t), curr_pos); % Get current position and rotate according to 'quaternion' step
            % Collect new datapoints 
            c.UData = curr_pos(1);
            c.VData = curr_pos(2);
            c.WData = curr_pos(3);
            
            % Update figure.
            drawnow;
                title(txt);

            % Get the frame and write it to video.
            if record
                F = getframe;           % Capture the frame
                writeVideo(writerObj,F) % add the frame to the movie
            end
            pause(framerate) % pause the run to an approx. framerate.
        end
    
        if record
            % Close video file.
            close(writerObj);
        end
    end
end

end