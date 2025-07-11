function [heading, position] = Deuteron_estimateheading(Accelerometer, Gyroscope, Magnetometer)
%
%
%
%
%
%

    % sample rate for motion sensors is 1000Hz. GyroscopeNoise and AccelerometerNoise
    % are determined from the hardware specifications .
    %  TODO, find out the real noise of our sensors, if different model.
    fs          = 1000;         % Sample Rate of the feeded data (Hz)
    Gyro_Noise  = 3.0462e-06;   % Gyroscope Noise (variance value) in units of rad/s. (MPU-9250)
    Accel_Noise = 0.0061;       % Accelerometer Noise(variance value)in units of m/s^2. (MPU-9250)
    
    Magnet = [ Magnetometer.X;   Magnetometer.Y;   Magnetometer.Z]';
    gyr = [Gyroscope.X',       Gyroscope.Z',     -Gyroscope.Y'];
    acc = [-Accelerometer.X', -Accelerometer.Z', Accelerometer.Y'];

    % Use the magcal function to obtain the correction coefficients for the
    % magnetometer. This helps with the typical soft/hard iron effect on the
    % magnetic field (distortion of the magnetic sphere).
    % First we feed the Magnetometer matrix as raw, to obtain the coefficients.
    [A, b, Mfield] = magcal(Magnet);  % A = 3x3 matrix for soft iron correction 
                                   % b = 3x1 vector for hard iron correction
    
    % Bochum Magnetic Field Horizontal Intensity. According to 
    % https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
    MField_Bochum = 19.7; % uTesla.
    
    % Display measured and expected Magnetic field in uTesla. Only informative
    % of not big difference to what it would be expected.
    disp(['The magnetic field for Bochum is ~', int2str(MField_Bochum), ...
        ' uTesla. The measured field seems to be ~', int2str(Mfield*1000000), ' uTesla.']);
    
    % Then we apply the corrections and create the Magnetometer matrix again.
    mag = [Magnetometer.X-b(1);  Magnetometer.Y-b(2);  Magnetometer.Z-b(3)]' * A;
    clear Magnet b A

    % AHRS filter initialization
    ahrs_filter = ahrsfilter('SampleRate', fs, 'GyroscopeNoise', Gyro_Noise, 'AccelerometerNoise', Accel_Noise);    

    % Run AHRS filter on sensor data
    orientation = ahrs_filter(gyr, acc, mag);
    
    % Calculate heading direction in radians
    heading = rotvec(orientation);
    
    % Calculate position using dead reckoning
    dt = 1/fs; % Sampling interval
    v = sqrt(Accelerometer.X.^2 + Accelerometer.Y.^2 + Accelerometer.Z.^2); % Total acceleration
    d = cumsum(v*dt); % Displacement
    dx = d'.*cos(heading); % Displacement in x direction
    dy = d'.*sin(heading); % Displacement in y direction
    

    %% Generate dynamic/static Agent position estimation
    if visual
        stopTimer = (Accelerometer.t(end)-Accelerometer.t(1))/fs; % seconds to run simulation
        framerate = 1/50; % As 1/Hz of pause for next measurement. Default to 50Hz
        
        if record
            % initialize the VideoWriter object.
            writerObj = VideoWriter('Agent_position_estim.avi','Motion JPEG AVI'); 
            open(writerObj); % Opens the file.
        end
    end

    % Plot position on cenital view of environment
    quiver(0,0,0,0);
    % Starting point
    axis equal;
    xlim([-1,1]);
    ylim([-1,1]);
    title('Agent Position');
    xlabel('X (m)');
    ylabel('Y (m)');
    hold on

    if visual
        % Prepare for non visualization purpouses TODO 
        ts = tic; % start timer
        pause(0.001) % Let clock tic to a first milisecond
        
        % run until elapsed time reaches set 'stopTimer' minus 1 msec to avoid breaking it
        while(toc(ts) < stopTimer-0.001) 
            t = round(toc(ts)*1000); % takes the approximated milisecond of the run.
            position = [sum(dx(t,:)), sum(dy(t,:))];
            
            h.XData = position(1); % Current position
            h.YData = position(2); % Current position
            drawnow;

            % If desired, get the frame and write it to video.
            if record
                F = getframe;           % Capture the frame
                writeVideo(writerObj,F) % add the frame to the movie
            end

        pause(framerate) % pause the run for 1/Hz msec, to an approx framerate.
        end

        if record
            % Close video file.
            close(writerObj);
        end
    else
%         % Run the thing for every single timepoint
%         rotators = FUSE(Acc,Gyr,Mag);
    end

end