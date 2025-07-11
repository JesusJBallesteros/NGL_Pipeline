function Deuteron_GetMotionSensors(opt)
% This function allows for Motion Data extraction from Deuteron's raw DT1 data
% format. Uses a few fixed parameters to give proper units to the extracted
% timeseries and sorts each sensor's data to its respective variable.
% It plots the extracted raw data and saves it into separated file.
% It processes the data using an attitude and heading reference system
% (AHRS) to hopefully put the data in a meaningful reference system that
% can be used to predict/estimate the animal's position/heading.
% 
% Last: Added warning for deadtime
% Jesus 22.04.2025

%% Some local Parameters
numFiles        = length(opt.myFiles);
opt.stream      = 2;
  
% Bochum Magnetic Field Horizontal Intensity. According to 
% https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
MField_Bochum = 19.7; % uTesla.

% Sample rate for motion sensors is 1000Hz. 
fsmot          = 1000;         % Sample Rate of the feeded data (Hz)

% Gyro/Accel_Noise are determined from the hardware datasheets.
%Gyro_Noise  = 1.7453e-04;   % Gyroscope Noise (variance value) in units of rad/s. (MPU-9250: 0.01 deg/sec)
%Accel_Noise = .008;        % Accelerometer Noise(variance value) in units of m/s^2 (g). (MPU-9250: 8 mg)

% The values for acclMax and gyroMax are chosen by the user. They can be found using the Event
% File Viewer in the file started event. If not activelly changed, they should stay as follows:
opt.acclMax = 2*MotionSensorConstants.G; % m/s^2, max value of selected range
opt.gyroMax = 250;                       % degrees/s, max value of selected range
opt.magMax  = MotionSensorConstants.Magnetometer9250Range; % Teslas, max value of selected range
    
if ~isfile(fullfile(opt.FolderProcDataMat, strcat('MotionData.mat')))
    %% Sort motion sensor data by data type.
    % Create structs
    data.acc    = struct('X', [], 'Y', [], 'Z', [], 'max', opt.acclMax);
    data.gyr    = struct('X', [], 'Y', [], 'Z', [], 'max', opt.gyroMax);
    data.mag    = struct('X', [], 'Y', [], 'Z', [], 'max', opt.magMax);
    timestamps  = [];
                
    % Axes description. With board plugged on animal's head, and according to the sensor sheet:
    %   Magnetometer: Y for vertical, X for AP and Z for DL.
    %   Acc: +X for vertical up, +Y for AP forward and +Z for DL left.
    %   Gyro: around +X for yaw left (look around),
    %         around +Y for roll left (rolling, 'croqueta'),
    %         around +Z for pitch down (nodding).
    % note: Feels like X-Y axes in Acc and Gyro are interchanged with the ones in the Magnetometer
    % Therefore:
    %     Head movement is detected mostly by gyroscope yaw (X)
    %     The magnetometer needs the two horizontal (planar) axes for head orientation: Z and Y
    %     Acc.X detects gravity acceleration (points down-up axes!)
    for i = 1:numFiles         
        if ~strcmp(opt.myFiles(i).name(1:4),'NEUR')
            % Skips Event and other files (do not contain data)
            continue
        else
            fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name), 'r');
                mot = Deuteron_extractData(fid, opt);
            fclose(fid);
    
            % MPU-9250 has two devices, the magnetometer and the accelerometer-gyroscope, on the same board. 
            % The axes of these devices are different from each other. The magnetometer axis is aligned 
            % with the NED coordinates. The axis of the accelerometer-gyroscope is different from 
            % magnetometer in MPU-9250. The accelerometer and the gyroscope axis need to be swapped and/or 
            % inverted to match the magnetometer axis. For more information refer to the section 
            % "Orientation of Axes" section in MPU-9250 datasheet.        
            data.acc.X = [data.acc.X  mot.Accelerometer.Data.X']; 
            data.acc.Y = [data.acc.Y  mot.Accelerometer.Data.Y']; 
            data.acc.Z = [data.acc.Z  -mot.Accelerometer.Data.Z']; % Note the sign inversion, to match magnetic
    
            data.gyr.X = [data.gyr.X  -mot.Gyroscope.Data.X'];
            data.gyr.Y = [data.gyr.Y  -mot.Gyroscope.Data.Y'];
            data.gyr.Z = [data.gyr.Z  -mot.Gyroscope.Data.Z']; % Note the sign inversion, to match magnetic
    
            data.mag.X = [data.mag.X mot.Magnetometer.Data.X']; 
            data.mag.Y = [data.mag.Y mot.Magnetometer.Data.Y']; 
            data.mag.Z = [data.mag.Z mot.Magnetometer.Data.Z'];
            
            timestamps = [timestamps mot.Accelerometer.timestamps];
        end
    end
    
    % % Remove all timestamps where all readings are 0 (failsafe)
    idx0 = find(~data.mag.X & ~data.mag.Y & ~data.mag.Z & ...
                ~data.acc.X & ~data.acc.Y & ~data.acc.Z & ...
                ~data.gyr.X & ~data.gyr.Y & ~data.gyr.Z);
    if ~isempty(idx0)
        data.acc.X(idx0) = [];
        data.acc.Y(idx0) = [];
        data.acc.Z(idx0) = [];
        data.gyr.X(idx0) = [];
        data.gyr.Y(idx0) = [];
        data.gyr.Z(idx0) = [];
        data.mag.X(idx0) = [];
        data.mag.Y(idx0) = [];
        data.mag.Z(idx0) = [];
        timestamps(idx0) = [];
    end
    
    % Find deadtimes
    idxt = find(diff(timestamps)>3);
    if ~isempty(idxt)
        timestamps(idxt+1:end) = timestamps(idxt+1:end)-(timestamps(idxt+1)-timestamps(idxt)-1);
    end
    
    % Calculate timestamps in seconds
    tsec = timestamps/16000 + 12/fsmot; % since midnight
    tsec = tsec - tsec(1); % relativize to recording

    %% Plot sensors readings. RAW.
    % Run the plot function. 
    Deuteron_PlotMotionSensors(data, tsec, opt, 1); % 4th input == 1 (raw data)
    exportgraphics(gcf, fullfile(opt.FolderProcDataMat, strcat('motion_raw.png')), 'Resolution', 300)
    close gcf

    clear fid mot numFiles idx0 timestamps

    %% Use the magcal function to obtain the correction coefficients for the
    % magnetometer. This helps with the typical soft/hard iron effect on the
    % magnetic field (distortion of the magnetic sphere).
    
    % First we feed the Magnetometer matrix as raw, to obtain the coefficients.
    % Put the readings as the scripts like them (t x axis matrices). Use filtered time series
    MS_Mag = [data.mag.X; data.mag.Y; data.mag.Z]';
    
    % Magcal
    [A, b, Mfield] = magcal(MS_Mag);  % A = 3x3 matrix for soft iron correction 
                                      % b = 3x1 vector for hard iron correction
    
    % Display measured and expected Magnetic field in uTesla. Only informative.
    disp(['The magnetic field for Bochum should be ~', int2str(MField_Bochum), ' uTesla. ' ...
        'The measured field was ~', int2str(Mfield*1000000), ' uTesla.']);
    
    % Then, apply the corrections and re-create the Magnetometer matrix again,
    % using the correction coefficients obtained from magcal.
    MS_Mag = [data.mag.X-b(1); data.mag.Y-b(2); data.mag.Z-b(3)]' * A;
    
    % Put data back to a struct for later.
    data.magcorr.X = MS_Mag(:,1)';
    data.magcorr.Y = MS_Mag(:,2)';
    data.magcorr.Z = MS_Mag(:,3)';
    
    % Plot sensors readings. MAG CORRECTED.
    % Run the plot function, set 3rd input to 1 to plot corrected data.
    Deuteron_PlotMotionSensors(data, tsec, opt, 1);
    exportgraphics(gcf, fullfile(opt.FolderProcDataMat, strcat('motion_magcorr.png')), 'Resolution', 300)
    close gcf

    clear A b i

    %% Save data to matfile
    save(fullfile(opt.FolderProcDataMat, strcat('MotionData.mat')), "data", "tsec", '-mat');
    disp('Magnetic-Corrected motion data saved.')
else
    load(fullfile(opt.FolderProcDataMat, strcat('MotionData.mat')));
end

% %% Use ecompass to merge Acc and Mag only, for the first 500 samples. Average to
% % % get an estimate of initial heading. Output is in quaternions as 'rotators.ecomp'
% % rotators.ecomp = ecompass(MS_Acc(1:500,:), MS_Mag(1:500,:), "rotmat");
% % rotators.ecomp = mean(rotators.ecomp,3); 
% % % poseplot(rotators.ecomp(1:3,:)); % to check
% MS_Acc = [data.acc.X; data.acc.Y; data.acc.Z]';
% MS_Mag = [data.mag.X; data.mag.Y; data.mag.Z]';
% 
% % use ecompass to fuse acc and mag (corrected)
% rotators = ecompass(MS_Acc, MS_Mag, 'quaternion');
% 
% % The slerp function is used to steer the filter state towards the current input. 
% % It is steered more towards the input when the difference between the input and current
% % filter state has a large dist, and less toward the input when dist gives a small value.
% % The interpolation parameter to slerp is in the closed-interval [0,1], so the output
% % of dist must be re-normalized to this range. However, the full range of [0,1] for the
% % interpolation parameter gives poor performance, so it is limited to a smaller range
% % hrange centered at hbias.
% slerpf.hrange = 0.2;
% slerpf.hbias = 0.4;
% 
% % Limit low and high to the interval [0, 1].
% slerpf.low  = max(min(slerpf.hbias - (slerpf.hrange./2), 1), 0);
% slerpf.high = max(min(slerpf.hbias + (slerpf.hrange./2), 1), 0);
% slerpf.hrangeLimited = slerpf.high - slerpf.low;
% 
% % Initialize the filter and preallocate outputs.
% y = rotators(1); % initial filter state
% rot_filt = zeros(size(y), 'like', y); % preallocate filter output
% rot_filt(1) = y;
% 
% % Filter the noisy trajectory, sample-by-sample.
% for ii=2:numel(rotators)
%     x = rotators(ii);
%     d = dist(y, rotators(ii));
% 
%     % Renormalize dist output to the range [low, high]
%     hlpf = (d./pi).*slerpf.hrangeLimited + slerpf.low;
%     y = slerp(y, x, hlpf);
%     rot_filt(ii) = y;
% end
% clear x d y hlpf
% 
% %% Plot
% % Run the plot function.
% Deuteron_PlotMotionSensors(rot_filt, tsec, opt, 2, 1, 1); % rotators data
end
