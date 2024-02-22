function Deuteron_GetMotionSensors(opt)
% This function allows for Motion Data extraction from Deuteron's raw DT1 data
% format. Uses a few fixed parameters to give proper units to the extracted
% timeseries and sorts each sensor's data to its respective variable.
% It plots the extracted raw data and saves it into separated file.
% It processes the data using an attitude and heading reference system
% (AHRS) to hopefully put the data in a meaningful reference system that
% can be used to predict/estimate the animal's position/heading.
% 
% IN WORKING PROCESS
%
% Jesus 26.01.2024

%% Some local Parameters
numFiles        = length(opt.myFiles);
opt.stream      = 2;
  
% Bochum Magnetic Field Horizontal Intensity. According to 
% https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
MField_Bochum = 19.7; % uTesla.

% Sample rate for motion sensors is 1000Hz. 
fs          = 1000;         % Sample Rate of the feeded data (Hz)

% Gyro/Accel_Noise are determined from the hardware datasheets.
Gyro_Noise  = 1.7453e-04;   % Gyroscope Noise (variance value) in units of rad/s. (MPU-9250: 0.01 deg/sec)
Accel_Noise = .008;        % Accelerometer Noise(variance value) in units of m/s^2 (g). (MPU-9250: 8 mg)

% The values for acclMax and gyroMax are chosen by the user. They can be found using the Event
% File Viewer in the file started event. If not activelly changed, they should stay as follows:
opt.acclMax = 2*MotionSensorConstants.G; % m/s^2, max value of selected range
opt.gyroMax = 250;                       % degrees/s, max value of selected range
opt.magMax  = MotionSensorConstants.Magnetometer9250Range; % Teslas, max value of selected range

%% Sort motion sensor data by data type.
% Create structs
Accelerometer   = struct('X', [], 'Y', [], 'Z', [], 'max', opt.acclMax);
Gyroscope       = struct('X', [], 'Y', [], 'Z', [], 'max', opt.gyroMax);
Magnetometer    = struct('X', [], 'Y', [], 'Z', [], 'max', opt.magMax);
timestamps      = [];
            
% Axes description. With board plugged on animal's head, and according to the sensor sheet:
%   Magnetometer: Y for vertical, X for AP and Z for DL.
%   Acc: +X for vertical up, +Y for AP forward and +Z for DL left.
%   Gyro: around +X for yaw left (look around),
%         around +Y for roll left (rolling, 'croqueta'),
%         around +Z for pitch down (nodding).
% note: Feels like X-Y axes in Acc and Gyro are interchanged with the ones in the Magnetometer
%
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
         data = Deuteron_extractData(fid, opt);
        fclose(fid);

        Accelerometer.X   = [Accelerometer.X  data.Accelerometer.Data.Y']; % Note here comes Y axis, to make it X.
        Accelerometer.Y   = [Accelerometer.Y  data.Accelerometer.Data.X']; % Note here comes X axis, to make it Y.
        Accelerometer.Z   = [Accelerometer.Z -data.Accelerometer.Data.Z']; % Note the sign inversion, to match magnetic

        Gyroscope.X       = [Gyroscope.X  data.Gyroscope.Data.Y']; % Note here comes Y axis, to make it X.
        Gyroscope.Y       = [Gyroscope.Y  data.Gyroscope.Data.X']; % Note here comes X axis, to make it Y.
        Gyroscope.Z       = [Gyroscope.Z -data.Gyroscope.Data.Z']; % Note the sign inversion, to match magnetic

        Magnetometer.X    = [Magnetometer.X data.Magnetometer.Data.X']; % Note this stays, so above match magnetic X.
        Magnetometer.Y    = [Magnetometer.Y data.Magnetometer.Data.Y']; % Note this stays, so above match magnetic Y.
        Magnetometer.Z    = [Magnetometer.Z data.Magnetometer.Data.Z']; % Note no sign inversion.
        
        timestamps        = [timestamps data.Accelerometer.timestamps];
    end
end

clear data fid

%% Save data to matfile
save(fullfile(opt.FolderProcDataMat, strcat('MotionData_raw.mat')), "Magnetometer", "Gyroscope", "Accelerometer", "timestamps", '-mat');
disp('Raw motion data saved.')

%% Plot sensors readings. UNTREATED, UNCORRECTED.
% pass data into single variable
data.acc = Accelerometer;
data.gyr = Gyroscope;
data.mag = Magnetometer;

% Run the plot function. 3r input == 1 (only raw data)
Deuteron_PlotMotionSensors(data, opt, timestamps, 1);
exportgraphics(gcf, fullfile(opt.FolderProcDataMat, strcat('MotionRaw.png')),'Resolution',300)
close gcf
clear data

%% Transform readings to match real axes
% An attitude and heading reference system (AHRS) consist of a 9-axis system 
% that uses an accelerometer, gyroscope, and magnetometer to compute orientation 
% of the device. The 'ahrsfilter' produces a smoothly changing estimate of 
% orientation of the device, while correctly estimating the north direction. 
% The 'ahrsfilter' has the ability to remove gyroscope bias and can also detect 
% and reject mild magnetic jamming.
% The following code snippets use 'ahrsfilter' system object to determine 
% orientation of the sensor and creates a figure which gets updated as you 
% move the sensor. The sensor has to be stationary, before the start of this example.

% Put the readings as the scripts like them (t x axis matrices).
% Note the axes swapping in Accelerometer and Gyroscope, to match NED magnetometer 
% coordinates system.
MS_Mag = [ Magnetometer.X;   Magnetometer.Y;   Magnetometer.Z]';
MS_Acc = [-Accelerometer.X; -Accelerometer.Z;  Accelerometer.Y]'; % Z/Y swapped
MS_Gyr = [ Gyroscope.X;      Gyroscope.Z;     -Gyroscope.Y]';     % Z/Y swapped

%% Use the magcal function to obtain the correction coefficients for the
% magnetometer. This helps with the typical soft/hard iron effect on the
% magnetic field (distortion of the magnetic sphere).
% First we feed the Magnetometer matrix as raw, to obtain the coefficients.
[A, b, Mfield] = magcal(MS_Mag);  % A = 3x3 matrix for soft iron correction 
                                  % b = 3x1 vector for hard iron correction

% Display measured and expected Magnetic field in uTesla. Only informative.
disp(['The magnetic field for Bochum should be ~', int2str(MField_Bochum), ' uTesla. ' ...
    'The measured field was ~', int2str(Mfield*1000000), ' uTesla.']);

% Then, apply the corrections and re-create the Magnetometer matrix again.
MS_Mag = [Magnetometer.X-b(1); Magnetometer.Y-b(2); Magnetometer.Z-b(3)]' * A;

% Bring back to a struc just for comparison 
Magnetometer.X = MS_Mag(:,1)';
Magnetometer.Y = MS_Mag(:,2)';
Magnetometer.Z = MS_Mag(:,3)';

% pass data into single variable
data.acc = Accelerometer;
data.gyr = Gyroscope;
data.mag = Magnetometer;

%% Plot sensors readings. TREATED, CORRECTED.
% Run the plot function, set 3rd input to 1 to plot corrected data.
Deuteron_PlotMotionSensors(data, opt, timestamps, 1); % corrected data
exportgraphics(gcf, fullfile(opt.FolderProcDataMat, strcat('MotionRaw_corrected.png')),'Resolution',300)
close gcf

clear A b Accelerometer Gyroscope Magnetometer

%% Create AHRS filter using matlab tools. TO VERIFY.
% Needs the sample rate and the sensor noise levels. The output once the 
% FUSE object is applied will be in 'quaternions' a complex expression of 
% 3D rotations. Not sure if themost useful, but let's see. Could be substituted 
% by 'Rotation matrices'which I think is the translation of the quaternions as 
% point rotations in a 3D coordinate system.

% FUSE = ahrsfilter('SampleRate',                     fs, ...
%                   'DecimationFactor',               1, ...
%                   'AccelerometerNoise',             Accel_Noise, ...
%                   'GyroscopeNoise',                 Gyro_Noise, ...
%                   'ExpectedMagneticFieldStrength',  Mfield*1000000, ...
%                   'OrientationFormat',              'quaternion');
% 
% % Run the thing for every timepoint. Output is in quaternions
% [orientation, ~] = FUSE(MS_Acc,MS_Gyr,MS_Mag);
% 
% clear Gyro_Noise Accel_Noise MField_Bochum

%% Plot
% Run the plot function, set 3rd input to 2 to plot FUSE data.
% Deuteron_PlotMotionSensors(orientation, opt, timestamps, 2, 1, 1); % orientation data
% clear data
end
