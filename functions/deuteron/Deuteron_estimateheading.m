function [heading, position] = Deuteron_estimateheading(Accelerometer, Gyroscope, Magnetometer)
% Deuteron_estimateheading  Estimate animal heading and position from MPU-9250 data.
%
% PURPOSE:
%   Uses an AHRS (Attitude and Heading Reference System) filter to fuse
%   accelerometer, gyroscope, and magnetometer readings from the Deuteron
%   MPU-9250 sensor into a heading estimate (orientation quaternion → rotation
%   vector) and a dead-reckoning position estimate.
%   Magnetometer hard/soft-iron distortion is corrected via MATLAB's magcal
%   before AHRS fusion.
%
% USAGE:
%   [heading, position] = Deuteron_estimateheading(Accelerometer, Gyroscope, Magnetometer)
%   Called from getfrom_Deuteron inside GetMotionSensors.
%
% INPUTS:
%   Accelerometer  - struct with fields .X, .Y, .Z (m/s²) and .t (timestamps)
%   Gyroscope      - struct with fields .X, .Y, .Z (deg/s)
%   Magnetometer   - struct with fields .X, .Y, .Z (Tesla)
%   All sensor arrays must be the same length (1000 Hz sample rate assumed).
%
% OUTPUTS:
%   heading   - [N × 3 double] rotation vector in radians at each time step
%   position  - [N × 2 double] dead-reckoning [x, y] displacement in meters
%
% NOTES:
%   - Sensor axes are remapped to NED-like convention before AHRS fusion:
%       gyr = [X, Z, −Y], acc = [−X, −Z, Y]
%   - magcal corrects for hard-iron (offset b) and soft-iron (matrix A) effects.
%   - Noise parameters (Gyro_Noise, Accel_Noise) are MPU-9250 datasheet values.
%   - Dead-reckoning from acceleration is coarse; treat position as qualitative.
%
% CALLS:
%   magcal, ahrsfilter, rotvec (Sensor Fusion and Tracking Toolbox)
%
% Last modified 08.05.2026 (Jesus)

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
        
        if r