function [Accelerometer, Gyroscope, Magnetometer] = Deuteron_GetMotionSensors(opt)
%
%
%
%
% Version 01.03.2023 Jesus

%% Get already existing Parameters
numFiles        = length(opt.myFiles);
stream          = 2;
    
% Sort motion sensor data by data type.
% The values for acclMax and gyroMax are chosen by the user. They can be found using the Event
% File Viewer in the file started event. If not activelly changed, they
% should stay as follows:
opt.acclMax = 2*MotionSensorConstants.G; % m/s^2, max value of selected range
opt.gyroMax = 250;                       % degrees/s, max value of selected range
opt.magMax  = MotionSensorConstants.Magnetometer9250Range; % Teslas, max value of selected range

% Create structs
Accelerometer   = struct('X',[],'Y',[],'Z',[],'t',[],'max',[]);
Gyroscope       = struct('X',[],'Y',[],'Z',[],'t',[],'max',[]);
Magnetometer    = struct('X',[],'Y',[],'Z',[],'t',[],'max',[]);
            
% Axes description. With board plugged on animal's head, and according to
% the sensor sheet:
%   Magnetometer: X for vertical, Y horizontal AP and Z horiz DL.
%   Acc and Gyro: X for horiz AP, Y for vertical, and Z horiz DL.
%
% By testing it, appears that:
%   X: Dorso-Ventral, vertical. Shows Gravity!.
%   Y: Antero-Posterior, horizontal forward-backward
%   Z: Medio-Lateral, horizontal sideward
%
% Feel like the X-Y axes for Acc and Gyro are rotated, to match the
% Magnetometer? Ask jacob.

% Therefore:
%     Head movement is detected mostly by accelerometer Y,Z. (being -Y forward, +Y backward?)
%     The gyroscope X would determine head turning. (+X right, -X left?)
%     The magnetometer needs the two horizontal axes (planar axes) for head orientation.
% Being X gravity, we need Y and Z from Magnetometer.

for i = 1:numFiles         
    if ~strcmp(opt.myFiles(i).name(1:4),'NEUR')
        % Skips Event opt.myFiles (do not contain data)
        continue
    else
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name), 'r');
        data = Deuteron_extractData(stream, fid, opt);
        fclose(fid);

        Accelerometer.X   = [Accelerometer.X data.Accelerometer.Data.X'];
        Accelerometer.Y   = [Accelerometer.Y data.Accelerometer.Data.Y'];
        Accelerometer.Z   = [Accelerometer.Z data.Accelerometer.Data.Z'];
        Accelerometer.t   = [Accelerometer.t data.Accelerometer.timestamps];
        Accelerometer.max = opt.acclMax;

        Gyroscope.X       = [Gyroscope.X data.Gyroscope.Data.X'];
        Gyroscope.Y       = [Gyroscope.Y data.Gyroscope.Data.Y'];
        Gyroscope.Z       = [Gyroscope.Z data.Gyroscope.Data.Z'];
        Gyroscope.t       = [Gyroscope.t data.Gyroscope.timestamps]; 
        Gyroscope.max     = opt.gyroMax;

        Magnetometer.X    = [Magnetometer.X data.Magnetometer.Data.X'];
        Magnetometer.Y    = [Magnetometer.Y data.Magnetometer.Data.Y'];
        Magnetometer.Z    = [Magnetometer.Z data.Magnetometer.Data.Z'];
        Magnetometer.t    = [Magnetometer.t data.Magnetometer.timestamps];
        Magnetometer.max  = opt.magMax;
    end
end

%% Save data to matfile
save("MotionData.mat", "Magnetometer", "Gyroscope", "Accelerometer", '-mat');

end
