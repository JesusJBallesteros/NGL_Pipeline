function [Accelerometer, Gyroscope, Magnetometer, param] = Deuteron_GetMotionSensors(in, sessions, ss)
%
%
%
%
%

%% Get already existing Parameters
Files           = sessions.info{ss}.files;
numFiles        = length(sessions.info{ss}.files);
stream          = 2;
    
% Sort motion sensor data by data type.
% The values for acclMax and gyroMax are chosen by the user. They can be found using the Event
% File Viewer in the file started event.
param.acclMax = 2*MotionSensorConstants.G; % m/s^2, max value of selected range
param.gyroMax = 250;                       % degrees/s, max value of selected range
param.magMax  = MotionSensorConstants.Magnetometer9250Range; % Teslas, max value of selected range

Accelerometer   = struct('X',[],'Y',[],'Z',[],'t',[]);
Gyroscope       = struct('X',[],'Y',[],'Z',[],'t',[]);
Magnetometer    = struct('X',[],'Y',[],'Z',[],'t',[]);
            
for i = 1:numFiles         
    if ~strcmp(Files(i).name(1:4),'NEUR')
        % Skips Event files (do not contain data)
        continue
    else
        fid = fopen(fullfile(in.pathRaw, Files(i).name), 'r');
        data = Deuteron_extractData(stream, fid, param);
        fclose(fid);

        Accelerometer.X   = [Accelerometer.X data.Accelerometer.Data.X'];
        Accelerometer.Y   = [Accelerometer.Y data.Accelerometer.Data.Y'];
        Accelerometer.Z   = [Accelerometer.Z data.Accelerometer.Data.Z'];
        Accelerometer.t   = [Accelerometer.t data.Accelerometer.timestamps];

        Gyroscope.X       = [Gyroscope.X data.Gyroscope.Data.X'];
        Gyroscope.Y       = [Gyroscope.Y data.Gyroscope.Data.Y'];
        Gyroscope.Z       = [Gyroscope.Z data.Gyroscope.Data.Z'];
        Gyroscope.t       = [Gyroscope.t data.Gyroscope.timestamps]; 

        Magnetometer.X    = [Magnetometer.X data.Magnetometer.Data.X'];
        Magnetometer.Y    = [Magnetometer.Y data.Magnetometer.Data.Y'];
        Magnetometer.Z    = [Magnetometer.Z data.Magnetometer.Data.Z'];
        Magnetometer.t    = [Magnetometer.t data.Magnetometer.timestamps];
    end
end

end
