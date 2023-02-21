function Deuteron_PlotMotionSensors(Acc, Gyr, Mag, param)
%
%
%
%
%

% Axes description. All three sensors share the schema. 
% With board plugged on animal's head:
%     X: Dorso-Ventral, vertical. Gravity.
%     Y: Antero-Posterior, horizontal forward-backward
%     Z: Medio-Lateral, horizontal sideward
% Therefore:
%     Head movement is detected mostly by accelerometer Y,Z. (being -Y forward, +Y backward?)
%     The gyroscope X would determine head turning. (+X right, -X left?)
%     The magnetometer needs the two horizontal axes (planar axes) for head orientation.
% Being X gravity, we need Y and Z from Magnetometer.

% TODO
calib_check();

% This formula gives the direction of the [Z, Y] vector, counted clockwise from the Y axis.
% Since Y axis is AP, it gives degrees from the 0 angle (forward vector)
az = atan2(Mag.Y, Mag.Z) * 180/pi;

% Plots 
figure,
if ~isempty(Acc)
    MSData = Acc;

    subplot(3,1,1)
    title('Accelerometer');
    % General plot
    plot(MSData.t, MSData.X); hold on
    plot(MSData.t, MSData.Y); hold on
    plot(MSData.t, MSData.Z);
    ylim([-param.acclMax param.acclMax]); ylabel('m/s^2');
    xlabel('ms');
    legend({'X' 'Y' 'Z'}, 'Box','off');
     box("off")
end

if ~isempty(Gyr)
    MSData = Gyr;

    subplot(3,1,2)
    title('Gyroscope');
    % General plot
    plot(MSData.t, MSData.X); hold on
    plot(MSData.t, MSData.Y); hold on
    plot(MSData.t, MSData.Z);
    ylim([-param.gyroMax param.gyroMax]); ylabel('deg/s')
    xlabel('ms');
    legend({'X' 'Y' 'Z'}, 'Box','off');
    box("off")
end

if ~isempty(Mag)
    MSData = Mag;

    subplot(3,1,3)
    title('Magnetometer');
    % General plot
    plot(MSData.t, MSData.X); hold on
    plot(MSData.t, MSData.Y); hold on
    plot(MSData.t, MSData.Z);
%         ylim([-param.magMax param.magMax]);
    ylim([-1e-4 1e-4]);  ylabel('Tesla')
    xlabel('ms');
    legend({'X' 'Y' 'Z'}, 'Box','off');
    box("off")
end

end

function calib = calib_check()
%% Magnetometer details and calibration
% Planar axes magnetic field intensity readings. This better be in a STATIC
% situation, parallel to surface. An average of a period of time probably
% useful.
Bx = mean(Magnetometer.Y(100:601)); % uT
By = mean(Magnetometer.Z(100:601)); % uT

% BochumMagnetic Field Horizontal Intensity:
MFhi = 19.7; % uTesla. Acc. to https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm

% The measured MFhi by the sensor can be calculated as the square root of
% the sum of the intensity readings from the planar axes, squared. Or:
Bh = sqrt(Bx^2 + By^2);

% Result TODO
if Bh < MFhi, calib = 1;
else,         calib = 0; end

end