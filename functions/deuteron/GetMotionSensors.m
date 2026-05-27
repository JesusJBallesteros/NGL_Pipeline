function GetMotionSensors(opt, input)
% GetMotionSensors  Extract raw IMU / accelerometer data into a common .mat.
%
% PURPOSE:
%   Dispatches raw motion-sensor extraction based on recording format:
%   - INTAN (fileperch): reads AUX*.dat ADXL335 accelerometer channels
%   - Deuteron (DT2/DF1): reads NEUR*.DF1 stream 2 via Deuteron_extractData
%   Converts ADC values to physical units, removes dead samples, corrects for
%   timestamp gaps, and saves MotionData_raw.mat.
%
%   NO axis corrections, filtering, AHRS, magcal, plotting or event-locked
%   analysis are applied here — those all live in ProcessMotionSensors (moved 
%   to NGL02).
%
% USAGE:
%   GetMotionSensors(opt, input)
%   Called from INTAN_PipelineWrapper and Deuteron_PipelineWrapper
%   when opt.GetMotionSensors = true.
%
% INPUTS:
%   opt    - options struct; relevant fields:
%              .FolderProcDataMat  output folder
%              .PathRaw            raw data folder
%              .myFiles            dir-struct of session files (Deuteron)
%   input  - struct; relevant field:
%              .sessions(x).info.fileformat  'fileperch' | 'DT2' | 'DF1'
%
% OUTPUT:
%   MotionData_raw.mat in opt.FolderProcDataMat.
%   Contains struct `raw` with fields:
%     .acc.X/Y/Z   [1×N] m/s²   accelerometer,  native chip frame
%     .gyr.X/Y/Z   [1×N] deg/s  gyroscope,      native chip frame ([] INTAN)
%     .mag.X/Y/Z   [1×N] µT     magnetometer,   native chip frame ([] INTAN)
%     .tsec        [1×N] s      time from session start
%     .fs          Hz           sample rate after any decimation
%     .source      char         'deuteron-mpu9250' | 'intan-adxl335'
%     .dof         scalar       9 (Deuteron) | 3 (INTAN)
%     .chipFrame   struct       documents native axis convention
%
% CALLS:
%   Deuteron_extractData, MotionSensorConstants
%
% SEE ALSO:
%   ProcessMotionSensors (NGL02 analysis — axis alignment, AHRS, plotting)
%
% Last modified 27.05.2026 (Jesus)

rawMatPath = fullfile(opt.FolderProcDataMat, 'MotionData_raw.mat');
    if isfile(rawMatPath)
        disp('- Raw motion data already exists for this session. Skipping extraction.')
        return
    end

fmt = input.sessions(input.run(1)).info.fileformat;
    if strcmp(fmt, 'fileperch')
        getfrom_INTAN(opt)

    elseif strcmp(fmt, 'DT2') || strcmp(fmt, 'DF1')
        getfrom_Deuteron(opt)

    else
        warning('NGL:unknownFormat', ...
            'GetMotionSensors: unrecognised format ''%s''. Skipping.', fmt);
    end

end

% Actual Functions
function getfrom_Deuteron(opt)
    % Read 9-DoF MPU-9250 data (acc + gyro + mag) from Deuteron DF1 files.
    % Values are in physical units as delivered by Deuteron_extractData /
    % ScaleMotionSensorData. Native chip frame is preserved — no remapping.
    
    numFiles   = length(opt.myFiles);
    opt.stream = 2;   % motion-sensor stream index inside DF1
    
    % Hardware scaling limits (documented in Event File Viewer, "started" event)
    opt.acclMax = 2 * MotionSensorConstants.G;                    % m/s²
    opt.gyroMax = 250;                                            % deg/s
    opt.magMax  = MotionSensorConstants.Magnetometer9250Range;    % T → converted to µT below
    
    accX = {}; accY = {}; accZ = {};
    gyrX = {}; gyrY = {}; gyrZ = {};
    magX = {}; magY = {}; magZ = {};
    ts   = {};
    
    disp('- Extracting Deuteron motion sensor data (stream 2)...')
    for i = 1:numFiles
        if strlength(opt.myFiles(i).name) < 4 || ~startsWith(opt.myFiles(i).name, 'NEUR')
            continue   % skip event files and non-neural files
        end
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name), 'r');
            mot = Deuteron_extractData(fid, opt);
        fclose(fid);
    
        % Physical units are already applied by ScaleMotionSensorData inside
        % Deuteron_extractData: acc in m/s², gyro in deg/s, mag in T.
        % Store in native chip frame — axis remapping is ProcessMotionSensors' job.
        % MPU-9250 note: acc/gyro axes DIFFER from magnetometer axes (see datasheet
        % "Orientation of Axes"). Both are saved as-is; document this in chipFrame.
        accX{i} = mot.Accelerometer.Data.X';  %#ok<*AGROW>
        accY{i} = mot.Accelerometer.Data.Y';
        accZ{i} = mot.Accelerometer.Data.Z';
    
        gyrX{i} = mot.Gyroscope.Data.X';
        gyrY{i} = mot.Gyroscope.Data.Y';
        gyrZ{i} = mot.Gyroscope.Data.Z';
    
        magX{i} = mot.Magnetometer.Data.X';
        magY{i} = mot.Magnetometer.Data.Y';
        magZ{i} = mot.Magnetometer.Data.Z';
    
        ts{i}   = mot.Accelerometer.timestamps;
    end
    
    % Concatenate all files
    raw.acc.X = cat(2, accX{:});
    raw.acc.Y = cat(2, accY{:});
    raw.acc.Z = cat(2, accZ{:});
    
    raw.gyr.X = cat(2, gyrX{:});
    raw.gyr.Y = cat(2, gyrY{:});
    raw.gyr.Z = cat(2, gyrZ{:});
    
    raw.mag.X = cat(2, magX{:}) * 1e6;   % T → µT
    raw.mag.Y = cat(2, magY{:}) * 1e6;
    raw.mag.Z = cat(2, magZ{:}) * 1e6;
    
    timestamps = cat(2, ts{:});
    clear accX accY accZ gyrX gyrY gyrZ magX magY magZ ts mot
    
    % Remove samples where all nine axes read zero (transmission gaps/failures)
    idx0 = ~raw.acc.X & ~raw.acc.Y & ~raw.acc.Z & ...
           ~raw.gyr.X & ~raw.gyr.Y & ~raw.gyr.Z & ...
           ~raw.mag.X & ~raw.mag.Y & ~raw.mag.Z;
    if any(idx0)
        flds = {'acc','gyr','mag'};
        axes = {'X','Y','Z'};
        for f = 1:3
            for a = 1:3
                raw.(flds{f}).(axes{a})(idx0) = [];
            end
        end
        timestamps(idx0) = [];
        fprintf('  Removed %d zero-padded samples.\n', sum(idx0));
    end
    
    % Correct for battery-swap / recording-gap jumps in the timestamp counter
    idxt = find(diff(timestamps) > 3);
    if ~isempty(idxt)
        timestamps(idxt+1:end) = timestamps(idxt+1:end) - ...
            (timestamps(idxt+1) - timestamps(idxt) - 1);
        fprintf('  Corrected %d timestamp gap(s) in Deuteron counter.\n', numel(idxt));
    end
    
    fs_mot   = MotionSensorConstants.AccelerometerFrequency;   % 1000 Hz
    raw.tsec = timestamps / fs_mot;
    raw.tsec = raw.tsec - raw.tsec(1);   % relativise to recording start
    
    raw.fs      = fs_mot;
    raw.source  = 'deuteron-mpu9250';
    raw.dof     = 9;
    raw.acclMax = opt.acclMax;
    raw.gyroMax = opt.gyroMax;
    raw.magMax  = opt.magMax * 1e6;   % µT
    
    % Document native chip frame so ProcessMotionSensors can apply the right rotation
    raw.chipFrame.acc_gyro = ...
        'MPU-9250: acc/gyro share one frame — X-forward, Y-right, Z-up (ENU-like).';
    raw.chipFrame.mag = ...
        'MPU-9250: magnetometer has a DIFFERENT frame — Z-down (NED-aligned). See datasheet §Orientation of Axes.';
    raw.chipFrame.note = ...
        'Acc/gyro must be remapped (permute+sign) to align with magnetometer before AHRS fusion. Do this in ProcessMotionSensors.';
    
    savePath = fullfile(opt.FolderProcDataMat, 'MotionData_raw.mat');
    save(savePath, 'raw', '-v7.3');
    fprintf('- Raw Deuteron motion data saved: %s\n', savePath)
end

function getfrom_INTAN(opt)
    % Read 3-DoF ADXL335 accelerometer from Intan AUX channels.
    % Converts counts → V → m/s² with bias calibration and decimates to 1 kHz.
    % Native sensor frame is preserved — no axis remapping.
    
    sessionDir = opt.PathRaw;
    assert(isfolder(sessionDir), 'Invalid session directory: %s', sessionDir);
    
    targetCh = dir(fullfile(sessionDir, '*AUX*.dat'));
    assert(~isempty(targetCh), 'No AUX*.dat files found in: %s', sessionDir);
    
    % Hardware constants
    par.fsRaw         = 10000;    % effective AUX sample rate after sub-sampling (Hz)
                                  % (Intan AUX runs at 1/3 amplifier rate; ~30kHz→10kHz)
    par.targetFs      = 1000;    % decimation target (Hz)
    par.voltsPerCount = 3.74e-5; % Intan RHD2000 ADC scaling (V/count)
    par.sensitivity   = 0.340;   % ADXL335 typical sensitivity (V/g at 3.3 V supply)
    par.g0            = 9.81;    % m/s² per g
    par.calibDur_s    = 1;       % duration of initial segment used for bias estimate (s)
    
    nAxes  = min(3, numel(targetCh));
    nCal   = round(par.calibDur_s * par.fsRaw);
    rawMat = [];   % [nSamples × nAxes] — filled per axis below
    
    disp('- Extracting INTAN AUX accelerometer data...')
    for k = 1:nAxes
        fid   = fopen(fullfile(sessionDir, targetCh(k).name), 'r');
            counts = fread(fid, inf, 'int16');
        fclose(fid);
        counts = counts(1:3:end);   % take every 3rd sample → 10 kHz
        rawMat(1:numel(counts), k) = counts;
    end
    
    % Counts → volts → m/s² with per-axis bias removal
    biasV   = median(rawMat(1:min(nCal,size(rawMat,1)), :), 1, 'omitnan');
    accel   = par.g0 * ((rawMat * par.voltsPerCount - biasV) / par.sensitivity);
    nSampl  = size(accel, 1);
    
    % Decimate from 10 kHz to 1 kHz using MATLAB's decimate (LP filter + downsample)
    dsFactor = round(par.fsRaw / par.targetFs);
    fprintf('  Decimating by factor %d (%d Hz → %d Hz)...\n', ...
        dsFactor, par.fsRaw, par.targetFs);
    accel_ds = zeros(floor(nSampl / dsFactor), nAxes);
    for k = 1:nAxes
        accel_ds(:, k) = decimate(accel(1:nSampl, k), dsFactor);
    end
    accel_ds = accel_ds(1:floor(nSampl / dsFactor), :);   % trim to exact length
    
    raw.acc.X = accel_ds(:, 1)';
    raw.acc.Y = accel_ds(:, 2)';
    raw.acc.Z = accel_ds(:, 3)';
    raw.gyr   = [];   % not available on ADXL335
    raw.mag   = [];   % not available on ADXL335
    
    nDs      = size(accel_ds, 1);
    raw.tsec = (0:nDs-1) / par.targetFs;
    raw.fs   = par.targetFs;
    
    raw.source = 'intan-adxl335';
    raw.dof    = 3;
    
    % Document native chip frame (ADXL335 datasheet)
    raw.chipFrame.acc  = ...
        'ADXL335: Z-up when flat. Axis-to-head mapping depends on headstage mounting orientation.';
    raw.chipFrame.note = ...
        'Apply headstage-specific mounting rotation in ProcessMotionSensors before pitch/roll estimation.';
    
    raw.meta.sessionDir = sessionDir;
    raw.meta.targetCh   = {targetCh.name};
    raw.meta.parameters = par;
    
    savePath = fullfile(opt.FolderProcDataMat, 'MotionData_raw.mat');
    save(savePath, 'raw', '-v7.3');
    fprintf('- Raw INTAN accelerometer data saved: %s\n', savePath)
end
