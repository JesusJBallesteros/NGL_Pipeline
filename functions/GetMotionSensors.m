function GetMotionSensors(opt, input)
% GetMotionSensors  Extract and process head-direction / accelerometer data.
%
% PURPOSE:
%   Dispatches motion-sensor extraction based on recording format:
%   - INTAN (fileperch): reads AUX*.dat accelerometer channels
%   - Deuteron (DT2/DF1): reads NEUR*.dat via Deuteron_extractData
%   Converts raw ADC values to physical units, applies AHRS processing to
%   estimate animal heading/position, plots the raw timeseries, and saves
%   a MotionData.mat file to opt.FolderProcDataMat.
%
% USAGE:
%   GetMotionSensors(opt, input)
%   Only called when opt.GetMotionSensors = true (gated in INTAN_PipelineWrapper
%   and Deuteron_PipelineWrapper).
%
% INPUTS:
%   opt    - options struct; relevant fields:
%              .FolderProcDataMat  output folder for MotionData.mat
%              .PathRaw            raw data folder (INTAN AUX*.dat location)
%   input  - struct; relevant field:
%              .sessions.info.fileformat  ('fileperch' | 'DT2' | 'DF1')
%
% OUTPUT:
%   MotionData.mat saved to opt.FolderProcDataMat; contains per-sensor
%   timeseries and AHRS-processed heading estimates.
%
% NOTES:
%   - Head-direction interpretation from AHRS requires careful calibration;
%     results should be verified manually before use in analysis.
%
% Jesus 08.05.2026
    
% options to parse into estimation functions
    % User decided
    peak_opt.limitEvents = {'bhv'}; % event of choice
        peak_opt.eventdef.(peak_opt.limitEvents{1}) = 3; % and its decimal value
        peak_opt.ev_minus = 2; % time to get before event
        peak_opt.ev_plus = 1; % time to get after event
    
    % Mostly settled
    peak_opt.fs = 1000;  % Sampling frequency, Hz
    peak_opt.hpass = 20; % standard high-cut, Hz
    peak_opt.filedir = fullfile(opt.FolderProcDataMat);
    peak_opt.stpSize = 50; % for PSH calculation
    peak_opt.binSize = 100; % for PSH calculation 
    peak_opt.s_around = 100/peak_opt.fs; % 100ms around peaks seem OK

if strcmp(input.sessions(input.run(1)).info.fileformat, 'fileperch')
    
    peak_opt.AlignMode = 'intan'; % for Intan
    getfrom_INTAN(opt, peak_opt)

elseif strcmp(input.sessions(input.run(1)).info.fileformat, 'DT2') || ...
       strcmp(input.sessions(input.run(1)).info.fileformat, 'DF1')

    peak_opt.AlignMode = 'mpu9250_nedlike'; % for Deuteron
    getfrom_Deuteron(opt, peak_opt)
end

end

% System dependent function to run
function getfrom_Deuteron(opt, peak_opt)
    %% Some local Parameters
    numFiles        = length(opt.myFiles);
    opt.stream      = 2;
      
    % Bochum Magnetic Field Horizontal Intensity. According to 
    % https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
    MField_Bochum = 19.7; % uTesla.
    
    % Sample rate for motion sensors is 1000Hz. 
    opt.fsmot          = 1000;         % Sample Rate of the feeded data (Hz)
    
    % % Gyro/Accel_Noise are determined from the hardware datasheets.
    opt.Gyro_Noise  = .01;   % Gyroscope Noise (variance value) in units of rad/s. (MPU-9250: 0.01 deg/sec)
    opt.Accel_Noise = .01;        % Accelerometer Noise(variance value) in units of m/s^2 (g). (MPU-9250: 8 mg)
    
    % The values for acclMax and gyroMax are chosen by the user. They can be found using the Event
    % File Viewer in the file started event. If not activelly changed, they should stay as follows:
    opt.acclMax = 2*MotionSensorConstants.G; % m/s^2, max value of selected range
    opt.gyroMax = 250;                       % degrees/s, max value of selected range
    opt.magMax  = MotionSensorConstants.Magnetometer9250Range; % Teslas, max value of selected range
        
    if ~isfile(fullfile(opt.FolderProcDataMat, strcat('MotionData.mat')))
        % Sort motion sensor data by data type.
        % Create structs
        data.acc    = struct('X', [], 'Y', [], 'Z', [], 'max', opt.acclMax);
        data.gyr    = struct('X', [], 'Y', [], 'Z', [], 'max', opt.gyroMax);
        data.mag    = struct('X', [], 'Y', [], 'Z', [], 'max', opt.magMax);
                    
        % Axes description. With board plugged on animal's head, and according to the sensor sheet:
          % Magnetometer: Y for vertical, X for AP and Z for DL.
          % Acc: +X for vertical up, +Y for AP forward and +Z for DL left.
          % Gyro: around +X for yaw left (look around),
          %       around +Y for roll left (rolling, 'croqueta'),
          %       around +Z for pitch down (nodding).
        % note: Feels like X-Y axes in Acc and Gyro are interchanged with the ones in the Magnetometer
        % Therefore:
        %     Head movement is detected mostly by gyroscope yaw (X)
        %     The magnetometer needs the two horizontal (planar) axes for head orientation: Z and Y
        %     Acc.X detects gravity acceleration (points down-up axes!)
        for i = 1:numFiles         
            if strlength(opt.myFiles(i).name) < 4 || ~startsWith(opt.myFiles(i).name,"NEUR")
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
                % NO MODIFICATION FROM RAW DATA JUST YET
                acc.X{i} = mot.Accelerometer.Data.X';
                acc.Y{i} = mot.Accelerometer.Data.Y';
                acc.Z{i} = mot.Accelerometer.Data.Z';
        
                gyr.X{i} = mot.Gyroscope.Data.X';
                gyr.Y{i} = mot.Gyroscope.Data.Y';
                gyr.Z{i} = mot.Gyroscope.Data.Z';
        
                mag.X{i} = mot.Magnetometer.Data.X';
                mag.Y{i} = mot.Magnetometer.Data.Y';
                mag.Z{i} = mot.Magnetometer.Data.Z';
 
                ts{i} = mot.Accelerometer.timestamps;
            end
        end

        data.acc.X = cat(2,acc.X{:});
        data.acc.Y = cat(2,acc.Y{:});
        data.acc.Z = cat(2,acc.Z{:});

        data.gyr.X = cat(2,gyr.X{:});
        data.gyr.Y = cat(2,gyr.Y{:});
        data.gyr.Z = cat(2,gyr.Z{:});

        data.mag.X = cat(2,mag.X{:});
        data.mag.Y = cat(2,mag.Y{:});
        data.mag.Z = cat(2,mag.Z{:});

        timestamps = cat(2,ts{:});
        clear ts acc gyr mag
        
        % Remove all timestamps where all readings are 0 (failsafe)
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
        
        % Find deadtimes, as changes of battery for long recordings will
        % lead to deadtimes in between the actual 
        idxt = find(diff(timestamps)>3);
        if ~isempty(idxt)
            timestamps(idxt+1:end) = timestamps(idxt+1:end)-(timestamps(idxt+1)-timestamps(idxt)-1);
        end
        
        % Calculate timestamps in seconds
        tsec = timestamps/opt.fsmot; % since midnight
        tsec = tsec - tsec(1); % relativize to recording
    
        % Plot sensors readings. RAW.
        % Run the plot function. 
        Deuteron_PlotMotionSensors(data, tsec, opt, 1);
        exportgraphics(gcf, fullfile(opt.FolderProcDataMat, strcat('motion_raw.png')), 'Resolution', 300)
        close gcf
    
        clear fid mot numFiles idx0 timestamps
    
        % Use the magcal function to obtain the correction coefficients for the
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
        Deuteron_PlotMotionSensors(data, tsec, opt, 1, 1);
        exportgraphics(gcf, fullfile(opt.FolderProcDataMat, strcat('motion_magcorr.png')), 'Resolution', 300)
        close gcf
    
        clear A b i
    
        % Save data to matfile
        save(fullfile(opt.FolderProcDataMat, strcat('MotionData.mat')), "data", "tsec", '-mat');
        disp('Magnetic-Corrected motion data saved.')
    else
        load(fullfile(opt.FolderProcDataMat, strcat('MotionData.mat')));
    end

    % get eventRecord
    load(fullfile(opt.FolderProcDataMat, 'EventRecord.mat'), 'EventRecord');

    % Peak estimation
    estimate_pecking(data, EventRecord, peak_opt, []);

    %% Estimate orientation and render a 3-D video.
    % makeOrientationVideoFromMotionData(data, tsec, opt, EventRecord)
    
   
end

function getfrom_INTAN(opt, peak_opt)
    % Finds aux-*-AUX*.dat files, converts to volts, creates accel matrix,
    % estimates pitch/roll, projects gaze to a screen plane, outputs 2D screen
    % coordinates (meters + pixels) and heatmaps.
    sessionDir      = opt.PathRaw;
    assert(isfolder(sessionDir), 'Invalid session directory.');

    % Identify which AUX channels correspond to accelerometer axes
    targetCh = dir('*AUX*.dat');

    %% Parameters
    par.fsAux_Hz        = {10000};        % AUX sampling rate (Hz)
    par.targetFs        = 1000;           % downsampled Hz, plenty for head motion
    par.voltsPerCount   = 3.74e-5;        % Intan documentation constant
    par.sensitivity_V_per_g = 0.340;      % ADXL335 typical ~0.300 V/g (set per-axis if available)
    par.g0 = 9.81;
    par.calibDur_s = 1;                   % initial calibration duration for bias (s)

    % Orientation & mounting: sensor axes -> head coordinate frame transform.
    % Default: assume accel columns are [X_forward, Y_right, Z_up] in sensor frame.
    par.sensorToHeadRot = eye(3);
    % % If your mounting is different, modify sensorToHeadRot (3x3 rotation).
    % par.sensorToHeadRot = [1 0 0; 0 1 0; 0 0 1]; % X for/backwards, Y up/down, Z left/right
          
    % Get data
    if ~isfile(fullfile(opt.FolderProcDataMat,'acceldata_raw.mat'))
        m = matfile(fullfile(opt.FolderProcDataMat, 'acceldata_raw.mat'),'Writable',true);
        for k = 1:3 % X, Y, Z Being +Y the G-axis when mounted in bird (Y = +1g, 9.81m/s2,  )
            auxN = str2double(targetCh(k).name(end-4));
            fprintf('Getting AUX channel %d ... \n', auxN);
            
            % Build filename (INTAN convention: aux-<prefix>-<native_channel_name>.dat)
            fname = fullfile(sessionDir, targetCh(k).name);
            if ~isfile(fname)
                error('File not found: %s', fname);
            end
            
            % Read raw (30KHz) int16 data, every third sample (10KHz)
            fid = fopen(fname, 'r');
                raw = fread(fid, inf, 'int16');
                nSampl = length(raw(1:3:end));
            fclose(fid);

            % Convert to volts (from documentation)
            m.accelData(1:nSampl,k) = par.voltsPerCount * double(raw(1:3:end));
        end
        clear raw

        tsec = (0:size(m.accelData,1)-1)/par.fsAux_Hz{1}; % in sec
    
        % RAW, uncorrected, 10 KHz figure
        f1 = tiledlayout;
            t1 = nexttile;
            plot(tsec(10000:110000), m.accelData(10000:110000,:)); % sample RAW data
            title(t1, 'RAW, uncorrected, 10 KHz');
            ylabel('V'), xlabel('session time (s)'), xlim(t1, [1 11]);

        % if numel(targetCh) > 3 % More than one headstage fix
        %     % Average axes 1,2,3 with 4,5,6 respectively
        %     accel(:,1) = mean(m.accelData(1:nSampl/3,[1 3]),2);
        %     accel(:,2) = mean(m.accelData(1:nSampl/3,[2 4]),2);
        %     accel(:,3) = mean(m.accelData(1:nSampl/3,[4 6]),2);
        %     targetCh(4:6) = [];
        % end

        % Volts -> accel (m/s^2) + bias calib
        disp('Processing Accelerometer data. It could take a moment...')
        nCal = max(1, round(par.calibDur_s * par.fsAux_Hz{1}));
        biasV = median(m.accelData(1:nCal,1:k), 1, 'omitnan');    % robust bias estimate
        accel = par.g0 * ((m.accelData - biasV) ./ par.sensitivity_V_per_g); % m/s^2
        
        % Rotate sensor frame to head frame if requested 
        % (DONE at estimate peaking)
        % m.accel_head = (par.sensorToHeadRot * accel.').';  % rows = samples
        m.accel_head = accel;
        clear accel
        
        % RAW, bias corrected, 10 KHz figure
            t2 = nexttile;
            plot(tsec(10000:110000), m.accel_head(50000:150000,:)); % sample UNBIAS data
            title(t2, 'RAW, bias corrected, 10 KHz');
            ylabel('m/s^2'), xlabel('session time (s)'), xlim(t2, [1 11]);
        
        % Save intermediate metadata
        meta = struct('sessionDir', sessionDir, 'parameters', par, ...
            'axesChannels', {targetCh});
    
        m.meta = meta;
    
    % Downsampling    
    dsFactor = round(par.fsAux_Hz{1} / par.targetFs);
    if dsFactor > 13
        dsFactor = factor(dsFactor);
    end

    % built-in filter+downsample
    if length(dsFactor) > 1
        % for f = 2:length(dsFactor)+1
        %     fprintf('Decimating by a factor of %d ... \n', dsFactor(f-1));
        %     for k = numel(targetCh):-1:1
        %         accel_head{f}(:,k) = decimate(m.accel_head{f-1}(1:nSampl,k), dsFactor(f-1));
        %         accel_head{f-1}(:,k) = [];
        %     end
        %     par.fsAux_Hz{f} = par.fsAux_Hz{f-1} / dsFactor(f-1);
        % end
    else
        fprintf('Decimating by a factor of %d ... \n', dsFactor);
        for k = 1:3
            accel_dwnsmpl{1}(:,k) = decimate(m.accel_head(1:nSampl,k), dsFactor); 
        end
        par.fsAux_Hz{1} = par.fsAux_Hz{1} / dsFactor;
    end

    % Keep last results only
    accel_dwnsmpl = accel_dwnsmpl{end};
    par.fsAux_Hz = par.fsAux_Hz{end};

    % Decimated time vector
    tsec = (0:size(accel_dwnsmpl,1)-1)/par.fsAux_Hz; % in sec

        % Downsampled, bias corrected, 1 KHz
        t3 = nexttile;
        plot(tsec(1000:11000),accel_dwnsmpl(1000:11000,:)); % sample DWNSMP data
        title(t3, 'Downsampled, bias corrected, 1 KHz');
        ylabel('m/s^2'), xlabel('session time (s)'), xlim([1 11]);

    m.accel_dwnsmpl = accel_dwnsmpl;
    m.tsec = tsec;
    meta = struct('sessionDir', sessionDir, 'parameters', par, ...
                  'axesChannels', {targetCh});
    m.meta = meta;
    peak_opt.filename =  fullfile(opt.FolderProcDataMat, 'Accel_data_sample.pdf');
    exportgraphics(f1, peak_opt.filename, 'ContentType', 'vector');

    else
        load(fullfile(opt.FolderProcDataMat,'acceldata_raw.mat'), 'accel_dwnsmpl');
    end

% get eventRecord
load(fullfile(opt.FolderProcDataMat, 'EventRecord.mat'), 'EventRecord');

% Peak Estimation
peak_opt.SDs = 10;  % Detection threshold, to adjust
peak_opt.thr = peak_opt.SDs*100;
peak_opt.peaks_bfEvent = 2; % N-peak-before-event to use as ground truth
peak_opt.usetemplatedetection = 1; % proceed with template matching and PCA
estimate_pecking(accel_dwnsmpl, EventRecord, peak_opt);

%% ON THE WORKS

    % Gaze projection geometry (meters). Edit to match your experiment:
    par.headPos         = [0.0, 0.00, 0.0]; % [x,y,z] head sensor origin (m) in world coords
    % Define screen plane using a point and a unit normal (screen-facing direction)
    par.screenCenter    = [0.00, 0.00, 0.00]; % screen center point in world coords (m)
    par.screenNormal    = [-1, 0, 0];        % screen normal (points toward subject). Unit-ish; will normalize
    par.screenWidth_m   = 0.50;              % physical screen width (m)
    par.screenHeight_m  = 0.30;              % physical screen height (m)
    par.screenResPx     = [1920, 1080];      % [width_px, height_px] - used to convert to pixels
    
    % Forward axis of sensor in sensor coords (direction that indicates "looking forward")
    par.fwdAxis_sensor  = [1; 0; 0];       % adjust if forward is a different axis
    
    % Yaw handling: 'fixed' or 'auto_global_sweep'
    par.yawMode = 'auto_global_sweep';    % 'fixed' or 'auto_global_sweep'
    par.yaw_fixed_deg = 0;                % used if yawMode == 'fixed'
    par.yawSweepRange_deg = 25;           % +/- range for auto sweep
    par.yawSweepStep_deg  = 2;            % step for sweep
    
    % Fixation detection for auto-sweep: dynamic accel threshold (m/s^2)
    par.fix_dyn_thresh_ms2 = 0.2;         % segments with dyn magnitude < threshold considered fixation
    par.fix_min_dur_s = 0.10;             % minimum fixation duration (s) to count
    
    % Filtering / gravity extraction
    par.gravityLP_Hz = 0.30;              % low-pass cutoff for gravity estimation (Hz)

    %% Gravity estimation (pitch & roll)
    if ~isfile(fullfile(opt.FolderProcDataMat,'grav_estimation.mat'))
        % low-pass filter to estimate gravity vector
        [bLP,aLP] = butter(2, par.gravityLP_Hz/(par.fsAux_Hz/2), 'low'); % apply along each column (use filtfilt for zero-phase)
        grav_ms2 = zeros(size(accel_dwnsmpl));
        for k = 1:numel(targetCh)
            grav_ms2(:,k) = filtfilt(bLP, aLP, accel_dwnsmpl(:,k));
        end
        dyn_ms2 = accel_dwnsmpl - grav_ms2;
        
        % compute pitch & roll from gravity (assumes head axes: x forward, y right, z up)
        gx = grav_ms2(:,1) / par.g0;
        gy = grav_ms2(:,2) / par.g0;
        gz = grav_ms2(:,3) / par.g0;
    
        roll_rad  = atan2(gy, gz);
        pitch_rad = atan2(-gx, sqrt(gy.^2 + gz.^2));
        roll_deg  = rad2deg(roll_rad);
        pitch_deg = rad2deg(pitch_rad);
        clear gx gy gz

        % Quick sanity plots
        figure('Name','Raw accel per axis','Color','w');
        plot(t, accel_dwnsmpl);
        xlabel('Time (s)'); ylabel('Accel (m/s^2)');
        legend('X','Y','Z'); title('Bias-corrected accel (downsampled)');

        figure('Name','Gravity vs raw accel','Color','w');
        subplot(3,1,1); plot(t, accel_dwnsmpl(:,1),'k',t,grav_ms2(:,1),'r'); ylabel('X');
        legend('raw','gravity');
        subplot(3,1,2); plot(t, accel_dwnsmpl(:,2),'k',t,grav_ms2(:,2),'r'); ylabel('Y');
        subplot(3,1,3); plot(t, accel_dwnsmpl(:,3),'k',t,grav_ms2(:,3),'r'); ylabel('Z');
        xlabel('Time (s)'); sgtitle('Accel vs gravity estimate');

        figure('Name','Pitch & Roll','Color','w');
        plot(t,rad2deg(pitch_rad),'b',t,rad2deg(roll_rad),'r');
        xlabel('Time (s)'); ylabel('deg');
        legend('Pitch','Roll');
        title('Head orientation');
    
        % Save
        disp('Saving gravity vectors and roll and pitch estimations ...');
        save(fullfile(opt.FolderProcDataMat,'grav_estimation.mat'), 't', 'grav_ms2','dyn_ms2','roll_rad','pitch_rad','roll_deg','pitch_deg','meta');
        clear grav_ms2 accel_head pitch_deg roll_deg

    else
        load(fullfile(opt.FolderProcDataMat,'grav_estimation.mat'), 'dyn_ms2','roll_rad','pitch_rad','meta');    
    end

    %% Fixation detection (for auto yaw sweep)
    if iscell(par.fsAux_Hz)
        par.fsAux_Hz = par.fsAux_Hz{end};
    end
    amag = sqrt(sum(dyn_ms2.^2,2));
    isFix = amag < par.fix_dyn_thresh_ms2;
    % enforce minimum duration
    minSamples = max(1, round(par.fix_min_dur_s * par.fsAux_Hz));
    isFix = movmax(isFix, minSamples) > 0;

    %% Yaw handling / selection
    % prepare screen geometry
    par.screenNormal = par.screenNormal(:) / norm(par.screenNormal(:));

    % build orthonormal axes on the screen plane for coordinate mapping:
    % choose screen_u (width direction) and screen_v (height direction)
    % naive: pick an arbitrary up vector (world Z) and orthonormalize
    worldUp = [0; 0; 1];
    screen_u = cross(worldUp, par.screenNormal(:));
    if norm(screen_u) < 1e-6
        % screen normal aligned with worldUp -> pick worldY as tangent
        screen_u = [0; 1; 0];
    else
        screen_u = screen_u / norm(screen_u);
    end
    screen_v = cross(par.screenNormal(:), screen_u);
    screen_v = screen_v / norm(screen_v);
    
    % screen corner extents in meters relative to screenCenter:
    halfW = par.screenWidth_m/2;
    halfH = par.screenHeight_m/2;

    % pick yaw according to selected mode
    switch lower(par.yawMode)
        case 'fixed'
            yaw_offset_deg = par.yaw_fixed_deg;
        case 'auto_global_sweep'
            yawVals = -par.yawSweepRange_deg:par.yawSweepStep_deg:par.yawSweepRange_deg;
            fracs = zeros(size(yawVals));
            for ii = 1:numel(yawVals)
                fracs(ii) = evaluate_yaw_fraction(yawVals(ii), ...
                    pitch_rad, roll_rad, par, halfW, halfH, isFix, screen_u, screen_v);
            end
            [~, idxMax] = max(fracs);
            yaw_offset_deg = yawVals(idxMax);
            fprintf('Auto yaw sweep selected yaw_offset = %.2f deg (max fraction on-screen = %.3f)\n', yaw_offset_deg, fracs(idxMax));
        otherwise
            error('Unknown yawMode: %s', yawMode);
    end
    clear isFix fracs idxMax

    %% Compute final gaze intersections
    yaw_rad = deg2rad(yaw_offset_deg);
    R_all = composeRotationMatrix(yaw_rad, pitch_rad, roll_rad); % 3x3xN
    N = size(R_all,3);
    dir_world = zeros(3,N);
    for ii = 1:N
        dir_world(:,ii) = squeeze(R_all(:,:,ii)) * par.fwdAxis_sensor(:);
    end
    clear R_all yaw_rad pitch_rad roll_rad
    
    % intersect ray P + t*d with screen plane dot((screenCenter-P), n) / dot(d,n)
    numer = dot((par.screenCenter(:)-par.headPos(:)), par.screenNormal(:));
    denom = sum(dir_world .* repmat(par.screenNormal(:),1,N), 1);
    t_inter = numer ./ denom;               % 1xN
    clear denom numer

    inter3 = par.headPos(:) + dir_world .* t_inter; % 3 x N
    inter3 = inter3.';                     % Nx3
    clear dir_world 
    
    % flag points in front of sensor and on-screen
    u_coords = (inter3 - par.screenCenter) * screen_u;  % Nx1
    v_coords = (inter3 - par.screenCenter) * screen_v;  % Nx1
    onScreen = (t_inter' > 0) & (u_coords >= -halfW) & (u_coords <= halfW) & (v_coords >= -halfH) & (v_coords <= halfH);
    clear inter3 t_inter

    % Convert to screen-local meters: origin at screen center, +u right, +v up
    screenXY_m = [u_coords, v_coords]; % Nx2 (meters)
    clear u_coords v_coords
    
    % Convert meters to pixels: map u ∈ [-halfW, halfW] → [1, screenResPx(1)]
    u_px = ((screenXY_m(:,1) + halfW) / (2*halfW)) * (par.screenResPx(1)-1) + 1;
    % map v ∈ [ -halfH, halfH ] to pixel rows. Pixel y: 1 = top; we'll put +v up -> smaller row index
    v_px = ((halfH - screenXY_m(:,2)) / (2*halfH)) * (par.screenResPx(2)-1) + 1;
    
    screenXY_px = [u_px, v_px];
    clear u_px v_px

    % Mask invalid points (off-screen) to NaN
    screenXY_m(~onScreen, :) = NaN;
    screenXY_px(~onScreen, :) = NaN;

    % Save
    save(fullfile(opt.FolderProcDataMat,'gaze_stimation.mat'), 'screenXY_px', 'screenXY_m', 'meta');
    clear onScreen

    %% Produce heatmaps
    % 1) Heatmap in meters (coarse grid)
    gridResX = 100;
    gridResY = 60; % bins in width/height for meters heatmap
    xedges = linspace(-halfW, halfW, gridResX);
    yedges = linspace(-halfH, halfH, gridResY);
    validIdx = ~any(isnan(screenXY_m),2);
    counts_m = histcounts2(screenXY_m(validIdx,1), screenXY_m(validIdx,2), xedges, yedges);
    counts_m = counts_m'; % for imagesc orientation

    figure('Name','Gaze Heatmap (meters)','Color','w','Position',[50 50 700 500]);
        imagesc(xedges, yedges, counts_m);
        axis xy;
        hold on;
        xlabel('u (m)');
        ylabel('v (m)');
        title('Gaze heatmap (meters)');
        colormap parula;
        colorbar;
        % draw screen rectangle
        rectangle('Position',[-halfW, -halfH, par.screenWidth_m, par.screenHeight_m], 'EdgeColor','w','LineWidth',1.5);

    % 2) Heatmap in pixels (full resolution but binned)
    figure('Name','Gaze Heatmap (pixels)','Color','w','Position',[800 50 700 500]);
        % Use coarser pixel bin to speed plotting
        binPx = 8;
        xedges_px = 1:binPx:par.screenResPx(1)+1;
        yedges_px = 1:binPx:par.screenResPx(2)+1;
        validIdx_px = ~any(isnan(screenXY_px),2);
        counts_px = histcounts2(screenXY_px(validIdx_px,1), screenXY_px(validIdx_px,2), xedges_px, yedges_px);
        counts_px = counts_px';
        imagesc(xedges_px, yedges_px, counts_px);
        axis xy;
        hold on;
        xlabel('pixel x');
        ylabel('pixel y');
        title('Gaze heatmap (pixels, binned)');
        colormap hot;
        colorbar;
        clear validIdx_px

    % Overlay sample gaze points (downsampled)
    % ds = max(1, round(N/1000));
    figure('Name','Sample Gaze Overlay','Color','w');
        plot(screenXY_m(validIdx,1), screenXY_m(validIdx,2), '.', 'MarkerSize', 4); axis equal;
        xlabel('u (m)');
        ylabel('v (m)');
        title('Sample gaze points (meters)');
        grid on;
        xlim([-halfW halfW]);
        ylim([-halfH halfH]);
    clear validIdx

end    

% Helpers
% given yaw_deg returns fraction of gaze points on-screen during fixations
function [frac, R] = evaluate_yaw_fraction(yaw_deg, pitch_rad, roll_rad, par, halfW, halfH, isFix, screen_u, screen_v)
    yaw_rad = deg2rad(yaw_deg);
    
    % build rotation matrices per sample using yaw + pitch + roll
    % Note: rotation composition uses R = R_yaw * R_pitch * R_roll with intrinsic rotations
    R = composeRotationMatrix(yaw_rad, pitch_rad, roll_rad); % 3x3xN
    dir_world = zeros(size(R,1), size(R,3));
    for ii = 1:size(R,3)
        d = squeeze(R(:,:,ii)) * par.fwdAxis_sensor(:);
        dir_world(:,ii) = d(:);
    end

    clear yaw_rad pitch_rad roll_rad yaw_deg

    % intersection with plane for fixation times
    P = par.headPos(:);
    numer = dot((par.screenCenter(:)-P), par.screenNormal(:));
    denom = sum(dir_world .* repmat(par.screenNormal(:),1,size(dir_world,2)), 1);
    t_inter = numer ./ denom;  % scalar / vector division will broadcast
    inter = P + dir_world .* t_inter; % 3xN
    clear dir_world numer denom
    
    % convert to screen u/v coordinates
    u = (inter' - par.screenCenter) * screen_u; % Nx1
    clear screen_u
    v = (inter' - par.screenCenter) * screen_v; % Nx1
    clear screen_v inter
    
    % check on-screen and in front (t_inter>0)
    onScreen = (t_inter' > 0) & (u >= -halfW) & (u <= halfW) & (v >= -halfH) & (v <= halfH);
    
    % count only fixation frames
    onScreen_fix = onScreen & isFix;
    frac = sum(onScreen_fix) / max(1, sum(isFix