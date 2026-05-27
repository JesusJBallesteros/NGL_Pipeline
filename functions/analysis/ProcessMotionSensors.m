function ProcessMotionSensors(opt, input, varargin)
% ProcessMotionSensors  Analyse raw IMU data and estimate head orientation.
%
% =========================================================================
% *** DRAFT — pseudocode skeleton for NGL02. Not yet executable. ***
% =========================================================================
%
% PURPOSE:
%   Loads MotionData_raw.mat (written by GetMotionSensors in NGL01), applies
%   a configurable body-frame rotation to align the native chip axes with the
%   animal's head coordinate system, then branches on sensor DOF:
%
%   9-DoF (Deuteron MPU-9250):
%     magnetometer hard/soft iron correction (magcal) →
%     AHRS sensor fusion (ahrsfilter or Madgwick) →
%     heading quaternion + Euler angles →
%     plots + MotionData.mat
%
%   3-DoF (INTAN ADXL335):
%     low-pass gravity extraction → pitch + roll from atan2 →
%     optional gaze projection onto screen plane →
%     plots + MotionData.mat
%
%   Event-locked analysis (both platforms, optional):
%     estimate_pecking (pecking / head-jerk detection)
%
% USAGE:
%   ProcessMotionSensors(opt, input)
%   ProcessMotionSensors(opt, input, 'doGaze', true)   % INTAN + screen setup
%   Called from NGL02_postPhy when opt.GetMotionSensors = true.
%
% INPUTS:
%   opt    - options struct; relevant fields:
%              .FolderProcDataMat  folder containing MotionData_raw.mat
%              .motionOpt          (optional sub-struct, see §OPTIONS below)
%   input  - struct from set_default
%   varargin - name-value pairs passed to the gaze-projection block
%
% OUTPUT:
%   MotionData.mat in opt.FolderProcDataMat, containing:
%     processed.acc/gyr/mag   axis-aligned, filtered sensor data
%     processed.heading        (9-DoF only) heading quaternion [N×4]
%     processed.euler          (9-DoF only) [roll, pitch, yaw] in degrees [N×3]
%     processed.pitch/roll     (3-DoF only) from gravity decomposition
%     processed.tsec           time vector
%     processed.R_body         3×3 body-frame rotation matrix applied
%     processed.source         forwarded from raw struct
%
% =========================================================================
%
% §OPTIONS  (set in opt.motionOpt or as defaults here)
%
%   Axis alignment:
%     motionOpt.R_acc_to_body   3×3 rotation, aligns acc/gyro native frame
%                               to body frame (X-forward, Y-right, Z-down).
%                               Default: eye(3) (no rotation — MUST be set).
%     motionOpt.R_mag_to_body   3×3 rotation for magnetometer frame.
%                               Default: eye(3).
%
%   Magnetometer (9-DoF only):
%     motionOpt.MField_ref      Reference field magnitude (µT) for the
%                               recording site. Bochum: 19.7 µT.
%                               Used to validate magcal output.
%     motionOpt.doMagcal        true (default) — run magcal hard/soft iron fit.
%
%   AHRS (9-DoF only):
%     motionOpt.ahrsType        'ahrsfilter' (default) | 'madgwick'
%     motionOpt.GyroNoise       Gyroscope noise variance (rad/s). Default: 0.01.
%     motionOpt.AccelNoise      Accelerometer noise variance (m/s²). Default: 0.01.
%
%   Gravity / pitch-roll (3-DoF only):
%     motionOpt.gravityLP_Hz    Low-pass cutoff for gravity estimate (Hz). Default: 0.3.
%
%   Gaze projection (3-DoF + screen only):
%     motionOpt.doGaze          false (default) | true
%     motionOpt.screenCenter    [x,y,z] screen centre in world coords (m).
%     motionOpt.screenNormal    [nx,ny,nz] unit normal pointing toward subject.
%     motionOpt.screenWidth_m   Screen width (m).
%     motionOpt.screenHeight_m  Screen height (m).
%     motionOpt.screenResPx     [width_px, height_px].
%
%   Event-locked analysis:
%     motionOpt.doEthology      false (default) | true — run estimate_pecking.
%
% =========================================================================
%
% CALLS (when complete):
%   magcal, ahrsfilter (or imufilter + Madgwick), Deuteron_PlotMotionSensors,
%   estimate_pecking
%
% SEE ALSO:
%   GetMotionSensors (NGL01 raw extraction)
%
% Last modified 15.05.2026 (Jesus) — DRAFT

% =========================================================================
%% 0. Guard: skip if already processed
% =========================================================================
% PSEUDOCODE:
%   outPath = fullfile(opt.FolderProcDataMat, 'MotionData.mat');
%   if isfile(outPath) && <file is non-empty>
%       disp('MotionData.mat already exists. Loading and returning.')
%       load(outPath); return
%   end

% =========================================================================
%% 1. Load raw data
% =========================================================================
% PSEUDOCODE:
%   rawPath = fullfile(opt.FolderProcDataMat, 'MotionData_raw.mat');
%   assert(isfile(rawPath), 'MotionData_raw.mat not found. Run NGL01 first.');
%   raw = load(rawPath).raw;
%   fprintf('Loaded %s | source: %s | DoF: %d | fs: %d Hz\n', ...
%       rawPath, raw.source, raw.dof, raw.fs);

% =========================================================================
%% 2. Parse options / set defaults
% =========================================================================
% PSEUDOCODE:
%   if isfield(opt, 'motionOpt'), mopt = opt.motionOpt; else, mopt = struct(); end
%
%   % Axis alignment rotations — MUST be set by user for correct results
%   R_acc = getdefault(mopt, 'R_acc_to_body', eye(3));   % acc/gyro → body
%   R_mag = getdefault(mopt, 'R_mag_to_body', eye(3));   % mag → body
%
%   if isequal(R_acc, eye(3))
%       warning('NGL:identityRotation', ...
%           'R_acc_to_body is identity. Axis alignment may be incorrect. ' ...
%           'Set opt.motionOpt.R_acc_to_body to match your headstage mounting.');
%   end
%
%   % 9-DoF options
%   MField_ref  = getdefault(mopt, 'MField_ref',  19.7);       % µT Bochum
%   doMagcal    = getdefault(mopt, 'doMagcal',    true);
%   ahrsType    = getdefault(mopt, 'ahrsType',    'ahrsfilter');
%   GyroNoise   = getdefault(mopt, 'GyroNoise',   0.01);
%   AccelNoise  = getdefault(mopt, 'AccelNoise',  0.01);
%
%   % 3-DoF options
%   gravLP      = getdefault(mopt, 'gravityLP_Hz', 0.3);
%   doGaze      = getdefault(mopt, 'doGaze',       false);
%
%   % Ethology
%   doEthology  = getdefault(mopt, 'doEthology',   false);

% =========================================================================
%% 3. Apply body-frame axis rotation
% =========================================================================
%
% MPU-9250 SPECIFIC NOTE (9-DoF):
%   The acc/gyro chip and magnetometer chip have DIFFERENT axis orientations
%   on the MPU-9250 package (documented in datasheet "Orientation of Axes").
%   Two separate rotations are required — one for acc/gyro, one for mag.
%
%   Default remapping confirmed to work with current Deuteron board mounting:
%     acc_body = R_acc * [acc.X; acc.Y; acc.Z]
%     where R_acc maps native [X,Y,Z] → body [X_fwd, Y_right, Z_down]:
%     From Deuteron_estimateheading: gyr_body = [gyr.X, gyr.Z, -gyr.Y]
%                                    acc_body = [-acc.X, -acc.Z, acc.Y]
%     Expressed as a rotation matrix:
%       R_acc_default = [-1  0  0;   % body-X  = -native-X
%                         0  0  1;   % body-Y  =  native-Z
%                         0 -1  0];  % body-Z  = -native-Y  ← points down
%
%   Magnetometer is already NED-aligned (Z-down), so for single-shank, single-
%   area standard mounting R_mag = eye(3) may be sufficient.
%   Verify by checking magcal output matches MField_ref.
%
% ADXL335 SPECIFIC NOTE (3-DoF):
%   Native frame is Z-up. To convert to body Z-down: negate Z.
%   X and Y assignment depends on physical headstage mounting orientation
%   (which channel is forward, which is lateral). Must be set per-experiment.
%
% PSEUDOCODE:
%   acc_body = (R_acc * [raw.acc.X; raw.acc.Y; raw.acc.Z])';  % [N×3]
%   if raw.dof == 9
%       gyr_body = (R_acc * [raw.gyr.X; raw.gyr.Y; raw.gyr.Z])';
%       mag_body = (R_mag * [raw.mag.X; raw.mag.Y; raw.mag.Z])';
%   end

% =========================================================================
%% 4a. 9-DoF branch — magnetometer calibration + AHRS
% =========================================================================
% if raw.dof == 9

    %% 4a-i. Hard/soft iron correction (magcal)
    % PSEUDOCODE:
    %   if doMagcal
    %       [A, b, Mfield] = magcal(mag_body);
    %       fprintf('magcal: measured field = %.1f µT | reference = %.1f µT\n', ...
    %           Mfield*1e6, MField_ref);
    %       if abs(Mfield*1e6 - MField_ref) > 5
    %           warning('magcal result deviates >5 µT from reference. Check mounting.');
    %       end
    %       mag_cal = (mag_body - b) * A;   % [N×3] corrected
    %   else
    %       mag_cal = mag_body;
    %   end

    %% 4a-ii. AHRS sensor fusion → orientation quaternion
    % PSEUDOCODE:
    %   if strcmp(ahrsType, 'ahrsfilter')
    %       % Requires Sensor Fusion and Tracking Toolbox
    %       fuse = ahrsfilter('SampleRate', raw.fs, ...
    %                         'GyroscopeNoise', GyroNoise, ...
    %                         'AccelerometerNoise', AccelNoise);
    %       % Input: acc [m/s²], gyro [rad/s], mag [µT]
    %       quat = fuse(acc_body, deg2rad(gyr_body), mag_cal);   % [N×1 quaternion]
    %
    %   elseif strcmp(ahrsType, 'madgwick')
    %       % Pure-MATLAB Madgwick filter (no toolbox needed) — implement or bundle
    %       quat = madgwick_ahrs(acc_body, deg2rad(gyr_body), mag_cal, raw.fs, beta);
    %   end

    %% 4a-iii. Quaternion → Euler angles
    % PSEUDOCODE:
    %   euler_deg = eulerd(quat, 'ZYX', 'frame');   % [N×3] yaw, pitch, roll
    %   heading_deg = wrapTo360(euler_deg(:,1));     % compass heading

    %% 4a-iv. Plots
    % PSEUDOCODE:
    %   Deuteron_PlotMotionSensors(data_struct, raw.tsec, opt, 1);     % raw timeseries
    %   Deuteron_PlotMotionSensors(data_struct, raw.tsec, opt, 1, 1);  % magcorr overlay
    %   exportgraphics(gcf, fullfile(opt.FolderProcDataMat, 'motion_raw.png'), 'Resolution', 300);
    %   Deuteron_PlotMotionSensors(quat, raw.tsec, opt, 2);            % Euler angles
    %   exportgraphics(gcf, fullfile(opt.FolderProcDataMat, 'motion_ahrs.png'), 'Resolution', 300);

% end  % 9-DoF branch

% =========================================================================
%% 4b. 3-DoF branch — gravity decomposition → pitch + roll
% =========================================================================
% if raw.dof == 3

    %% 4b-i. Separate gravity (static) from dynamic acceleration
    % PSEUDOCODE:
    %   fs = raw.fs;
    %   [bLP, aLP] = butter(2, gravLP / (fs/2), 'low');
    %   grav = filtfilt(bLP, aLP, acc_body);   % [N×3] gravity estimate
    %   dyn  = acc_body - grav;                % [N×3] dynamic (head movements)
    %
    %   g0   = 9.81;
    %   gx   = grav(:,1) / g0;
    %   gy   = grav(:,2) / g0;
    %   gz   = grav(:,3) / g0;
    %
    %   roll_rad  = atan2(gy, gz);
    %   pitch_rad = atan2(-gx, sqrt(gy.^2 + gz.^2));
    %   roll_deg  = rad2deg(roll_rad);
    %   pitch_deg = rad2deg(pitch_rad);

    %% 4b-ii. Optional gaze projection onto screen plane
    % (Only meaningful for screen-based paradigms)
    % PSEUDOCODE:
    %   if doGaze
    %       <see existing getfrom_INTAN gaze block — move here verbatim>
    %       % Outputs: screenXY_m [N×2], screenXY_px [N×2], heatmap figures
    %       % Save to gaze_estimation.mat
    %   end

    %% 4b-iii. Plots
    % PSEUDOCODE:
    %   figure; plot(raw.tsec, acc_body);
    %   xlabel('Time (s)'); ylabel('Accel (m/s²)'); legend('X','Y','Z');
    %   title('Bias-corrected accel — body frame');
    %   exportgraphics(gcf, fullfile(opt.FolderProcDataMat, 'motion_accel.png'), 'Resolution', 300);
    %
    %   figure; plot(raw.tsec, pitch_deg, 'b', raw.tsec, roll_deg, 'r');
    %   xlabel('Time (s)'); ylabel('deg'); legend('Pitch','Roll');
    %   exportgraphics(gcf, fullfile(opt.FolderProcDataMat, 'motion_pitchroll.png'), 'Resolution', 300);

% end  % 3-DoF branch

% =========================================================================
%% 5. Event-locked ethology analysis (optional, both platforms)
% =========================================================================
% PSEUDOCODE:
%   if doEthology
%       load(fullfile(opt.FolderProcDataMat, 'EventRecord.mat'), 'EventRecord');
%       % estimate_pecking requires the signal and the EventRecord struct.
%       % peak_opt should be set via opt.motionOpt.peak_opt or a sub-struct.
%       if raw.dof == 9
%           estimate_pecking(data_struct, EventRecord, peak_opt, []);
%       else
%           estimate_pecking(acc_body, EventRecord, peak_opt);
%       end
%   end

% =========================================================================
%% 6. Assemble and save MotionData.mat
% =========================================================================
% PSEUDOCODE:
%   processed.source  = raw.source;
%   processed.dof     = raw.dof;
%   processed.tsec    = raw.tsec;
%   processed.fs      = raw.fs;
%   processed.R_body  = R_acc;
%   processed.acc     = acc_body;
%
%   if raw.dof == 9
%       processed.gyr     = gyr_body;
%       processed.mag     = mag_body;
%       processed.mag_cal = mag_cal;
%       processed.quat    = quat;
%       processed.euler   = euler_deg;
%       processed.heading = heading_deg;
%   else
%       processed.pitch   = pitch_deg;
%       processed.roll    = roll_deg;
%       processed.dyn     = dyn;
%       processed.grav    = grav;
%       if doGaze, processed.gaze = gaze_struct; end
%   end
%
%   save(fullfile(opt.FolderProcDataMat, 'MotionData.mat'), 'processed', '-v7.3');
%   disp('MotionData.mat saved.')

end


% =========================================================================
% Local helpers (to be implemented)
% =========================================================================

% function val = getdefault(s, field, default)
% % Return s.(field) if it exists, otherwise return default.
%     if isfield(s, field) && ~isempty(s.(field))
%         val = s.(field);
%     else
%         val = default;
%     end
% end

% function quat = madgwick_ahrs(acc, gyro, mag, fs, beta)
% % Pure-MATLAB Madgwick AHRS filter.
% % acc  [N×3] m/s², gyro [N×3] rad/s, mag [N×3] µT, fs Hz, beta tuning.
% % Returns quat [N×4] as [w, x, y, z].
%     <implementation or call to bundled Madgwick .m>
% end
