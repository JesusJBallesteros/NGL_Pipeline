function [magnitude_vector, valid_peaks, valid_peaks_wf] = estimate_pecking(MotionData, EventRecord, opt)
% The function takes the precessed data from the motion sensors, and keeps
% only that coming from accelerometers (3/9 from Deuteron, 3/3 fron INTAN).
% This data is bandpassed (1-50Hz) to remove drift and noise, then the
% combined magnitude of the three vectors is calculated, to simplify the
% orientation dependence. The detection is done via 'findpeaks' to locate
% candidate events, and then compared to the expected timing given in the
% EventRecord variable.
% INPUTS:
%   MotionData: struct with fields:
%       .acc.X,  data.acc.Y, data.acc.Z   (Nx1) accelerometer (m/s^2 or g; consistent units)
%       .gyr.X,  data.gyr.Y, data.gyr.Z   (Nx1) gyroscope (rad/s or deg/s; consistent units)
%       .t                                (optional Nx1 time vector)
%   EventRecord: as it is generated during preprocessing
%   opts: struct with fields
%       .AlignMode  = 'mpu9250_nedlike' or 'intan'
%       .SDs        = 16;  % Detection threshold, to adjust
%       .limitEvents = {'bhv'}; % event of choice
%       .eventdef.(peak_opt.limitEvents{1}) = 3; % and its decimal value
%       .ev_minus   = 3; % time to get before event
%       .ev_plus    = 2; % time to get after event
%       .fs         = 1000;  % Sampling frequency, Hz
%       .hpass      = 20; % standard high-cut, Hz
%       .filedir    = fullfile(opt.FolderProcDataMat);
%       .stpSize    = 50; % for PSH calculation
%       .binSize    = 200; % for PSH calculation 
%       .s_around   = 100/peak_opt.fs; % 100ms around peaks seem OK
%   magnitude_vector: magnitude vector calculated in another run or 
%                     from other sources other than the raw data.
% Jesus. 28.01.2026

if ~exist('opt','var') || ~isstruct(opt)
    error('opt must be provided as a struct with required fields');
end

% defaults
if ~isfield(opt,'fs'), opt.fs = 1000; end
if ~isfield(opt,'high_cut'), opt.high_cut = 50; end
if ~isfield(opt,'SDs'), opt.SDs = 16; end

if isstruct(MotionData)
    data = MotionData;
else
    data.acc.X = MotionData(:,1);
    data.acc.Y = MotionData(:,2);
    data.acc.Z = MotionData(:,3);
end

% Extract accelerometer data from MotionData
[a, info] = alignAxes(data, opt);
clear MotionData data

% High-pass accel magnitude to emphasize rapid impacts
[info.HPb,info.HPa] = butter(3, opt.hpass/(opt.fs/2), 'high');
aHP = filtfilt(info.HPb,info.HPa, a);

% jerk
j = zeros(size(aHP));
j(2:end-1,:) = (aHP(3:end,:) - aHP(1:end-2,:)) / (2*(1/opt.fs));

% Euclidean Magnitudes
jMag = sqrt(sum(j.^2,2));

% clean up
clear j aHP

% ensure column vectors
data = jMag(:);

% Add info
info.nSig = size(data,2);
info.halfW = 50;

%% Timelimits and event times
% By default look backwards from reward event
if isscalar(opt.limitEvents)
    end_idx = find(EventRecord.EventType == opt.eventdef.(opt.limitEvents{1}));
    tlim_idx(1:numel(end_idx),1) = EventRecord.TimeSecFromMidnight(end_idx)-opt.ev_minus; % Look backwards
    tlim_idx(1:numel(end_idx),2) = EventRecord.TimeSecFromMidnight(end_idx)+opt.ev_plus; % Look forward
    clear end_idx
end
% event times in ms
ev = EventRecord.TimeSecFromMidnight(EventRecord.EventType==opt.eventdef.(opt.limitEvents{1}));

% %% Calculate (positive) threshold based on magnitude mean and SD
% [locs, info] = locatePeaks(data, opt, info);
% 
% % Pick n last ones
% ptimes_peak = get_lastones(locs, ev, opt.peaks_bfEvent);

%% Outlier normalized energy detector. Perhaps to substitute the above?
[locs, info] =  detect_collisions(data, opt, info);

% Pick n last ones
ptimes_col = get_lastones(locs, ev, opt.peaks_bfEvent);

clear locs

%% Using the detected peaks/collisions, extract times ands wfs
[accelData_trial, jMag_trial, col_peaks, col_peaks_wf, accAxes_wf_col, accAxes_avg_col, info] = ...
    get_peaks_wfs(a, data, tlim_idx, opt, ptimes_col, info);

%% pseudo-KS4 detector
r = detect_jerkEvents(data, opt.fs, ...
        'ThresholdSD',  12, ...     % liberal first pass
        'ThresholdSD2', 16, ...     % tighter second pass
        'NumClusters',  'auto', ... % or 'auto'
        'CustomTemplates', ptimes_col, ...
        'PlotTrace',    true, ...
        'PlotWaveforms',true, ...
        'PlotPCs',      true);

%% Plot for Simple peak location
% With last valid peaks
f1 = plot_peakingEstimation(accelData_trial, jMag_trial, accAxes_wf, valid_peaks_wf, valid_peaks, accAxes_avg, opt, info);
    opt.filename = fullfile(opt.filedir, sprintf('Fig_findpeaks_%dSD.pdf', opt.SDs));
    exportgraphics(f1, opt.filename, 'ContentType','vector');

% With collision energy detector
f2 = plot_peakingEstimation(accelData_trial, jMag_trial, accAxes_wf_col, col_peaks_wf, col_peaks, accAxes_avg_col, opt, info);
    opt.filename = fullfile(opt.filedir, sprintf('Fig_collisions_%d.pdf', opt.thr));
    exportgraphics(f2, opt.filename, 'ContentType','vector');

close all

%% Detect events from templates, using either 'valid_peaks_wf' or 'col_peaks_wf' as templates.
if opt.usetemplatedetection
    % Options for 'detect_events'
    opt.kMAD = 10;
    opt.minAmplFraction = 0.5;
    opt.maxAmplFraction = 1.5;
    wf_to_use = col_peaks_wf; % or 'valid_peaks_wf'
    
    % Template matching and PCA analysis
    events_detected = detect_events(wf_to_use, data, opt, 1);
    ptimes_col = [events_detected.time]'/1000;

    % Then, re-run function 'get_peaks_wfs' to obtain the new detected times and wfs
    [accelData_trial, jMag_trial, ~, ~, col_peaks, col_peaks_wf, ~, accAxes_wf_col, ~, accAxes_avg_col, info] = ...
        get_peaks_wfs(a, data, tlim_idx, opt,ptimes_peak, ptimes_col, info);

    % Plot for events detected by template matching and PCA
    f3 = plot_peakingEstimation(accelData_trial, jMag_trial, accAxes_wf_col, col_peaks_wf, col_peaks, accAxes_avg_col, opt, info);
        opt.filename = fullfile(opt.filedir, sprintf('Fig_collisions_%d.pdf', opt.thr));
        exportgraphics(f3, opt.filename, 'ContentType','vector');
end

end

%% HELPER FUNCTIONS
function [a, info] = alignAxes(data, opt)
    axr = data.acc.X(:);
    ayr = data.acc.Y(:);
    azr = data.acc.Z(:);

    info = struct();
    info.mode = opt.AlignMode;
    
    if strcmpi(opt.AlignMode,'mpu9250_nedlike')
        % Common MPU-9250 reconciliation to a mag/NED-like convention:
        % [x;y;z]_world = [ y; x; -z ]_raw  (swap x/y, flip z)
        A2W = [0 1 0;
               1 0 0;
               0 0 -1];
        info.A2W = A2W;
    elseif strcmpi(opt.AlignMode,'intan')
        % X for/backwards, Y up/down, Z left/right
        A2W = [1 0 0;
               0 1 0;
               0 0 1]; 
        info.A2W = A2W;
    end
    
    a = (A2W * [axr ayr azr]')';
end

% function [locs, info] = locatePeaks(jerk, opt, info)
%     % Calculate (positive) threshold based on magnitude mean and SD
%     info.mu = mean(jerk);
%     info.sd = std(jerk);
%     info.sdThreshold = opt.SDs;
%     info.threshold = info.mu + opt.SDs*info.sd;
% 
%     % Locate candidate peaks of activity
%     [~,locs(:,1)] = findpeaks(jerk, ...
%                         opt.fs, ...
%                         'MinPeakHeight', info.threshold, ...
%                         'MinPeakDistance', opt.s_around);
% end

function [locs, info] = detect_collisions(jerk, opt, info)
    jerk = jerk(:); % Force Nx1
    evtime = 0.1*opt.fs;
    win = round(evtime/5); % 20 ms window
       
    % Calculate energy with sliding window
    E = conv(jerk.^2, ones(win,1), 'same');
    
    % Threshold on Z-Normalized energy
    E0 = median(E);
    Es = mad(E,1);
    Z = (E - E0) / max(Es,1e-9);

    % Peak detection on energy signal
    [~,locs(:,1)] = findpeaks(Z, ...
                        opt.fs, ...
                        'MinPeakHeight', opt.thr, ...
                        'MinPeakDistance', opt.s_around);

    info.medE = mean(jerk);
    info.madE = std(jerk);
    info.threshold = opt.thr;
end

function last_ptime = get_lastones(locs, ev, numlast)
    % For each 'ev' find the prev 'peaktimes'
    for i = 1:numel(ev)
        idx = find(locs < ev(i), numlast, "last"); % as many last ones
        if numlast>1
            idx(1:numlast-1) = []; % remove all but that n-last one
        end
        if isempty(idx)
            prevPeak(i) = NaN;
            peck_idx(i) = NaN;
        else
            % calculate time to event
            prevPeak(i) = ev(i) - locs(idx);
            if prevPeak(i) <= 0.25
                % Keep signal as possible peck
                peck_idx(i) = idx;
            else
                peck_idx(i) = NaN;
            end
        end
    end
    peck_idx(isnan(peck_idx)) = [];
    last_ptime = unique(locs(peck_idx));
end

function [accelData_trial, jMag_trial, col_peaks, col_peaks_wf, accAxes_wf_col, accAxes_avg_col, info] = ...
    get_peaks_wfs(a, data, tlim_idx, opt, ptimes_col, info)
    % Get waveforms and trial traces

    accelData_trial = cell(length(tlim_idx),1); % Accel data traces per trial
    jMag_trial      = cell(length(tlim_idx),info.nSig); % aMag & jMag traces per trial    
    % valid_peaks     = cell(length(tlim_idx),info.nSig); % valid peak times per trial
    col_peaks       = cell(length(tlim_idx),info.nSig); % valid peak times per trial
    % valid_peaks_wf  = NaN(1,info.halfW*2);
    col_peaks_wf    = NaN(1,info.halfW*2);        
    % accAxes_wf      = NaN(info.halfW*2,3,1);
    accAxes_wf_col  = NaN(info.halfW*2,3,1);
    % wf_count = 0;
    wf_count_col = 0;
    
    for t = 1:length(tlim_idx)
        if (tlim_idx(t,1)) < 1, continue, end
        
        % Whole trial traces
            seg_st  = round(tlim_idx(t,1)*opt.fs); % start sample
            seg_end = round(tlim_idx(t,2)*opt.fs)-1; % end sample
        
            % AccelData traces
            accelData_trial{t} = a(seg_st:seg_end,:);
        
            % Jerk magnitude trace
            jMag_trial{t} = data(seg_st:seg_end);
    
        % Find last-peak times falling within trial limits
            % % For first peak detector
            % if any(ptimes_peak >= tlim_idx(t,1) & ptimes_peak < tlim_idx(t,2)) % if any at all
            %     valid_peaks{t}(:,1) = ptimes_peak(ptimes_peak >= tlim_idx(t,1) & ptimes_peak < tlim_idx(t,2)); % absolute ts
            %     valid_peaks{t}(:,2) = valid_peaks{t}(:,1) - tlim_idx(t,2) + opt.ev_plus;  % relative to event time
            % else
            %     valid_peaks{t} = [NaN NaN];
            % end
            % valid_peaks = valid_peaks(:);
            % 
            % % Get all peak wfs
            % if ~isnan(valid_peaks{t}(1))
            %     wf_count = wf_count+1;
            %     wf_st  = round(valid_peaks{t}(1)*opt.fs) - info.halfW;
            %     wf_end = round(valid_peaks{t}(1)*opt.fs) + info.halfW-1;
            % 
            %     % retrieve waveforms from the acc magnitude
            %     valid_peaks_wf(wf_count,:) = data(wf_st:wf_end);
            %     % retrieve waveforms from the raw acc
            %     accAxes_wf(:,:,wf_count) = a(wf_st:wf_end,:);
            % end
            
            % For second peak detector            
            if any(ptimes_col >= tlim_idx(t,1) & ptimes_col < tlim_idx(t,2)) % if any at all
                col_peaks{t}(:,1) = ptimes_col(ptimes_col >= tlim_idx(t,1) & ptimes_col < tlim_idx(t,2)); % absolute ts
                col_peaks{t}(:,2) = col_peaks{t}(:,1) - tlim_idx(t,2) + opt.ev_plus;  % relative to event time
            else
                col_peaks{t} = [NaN NaN];
            end
            col_peaks = col_peaks(:);
        
            if ~isnan(col_peaks{t}(1))
                wf_count_col = wf_count_col+1;
                wf_st  = round(col_peaks{t}(1)*opt.fs) - info.halfW;
                wf_end = round(col_peaks{t}(1)*opt.fs) + info.halfW-1;
        
                % retrieve waveforms from the acc magnitude
                col_peaks_wf(wf_count_col,:) = data(wf_st:wf_end);
                % retrieve waveforms from the raw acc ax_fes
                accAxes_wf_col(:,:,wf_count_col) = a(wf_st:wf_end,:);
            end
    end
            
    % Averages of single axes data
    % accAxes_avg = mean(accAxes_wf,3); 
    accAxes_avg_col = mean(accAxes_wf_col,3); 
    
    % Recollect wf counts
    % info.wf_count = wf_count;
    info.wf_count_col = wf_count_col;
    % info.wf_count   = wf_count;
    info.wf_count_col = wf_count_col;
end