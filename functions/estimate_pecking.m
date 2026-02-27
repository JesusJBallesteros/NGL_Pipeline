function [magnitude_vector, valid_peaks, valid_peaks_wf] = estimate_pecking(MotionData, tsec, EventRecord, opt, magnitude_vector)
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
%       .mag.x,  data.mag.y, data.mag.z   (optional; not required for detection)
%       .t                                (optional Nx1 time vector)
%   EventRecord: as it is generated during preprocessing
%   opts: (optional) struct with fields (all optional)
%       .AlignMode = 'mpu9250_nedlike' (default) or 'none'
%
%   magnitude_vector: magnitude vector calculated in another run or 
%                     from other sources other than the raw data.
% Jesus. 28.01.2026

if ~exist('opt','var') || ~isstruct(opt)
    error('opt must be provided as a struct with required fields');
end

% defaults
if ~isfield(opt,'fs'), opt.fs = 1000; end
if ~isfield(opt,'filter'), opt.filter = 'filtfilt'; end
if ~isfield(opt,'low_cut'), opt.low_cut = 1; end
if ~isfield(opt,'high_cut'), opt.high_cut = 50; end
if ~isfield(opt,'sgolayOrder'), opt.sgolayOrder = 3; end
if ~isfield(opt,'frameLen'), opt.frameLen = 11; end
if ~isfield(opt,'AlignMode'), opt.AlignMode = 'mpu9250_nedlike'; end

if isempty(magnitude_vector)
% Extract accelerometer data from MotionData
    if isstruct(MotionData)
        % Axis alignment (accel/gyro)
        [ax, ay, az, gx, gy, gz, info] = alignAxes(MotionData, opt);
    else
        warning('MotionData must be provided as a struct with required fields');
        data.acc.X = MotionData(:,1);
        data.acc.Y = MotionData(:,2);
        data.acc.Z = MotionData(:,3);
        data.gyr.X = NaN;
        data.gyr.Y = NaN;
        data.gyr.Z = NaN;

        [ax, ay, az, gx, gy, gz, info] = alignAxes(data, opt);
        clear data
    end
    clear MotionData

    % Core signals
    a = [ax ay az];
    g = [gx gy gz];
    
    aMag = sqrt(sum(a.^2,2)); % includes gravity
    gMag = sqrt(sum(g.^2,2));

    % Filter Acc and Gyr data
    if strcmp(opt.filter,'filtfilt')
        % High-pass accel magnitude to emphasize rapid impacts (removes g)
        [info.HPb,info.HPa] = butter(3, 8/(opt.fs/2), 'high');
        aHPmag = filtfilt(info.HPb,info.HPa, aMag);
        
        % Band-pass accel magnitude for impact energy (removes g)
        [info.BPb,info.BPa] = butter(3, [20 80]/(opt.fs/2), 'bandpass');
        aBPmag = filtfilt(info.BPb,info.BPa, aMag);

    end

    % Calculate jerk magnitude
    jerk = [0; diff(aBPmag)] * opt.fs;

    % Smoothed HF energy (moving RMS) over a short window
    info.win = max(3, round(0.05*opt.fs)); % ~30 ms
    hfRMS = sqrt(movmean(aBPmag.^2, info.win, 'Endpoints','shrink'));
end

% Fourier tests
% Inputs: hfRMS, jerk, aHPmag (vectors), Fs (sampling frequency in Hz)
data = [hfRMS(:), jerk(:), aHPmag(:)]; % ensure column vectors
info.signals = {'hfRMS','jerk','aHPmag'};
nSig = size(data,2);
    
clear aBPmag aHPmag hfRMS jerk aMag gMag a g

% % PWELCH parameters
% winLen = round(opt.fs * 1);    % 1 second window in samples
% win = hann(winLen);            % window vector (Hann)
% noverlap = round(0.5 * winLen);% 50% overlap (adjust if desired)
% nfft = max(256, 2^nextpow2(winLen)); % nfft at least 256 or pow2 of win
% freqRange = [1 100];           % Hz
% 
% % Compute PSD with pwelch for each signal over 1-100 Hz
% Pxx = cell(1,nSig);
% f = cell(1,nSig);
% for k = 1:nSig
%     x = data(:,k) - mean(data(:,k));
%     [Pxx{k}, f{k}] = pwelch(x, win, noverlap, nfft, opt.fs, "onesided");
% end
% 
% % Plot linear PSD
% figure;
% for k = 1:nSig
%     subplot(nSig,1,k)
%     plot(f{k}, Pxx{k}, 'LineWidth', 1.2)
%     xlabel('Frequency (Hz)')
%     ylabel('PSD (V^2/Hz)')
%     title(['pwelch PSD (linear): ' signals{k}])
%     grid on
% end
% 
% % Optional: plot in dB/Hz
% figure;
% for k = 1:nSig
%     subplot(nSig,1,k)
%     plot(f{k}, 10*log10(Pxx{k}), 'LineWidth', 1.2)
%     xlabel('Frequency (Hz)')
%     ylabel('PSD (dB/Hz)')
%     title(['pwelch PSD (dB): ' signals{k}])
%     grid on
% end

% Calculate (positive) threshold based on magnitude mean and SD
for k = 1:nSig
    info.mu{k} = mean(data(:,k));
    info.sd{k} = std(data(:,k));
    threshold{k} = info.mu{k} + opt.SDs*info.sd{k};

    %% Locate candidate peaks of activity
    [peakVal{k}, ptimes{k}(:,1)] = findpeaks(data(:,k), opt.fs, ...
                                          'MinPeakHeight', threshold{k},...
                                          'MinPeakDistance', opt.s_around);
                                            %, 'MinPeakProminence', threshold*1.5); 
    % Convert to indices
    peakIdx{k} = round(ptimes{k} * opt.fs);
    peakIdx{k} = max(1, min(numel(data(:,k)), peakIdx{k})); % clamp
    
    % PeakVal distribution check
    big_peak_idx{k} = peakVal{k} > median(peakVal{k});
    
    peaktimes{k,1} = ptimes{k}(big_peak_idx{k});
    peaktimes{k,2} = ptimes{k}(~big_peak_idx{k});

    %% find the previous peak to an ev
    % For each 'ev' find the prev 'peaktimes'
    ev = EventRecord.TimeSecFromMidnight(EventRecord.EventType==7);
    for pp = 1:2
        for i = 1:numel(ev)
            idx = find(peaktimes{k,pp} < ev(i),1,"last");
            if isempty(idx)
                prevPeak{k,pp}(:,i) = NaN;
                peak_idx{k,pp}(:,i) = NaN;
            else
                % Store the previous peak time for the current event
                prevPeak{k,pp}(:,i) = ev(i) - peaktimes{k,pp}(idx);
                peak_idx{k,pp}(:,i) = idx;
            end
        end   
    % distr_prevPeak{k}(:,1) = prevPeak(k,:);
    % p_ev_prec_by_Speak = sum(dist_prevPeak(:,1) <= 1) / numel(ev); 
    end
end

%% Timelimits
% By default look backwards from reward event
if isscalar(opt.limitEvents)
    end_idx = find(EventRecord.EventType == opt.eventdef.(opt.limitEvents{1}));
    tlim_idx(1:numel(end_idx),1) = EventRecord.TimeSecFromMidnight(end_idx)-opt.ev_minus; % Look backwards
    tlim_idx(1:numel(end_idx),2) = EventRecord.TimeSecFromMidnight(end_idx)+opt.ev_plus; % Look backwards
    clear end_idx
end

% Keep events happening only at expected response times
info.halfW = 50;
wf_count=0;
SorB = 1; % 2 for Big
for s = 1:nSig
    %opt.filename = fullfile(opt.filedir, sprintf('Figure_%s_%dSD_%ds_%d.pdf', opt.limitEvents{1}, opt.SDs, opt.ev_minus+opt.ev_plus, s));

    for t = 1:length(tlim_idx)
        if (tlim_idx(t,1))*opt.fs < 1, continue, end
        % find peaks in signal between time limits
        valid_peaks{t,1} = peaktimes{s,SorB}(peaktimes{s,SorB}(:,1) >= tlim_idx(t,1) & peaktimes{s,SorB}(:,1) <= tlim_idx(t,2));
        % relativize them to the event time
        valid_peaks{t,2} = valid_peaks{t,1} - (tlim_idx(t, 2)-opt.ev_plus);
        
        if s==1
            % obtain accelData traces for the whole time considered
            accelData_trial{t,1}(:,1) = ax((tlim_idx(t,1))*opt.fs:((tlim_idx(t,2))*opt.fs)-1); 
            accelData_trial{t,1}(:,2) = ay((tlim_idx(t,1))*opt.fs:((tlim_idx(t,2))*opt.fs)-1); 
            accelData_trial{t,1}(:,3) = az((tlim_idx(t,1))*opt.fs:((tlim_idx(t,2))*opt.fs)-1); 
        
            % gyrData_trial{t,1}(:,1) = gx((tlim_idx(t,1))*opt.fs:((tlim_idx(t,2))*opt.fs)-1); 
            % gyrData_trial{t,1}(:,2) = gy((tlim_idx(t,1))*opt.fs:((tlim_idx(t,2))*opt.fs)-1); 
            % gyrData_trial{t,1}(:,3) = gz((tlim_idx(t,1))*opt.fs:((tlim_idx(t,2))*opt.fs)-1); 
        end

        magnitude_vector_trial{t,:} = data((tlim_idx(t,1))*opt.fs:((tlim_idx(t,2))*opt.fs)-1, s);
        if ~isempty(valid_peaks{t,1})
            if any(ismember(valid_peaks{t,1},peaktimes{s,SorB}))
                % for tt = 1:length(valid_peaks{t,1})
                %     wf_count=wf_count+1;
                %     % retrieve waveforms from the acc magnitude
                %     valid_peaks_wf(wf_count,:) = ...
                %         magnitude_vector(Speaktimes(Speaktimes==valid_peaks{t,1}(tt),1)*opt.fs-49:Speaktimes(Speaktimes==valid_peaks{t,1}(tt),1)*opt.fs+50);
                %     % retrieve waveforms from the raw acc ax_fes
                %     accAxes_wf(:,1,wf_count) = ...
                %         ax_f(Speaktimes(Speaktimes==valid_peaks{t,1}(tt),1)*opt.fs-49:Speaktimes(Speaktimes==valid_peaks{t,1}(tt),1)*opt.fs+50);
                %     accAxes_wf(:,2,wf_count) = ...
                %         ay_f(Speaktimes(Speaktimes==valid_peaks{t,1}(tt),1)*opt.fs-49:Speaktimes(Speaktimes==valid_peaks{t,1}(tt),1)*opt.fs+50);
                %     accAxes_wf(:,3,wf_count) = ...
                %         az_f(Speaktimes(Speaktimes==valid_peaks{t,1}(tt),1)*opt.fs-49:Speaktimes(Speaktimes==valid_peaks{t,1}(tt),1)*opt.fs+50);
                % 
                % end
                for k = 1:numel(peakIdx{s})
                    st = peakIdx{s}(k) - info.halfW + 1;
                    e = peakIdx{s}(k) + info.halfW;
                    if st < 1 || e > size(data(:,s),1)
                        continue
                    end
                    wf_count = wf_count + 1;
                    if s==1
                        accAxes_wf(:,1,wf_count) = ax(st:e);
                        accAxes_wf(:,2,wf_count) = ay(st:e);
                        accAxes_wf(:,3,wf_count) = az(st:e);
                        % gyrAxes_wf(:,1,wf_count) = gx(st:e);
                        % gyrAxes_wf(:,2,wf_count) = gy(st:e);
                        % gyrAxes_wf(:,3,wf_count) = gz(st:e);
                    end

                    valid_peaks_wf(wf_count,:) = data(st:e, s);
                end
            end
        end
    end

    %% Averages
    % of single axes data
    accAxes_avg = mean(accAxes_wf,3);
    AvWf_plot = normalize(accAxes_avg, 'scale')*wf_count/15;
    
    % gyrData_trial = cellfun(@(x,y) x/100, gyrData_trial, 'UniformOutput',false);
    % gyrAxes_avg = mean(gyrAxes_wf,3);
    % AvWfG_plot = normalize(gyrAxes_wf, 'scale')*wf_count/15;
    
    % and the acc. magnitude waveforms
    valid_aligned_ms = cellfun(@(x) x*opt.fs, valid_peaks(:,2), 'UniformOutput', false);
    
    %% Plots
    red_rate = 500;
    figure
    fig = tiledlayout('flow',TileSpacing='tight', Padding='tight');
    
        % sample filtered data
        r = randperm(numel(data(:,s))-10000,1);
        t0 = nexttile(1, [1 6]);
            plot(tsec(r:r+9999), data(r:r+9999,s)/red_rate, LineWidth=2, Color=[.5 .5 .5]); hold on
            plot(tsec(r:r+9999), [ax(r:r+9999), ay(r:r+9999), az(r:r+9999)]'); 
            str1 = sprintf('thr. (%i SD)', opt.SDs);
            yline(threshold{s}/red_rate,'--', str1, FontSize=10);
            title('Downsampled, bias corrected, 1 KHz, filtered');
            ylabel('m/s^3'), xlabel('session time (s)');
            xlim(t0, [tsec(r) tsec(r+9999)])
            box off
    
    
        % Plot combined examples of 3-axis signal 
        r = randperm(size(accAxes_wf,3),6);
        for i = 1:length(r)
            nexttile(i+6);
            plot(valid_peaks_wf(r(i),:)/red_rate, LineWidth=2, Color=[.7 .7 .7]); hold on
            plot(accAxes_wf(:,:,r(i))); 
            xlabel('ms'), xlim([0 100]), xticks(25:25:75), xticklabels(-25:25:25);
            ylim([-20 20]);
            if i==1, ylabel('m/s^2'), yticks(-20:20:80);
            else,   ylabel([]), yticks([]); end
            if i==3, title('examples X/Y/Z traces', FontSize=14, FontWeight='bold'); end
            yline(threshold{s}/red_rate,'--');
            set(gca, FontSize=14, LineWidth=1.5, TickDir='in'); box off
        end
        clear i r
    
        % Plot individual axes as imagesc, overlay average waveform
        tx = nexttile(13, [1 2]);
            imagesc(squeeze(accAxes_wf(:,1,:))'); hold on
                colormap cool
            plot((wf_count/2)+AvWf_plot(:,1), LineWidth=2, Color='k');
            xlabel('ms'), xlim([20 80]), xticks(25:25:75), xticklabels(-25:25:75);
            ylabel('Trace #');
            box off
    
        ty = nexttile(14, [1 2]);
            imagesc(squeeze(accAxes_wf(:,2,:))'); hold on
                colormap cool
            plot((wf_count/2)+AvWf_plot(:,2), LineWidth=2, Color='k');
            xlabel('ms'), xlim([0 100]), xticks(25:25:75), xticklabels(-25:25:75);
            ylabel([]), yticks([]), yticklabels([]);
            title('All XYZ-axis traces at detection');        
            box off
    
        tz = nexttile(15, [1 2]);
            imagesc(squeeze(accAxes_wf(:,3,:))'); hold on
                colormap cool
            plot((wf_count/2)+AvWf_plot(:,3), LineWidth=2, Color='k');
            xlabel('ms'), xlim([20 80]), xticks(25:25:75), xticklabels(-25:25:75);
            ylabel([]), yticks([]), yticklabels([]);
            box off
        
            c = colorbar(tz,"eastoutside","Box","off");
                c.Label.String = 'm/s^2';
    
        % Plot random trial 3-axis raw data, with overlaid Acc Magnitude, around Event
        r = randperm(length(gyrData_trial),10);
        t1 = nexttile(16, [2 6]);
            add_y = 1;
            for i = 1:length(r)
                val_p_scaled = normalize(magnitude_vector_trial{r(i)}, 'scale');
                plot(add_y + gyrData_trial{r(i)}); hold on
                plot(add_y + val_p_scaled, 'LineWidth', .2 , Color=[.7 .7 .7]);
                for w = 1:length(valid_peaks{r(i),2})
                    peakt = valid_peaks{r(i),2}(w)*opt.fs + opt.ev_minus*opt.fs;
                    text(peakt, add_y + 10, "*", FontSize=20)
                end
                add_y = add_y + 30;
            end
                xline(opt.ev_minus*opt.fs, '-.', LineWidth=2, FontSize=12);
                plotops = struct('xlabel', {''}, 'ylabel', {'Trial# & Accel values (m/s^2)'}, ...
                             'xticks',  0:opt.fs:(opt.ev_minus+opt.ev_plus)*opt.fs,  'yticks', 1:30:331, ...
                             'xticklabels', -opt.ev_minus:1:opt.ev_plus, 'yticklabels', r);
                prettify(plotops);
                ylim([-15 301]);
            title('Accel. & Magnitude Traces'); box off
    
        % Plot valid pecks raster around event of interest
        str1 = sprintf('%s trials', opt.limitEvents{1});
        t2 = nexttile(28, [1 6]);
            plotRaster(valid_peaks(:,2), 1);
                plotops = struct('xlabel', {''}, 'ylabel', {str1}, ...
                                 'xticks', -opt.ev_minus:1:opt.ev_plus,  'yticks', 0:floor(size(valid_peaks(:,2),1)/2):size(valid_peaks(:,2),1), ...
                                 'xticklabels', -opt.ev_minus:1:opt.ev_plus, 'yticklabels', 0:floor(size(valid_peaks(:,2),1)/2):size(valid_peaks(:,2),1));
                prettify(plotops);
            xline(t2, 0, '-.', LineWidth=2, FontSize=12);
            title('Valid pecks raster and rate')
    
        % Plot perievent peck rate
        t3 = nexttile(34, [1 6]);
            plotPSTH(valid_aligned_ms, opt.stpSize, opt.binSize, [-opt.ev_minus*opt.fs opt.ev_plus*opt.fs], opt.fs);
                plotops = struct('xlabel', {'time (s)'}, 'ylabel', {'Pecks/s'}, ...
                                 'xticks', 0:(opt.binSize/opt.stpSize)*(opt.fs/opt.binSize):(opt.ev_minus+opt.ev_plus)*(opt.binSize/opt.stpSize)*(opt.fs/opt.binSize), ...
                                 'yticks', 0:1:3, ...
                                 'xticklabels', -opt.ev_minus:1:opt.ev_plus, ...
                                 'yticklabels', 0:1:3);
                prettify(plotops);
                xline(t3, opt.ev_minus*(opt.binSize/opt.stpSize)*(opt.fs/opt.binSize),'-.', opt.limitEvents{1}, LineWidth=2, FontSize=16);
    
        % Set texts and font sizes
        text(tx, 21, wf_count*.9, '- Av. axis trace', FontSize=12, FontWeight='bold')
    
        set(t0, FontSize=14, LineWidth=1.5, TickDir='out');
        ylim(t0, [-20 20])
    
        set(tx, YDir='normal', FontSize=12, LineWidth=1.5, TickDir='in');
        set(ty, YDir='normal', FontSize=12, LineWidth=1.5, TickDir='in');
        set(tz, YDir='normal', FontSize=12, LineWidth=1.5, TickDir='in');
        
        set(t1, FontSize=14, LineWidth=1.5, TickDir='out');
        xlim(t1, [0 (opt.ev_minus+opt.ev_plus)*opt.fs]);
    
        set(t2, FontSize=14, LineWidth=1.5, TickDir='out');
        xlim(t2, [-opt.ev_minus opt.ev_plus]);
        
        set(t3, FontSize=14, LineWidth=1.5, TickDir='out');
        xlim(t3, [0 (opt.ev_minus+opt.ev_plus)*(opt.binSize/opt.stpSize)*(opt.fs/opt.binSize)]);
        ylim(t3, [0 3])
    
    exportgraphics(fig, opt.filename, 'ContentType', 'vector');
    
end
end

function [ax, ay, az, gx, gy, gz, info] = alignAxes(data, opt)
    axr = data.acc.X(:);
    ayr = data.acc.Y(:);
    azr = data.acc.Z(:);

    gxr = data.gyr.X(:);
    gyr = data.gyr.Y(:);
    gzr = data.gyr.Z(:);

    info = struct();
    info.mode = opt.AlignMode;
    
    if strcmpi(opt.AlignMode,'mpu9250_nedlike')
        % Common MPU-9250 reconciliation to a mag/NED-like convention:
        % [x;y;z]_world = [ y; x; -z ]_raw  (swap x/y, flip z)
        A2W = [0 1 0;
               1 0 0;
               0 0 -1];
        G2W = A2W;
        info.A2W = A2W;
        info.G2W = G2W;
    else
        A2W = eye(3);
        G2W = eye(3);
        info.A2W = A2W;
        info.G2W = G2W;
    end
    
    aW = (A2W * [axr ayr azr]')';
    gW = (G2W * [gxr gyr gzr]')';
        
    ax = aW(:,1);
    ay = aW(:,2);
    az = aW(:,3);

    gx = gW(:,1); 
    gy = gW(:,2);
    gz = gW(:,3);
end
