function [data] = intan2MAT_wrapper(sessions, opt)
% intan2MAT_wrapper  Read and preprocess INTAN data into pseudo-FieldTrip struct.
%
% PURPOSE:
%   Reads INTAN amp*.dat files, applies the LFP preprocessing chain (scale - CAR -
%   detrend - low-pass - line-noise filter - downsample), and assembles the
%   result into a pseudo-FieldTrip data struct. This is the input to
%   MAT2FieldTrip, which concludes the FieldTrip format.
%
% USAGE:
%   data = intan2MAT_wrapper(sessions, opt)
%   Called from INTAN_PipelineWrapper.
%
% INPUTS:
%   sessions  - struct from input.sessions(x), must contain:
%                 .info.files               dir-struct of amp*.dat files
%                 .info.nChannels           electrode count
%                 .info.nfiles              number of .dat files
%                 .info.INTAN_hdr           full RHD header (for channel labels)
%                 .info.amplifier_sample_rate (Hz)
%   opt       - options struct; relevant fields:
%                 .dwnsmplRate   target LFP sample rate (Hz); [] -> 937.5 Hz
%                 .lowpassFT     LFP low-pass cutoff (Hz, default 250)
%                 .linefilter    line-noise frequency (Hz, 0 = off)
%                 .CAR           common-average referencing flag (0 = off)
%
% OUTPUT:
%   data  - pseudo-FieldTrip struct:
%             .label      {nChannels × 1 cell} channel name strings
%             .trial      {1 × 1 cell} [nChannels × nSamples double] LFP data (µV)
%             .time       {1 × 1 cell} [1 × nSamples double] time vector (s, 0-indexed)
%             .sampleinfo [1 2] = [1, nSamples]
%
% PREPROCESSING (if opt.set_filter == 1):
%   1. Scale to uV, ×0.195
%   2. CAR (optional): ft_preproc_rereference with 'all'/'median'
%   3. Detrend: ft_preproc_detrend (removes DC offset per channel)
%   4. Low-pass: ft_preproc_lowpassfilter (Butterworth 6th order, twopass)
%   5. Line-noise: ft_preproc_bandstopfilter (opt.linefilter ±2 Hz, if > 0)
%   6. Downsample: downsampleVolt (integer factor only)
%
% NOTES:
%   - Time vector is 0-indexed: (0:N-1)/dwnsmplRate in seconds.
%   - For CAR with > 32 channels, applies per-bank (1:32, 33:64) referencing.
%   - Downsample default 937.5 Hz = 30000/32; must be integer divisor of Fs.
%
% Version 07.05.2026 (Jesus)

%% Collect parameters to proceed with file creation
% List all files (multiple or single depending on type). If No lowpass
% files found, we will use the raw data, and filtering will be applied.
disp('Will convert raw data to preprocessed pseudo-FT format.');

% opt.myFiles = dir('low*.dat'); % low files won;t need processing
% opt.myFiles = dir('amp*.dat'); % no need for re-reading
opt.myFiles = sessions.info.files; % use collected files

if ~isempty(opt.myFiles)
    opt.set_filter = 1; % it's raw, needs lowpassing and downsampling

    % Gather info to create and apply the lowpass filter
    opt.sampleRate  = sessions.info.amplifier_sample_rate;
    if isempty(opt.dwnsmplRate) % if downsample rate is set to 1, we won't downsample anything
        % The lowpass data (250Hz) can be sampled at 1KHz, or 937.5 to fit INTAN's numbers
        opt.dwnsmplRate = 937.5; % would match the INTAN lowpass files
    end
else
    error('No raw data to process found')
%     opt.set_filter = 0; % If the files are already lowpassed and downsampled
%     opt.dwnsmplRate = sessions.info.amplifier_sample_rate / sessions.info.lowpass_downsample;
end

nfiles = length(opt.myFiles);

%% Open INTAn file/s and bring to matlab temporal array
% Either at once 
if nfiles == 1
    disp('All channels are being read from single file.')
    % Read voltage data according to INTAN
    % Open file, read as 'int16' but store as double.
    fid = fopen(sessions.info.files.name, 'r');
        tmp = fread(fid, [sessions.info.nChannels inf], 'int16');
    fclose(fid);

    % Scale
    tmp = doScale(tmp);

    % Common Median referencing
    if opt.CAR, tmp = doCar(tmp, sessions);  end

    % Filtering
    if opt.set_filter, tmp = doFilters(tmp, opt); end
    
    % Downsample
    volt = doDownsample(tmp, opt);

% or channel by channel
else  
    if opt.CAR % All channels needed
        disp('Because CAR, all files will be opened one by one but treated at once.')

        % Open file by file
        for b = 1:nfiles
            % Read voltage as 'int16', store as double.
            fprintf('- Opening file %d.\n', b);
            fid = fopen(sessions.info.files(b).name, 'r');
                tmp(b,:) = fread(fid, [1 inf], 'int16');
            fclose(fid);
        end

        % Scale
        tmp = doScale(tmp);

        % Common Median referencing
        tmp = doCar(tmp, sessions);

        % Filtering
        if opt.set_filter, tmp = doFilters(tmp, opt); end

        % Downsample
        volt = doDownsample(tmp, opt);

    else
        disp('All files will be opened and treated one by one.')
        for b = 1:nfiles
            % Read voltage as 'int16', store as double.
            fprintf('- Opening file %d.\n', b);
            fid = fopen(sessions.info.files(b).name, 'r');
                tmp = fread(fid, [1 inf], 'int16');
            fclose(fid);

            % Scale
            tmp = doScale(tmp);

            % Filtering
            if opt.set_filter
                fprintf('- Filtering channel %d.\n', b);
                tmp = doFilters(tmp, opt);
            end

            % Downsample
            volt(b,:) = doDownsample(tmp, opt);
        end
        clear tmp

    end
end

%% Convert to pseudo-FieldTrip
% It's only pseudo until we run the proper FT tool to check for format and
% header info. Because we have not given any trial info so far, the data
% comes as a continous single trial. 
% data.label      % cell-array containing strings, Nchan*1
% data.trial      % cell-array containing a data matrix for each
%                 % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
% data.time       % cell-array containing a time axis for each
%                 % trial (1*Ntrial), each time axis is a 1*Nsamples vector
% data.sampleinfo % optional array (Ntrial*2) containing the start and end
%                 % sample of each trial

disp('Creating pseudo-FieldTrip structure...');
% Starting with labels as they have been extracted from the INTAN header
for i = 1:sessions.info.nChannels
    if i<=sessions.info.nfiles
        data.label{i,1} = convertStringsToChars(sessions.info.INTAN_hdr.amplifier_channels(i).native_channel_name);
    else
        continue
    end
end

% The only trial contains all channels*time info                
data.trial{1} = volt;

% The only trial is the whole time-series.
% We simply create a time-vector from samples and divide it by the sampling rate.
% data.time{1} = (1:length(volt)) / opt.dwnsmplRate; % Fixed 1-idexed, when FT expects 0-indexed
data.time{1} = (0:length(volt)-1) / opt.dwnsmplRate; % in Seconds

% Therefore, trial starts at first sample and ends at last sample
data.sampleinfo(1,:) = [1 length(volt)];

disp('Done.')
clear time volt
end

%% Helper functions
% Scale to uvolt
function tmp = doScale(tmp)
    tmp = tmp * 0.195;
end

% Re-reference channels
function tmp = doCar(tmp, sessions)
    if sessions.info.nChannels > 32 % Two banks, from two different regions. Hardcoded for 'chgDet' specific case
        % TODO generalize
        disp('Re-referencing by Common Average Referencing (CARing) 1/2.')
        [tmp(1:32,:), ~] = ft_preproc_rereference(tmp(1:32,:), 'all', 'median');
        disp('Re-referencing by Common Average Referencing (CARing) 2/2.')
        [tmp(33:64,:), ~] = ft_preproc_rereference(tmp(33:64,:), 'all', 'median');
    else
        % In principle, data from a single HS on a single region.
        disp('Re-referencing by Common Average Referencing (CARing).')
        tmp = ft_preproc_rereference(tmp, 'all', 'median');
    end

%     disp('Saving CARed file, will take a while.')
%     save(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_CARed.mat']), 'tmp', '-v7.3');
end

% Filtering, if required (preprocessing raw)
function tmp = doFilters(tmp, opt)  
    % Now, channel by channel
    for i = 1:size(tmp,1)    
        % if more than one
        if size(tmp,1) > 1
            fprintf('- Filtering channel %d of %d.\n', i, opt.numChannels);
        end

        % Detrend channel (remove DC)
        disp('Detrending...')
        tmp(i,:) = ft_preproc_detrend(tmp(i,:));

        % Lowpass filter channel (Butterwort, 6th order, back&forth)
        disp('Lowpassing...')
        [tmp(i,:), ~, ~] = ft_preproc_lowpassfilter(tmp(i,:), opt.sampleRate, opt.lowpassFT, 6, 'but', 'twopass');
                    
        % FT's bandstop filter (btw 50 +-2 Hz, Butterwort, 2nd order, back&forth)
        if opt.linefilter > 0
            disp('Line denoising...')
            [tmp(i,:), ~, ~] = ft_preproc_bandstopfilter(tmp(i,:), opt.sampleRate, [opt.linefilter-2 opt.linefilter+2], 2, 'but', 'twopass', 'split');
        end
    end
end

% Downsample. Get new time vector.
function volt = doDownsample(tmp, opt)
    disp('Downsampling.')
    [volt, ~, ~] = downsampleVolt(tmp, opt.sampleRate, opt.dwnsmplRate, 2);
    
end