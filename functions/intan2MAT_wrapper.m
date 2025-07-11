function [data] = intan2MAT_wrapper(sessions, opt)
% Obtain data files and read whole or file by file. In any case, process in
% a channel by channel basis through: detrending, lowpass filtering, line
% noise filtering, re-referencing (if requested) and downsampling to 1KHz.
% All that gets into a single variable and given a pseudo_FT format, ready
% to get a 'continuous', 'trial-parsed' or both treatments in the next
% step.
%
% Version 16.06.2023 Jesus

%% Collect parameters to proceed with file creation
% List all files (multiple or single depending on type). If No lowpass
% files found, we will use the raw data, and filtering will be applied.
disp('Will convert raw data to preprocessed pseudo-FT format.');
% opt.myFiles = dir('low*.dat');
opt.myFiles = dir('amp*.dat');

if ~isempty(opt.myFiles)
    opt.set_filter = 1; % it's raw, needs lowpassing and downsampling

    % Gather info to create and apply the lowpass filter
    opt.sampleRate  = sessions.info.amplifier_sample_rate;
    opt.dwnsmplRate = 937.5; % would match the INTAN lowpass files

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
    if opt.CAR, tmp = doCar(tmp);  end

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
data.time{1} = (1:length(volt)) / opt.dwnsmplRate; % in Seconds

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
        [tmp(i,:), ~, ~] = ft_preproc_lowpassfilter(tmp(i,:), opt.sampleRate, opt.lowpass, 6, 'but', 'twopass');
                    
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