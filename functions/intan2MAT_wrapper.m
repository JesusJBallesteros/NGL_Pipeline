function [data] = intan2mat_wrapper(sessions, ss, varargin)
% We need a mix of INTAN file reading tools to bring data into MATLAB.
% It def needs to know if we have one file per channel or one file per type.
% To complete description ...
%
% Dependencies: 'bandFilter'
%               'downsampleVolt'
%
% Version 05.04.2023 Jesus

%% Collect parameters to proceed with file creation
% List all files (multiple or single depending on type). If No lowpass
% files found, we will use the raw data, and filtering will be applied.
disp('Will convert session to pseudo-FT format.');
opt.myFiles = dir('low*.dat');

if isempty(opt.myFiles)
    opt.myFiles = dir('amp*.dat');
    opt.set_filter = 1;

    % Gather info to create and apply the lowpass filter
    opt.sampleRate  = sessions.info(ss).amplifier_sample_rate;
    opt.dwnsmplRate = 937.5; % Matches INTAN's 32x downsample factor
else
    opt.set_filter = 0;
    opt.dwnsmplRate = sessions.info(ss).amplifier_sample_rate / sessions.info(ss).lowpass_downsample;
end

nfiles = length(opt.myFiles);

%% Only one file. The raw will be big. Will proceed by filtering and 
% downsampling one channel at a time.
if nfiles == 1
    disp('All channels are being readed from single file.')
    
    % Read voltage data according to INTAN
    % Open file, read as 'int16' but store as double.
    fid = fopen(sessions.info(ss).files.name, 'r');
        tmp = fread(fid, [sessions.info(ss).nchannels inf], 'int16');
    fclose(fid);

    % Convert to microvolts
    tmp = tmp * 0.195;
    
    % If filtering is required (meaning, we are dealing with 'amp' files)
    if opt.set_filter
        % Go channel by channel.
        for b = 1:sessions.info(ss).nchannels
            fprintf('- Filtering channel %d of %d.\n', b, opt.numChannels);
            
            % Proceed with filter. 'bandFilter' likes double precision.
            % Check function for more options, arguments and doings.
            [tmp(b,:), ~, ~] = bandFilter(tmp(b,:), [], opt.lowpass, opt.sampleRate);

            % Proceed with downsampling. Down to a fix 937.5 Hz, matching
            % the possible 'low' files from INTAN, if downsampling at that
            % time selected.
            [volt(b,:), ~, ~] = downsampleVolt(tmp(b,:), opt.sampleRate, opt.dwnsmplRate);
            
        end

    % If filtering is not required, just pass the data through.
    else
        volt = tmp;
    end

elseif nfiles > 1
    %% Many files 
    disp('Multiple files will be opened and readed, one by one.')
    % Open file by file
    for b = 1:nfiles
        % Read voltage as 'int16', store as double.
        fid = fopen(sessions.info(ss).files(b).name, 'r');
            tmp = fread(fid, [1 inf], 'int16');
        fclose(fid);
        
        % Convert to microvolts
        tmp = tmp * 0.195;

        % If filtering is required.
        if opt.set_filter
            fprintf('Processing file %d of %d.\n', b, nfiles);
            % Proceed with filter. Likes double precision, check function
            % for more about it.
            [tmp, ~, ~] = bandFilter(tmp, [], opt.lowpass, opt.sampleRate);

            % Proceed with downsampling
            [volt(b,:), ~, ~] = downsampleVolt(tmp, opt.sampleRate, opt.dwnsmplRate);

        else
            % Filtering is not required
            volt(b,:) = tmp;
        end
    end
end

clear fid tmp tmp_filt

%% Get time series
% There is only one file for the time series: 'time.dat' at 30KHz. We don't want it.
% For LFP, having the voltage series already, we are going to use the number
% of samples there to create our own time-series. 

% We create a time-vector in samples and divide it by the sampling rate.
time = (1:length(volt)) / opt.dwnsmplRate; % in Seconds
        
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
for i = 1:sessions.info(ss).nchannels
    data.label{i,1} = convertStringsToChars(sessions.info(ss).INTAN_hdr.amplifier_channels(i).native_channel_name);
end

% The only trial contains all channels*time info                
data.trial{1}        = volt;

% The only trial is the whole time-series
data.time{1}         = time;

% The trial starts at 0 and ends at last sample
data.sampleinfo(1,:) = [1 length(volt)];

disp('Done.')

end