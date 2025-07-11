function [data] = Deuteron2Fieldtrip(opt)
% The function will extract data from Deuteron files, rearrange it as needed
% and filter one channel at a time. For FieldTrip, so far, we want the
% lowpass signal for LFP studies and we want to downsample to reduce
% the amount of data. We convert this into a flat .mat file that will be
% feeded into 'mat2FieldTrip'.

% Jesus 11.06.2024

%% Check existence of a FieldTrip file.
% If existing, load it instead and return to main script
if isfile(fullfile(opt.trialSorted, strcat(opt.SavFileName,'_FTcont.mat')))
    disp('A Fieldtrip-formatted file found in this directory, skipping.')
    data = [];
    return
end

%% Filter preparations
% Gather info to create and apply the lowpass filter
opt.dwnsmplRate = opt.sampleRate/32; % Matches INTAN's 32x downsample factor.

%% Find pre-processed data if any
if ~isfile(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_filt_dwn.mat']))
    %% Proceed for DT2 Format. Deprecating.
    if strcmp(opt.ext, 'DT2')
        disp('Converting DT2 files to pseudo-FieldTrip.')
        
        % Initiate matrix and sample index.
        tmp = int16([]);
        indexPos  = 0;
    
        disp('Obtaining data from Deuteron files.')
        for i = 1:length(opt.myFiles)
            % Neural data points are 16 bit words
            fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
                tempdata = fread(fid, 'uint16');
            fclose(fid);
            
            % Remove trailing zeroes in last file
            if i == length(opt.myFiles)
                tempdata(tempdata==0) = [];
            end
    
            % Reshape concatenated channels to channels x samples.
            tempdata = reshape(tempdata', opt.numChannels, []);
    
            % Convert ADC steps into microvolts, so conversion to int16 is
            % possible without loss.
            tempdata = int16((opt.voltageResolution * (tempdata - opt.offset)) * 1000000);
            
            % Get nSamples coming from this file. Shouls stay constant. 
            nSamples = size(tempdata,2);
    
            % Collect chunks into full matrix for further treatment.
            tmp(:,indexPos+1:indexPos+nSamples) = tempdata;
           
            % Get to next starting sample
            indexPos = indexPos+nSamples;
        end
    end
    %% Proceed for DF1 Format.
    if strcmp(opt.ext, 'DF1') % opt.ext = DF1
        tmp = int16([]); % Create empty variable to store all data (do not pre-allocate the whole matrix)
        opt.stream   = 1;     % Pass variable to read continuous neural signals.
    
        % Open each channel file and read it, resize data to fit Kilosort expectations,
        % and concatenate it to the matrix consecutively.
        disp('Converting DF1 files to pseudo-FieldTrip. It may take a moment.')
        for i = 1:length(opt.myFiles)
            fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name), 'r'); % Open file.
                tempdata = Deuteron_extractData(fid, opt); % read data with Deuteron's script.
            fclose(fid); % Close file
    
            % Convert ADC bit steps into microvolts to convert to int16 without loss.
            tempdata = int16((opt.voltageResolution * (tempdata - opt.offset)) * 1000000);
    
            tmp = [tmp; tempdata];
        end
        clear tempdata fid
    
        % Reshape to sort as channels x samples.
        tmp = reshape(tmp, opt.numChannels, []);
    end 

    %% Re-reference channels
    if opt.CAR == 2
        % Check for pre-treated files
        if ~isfile(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_CARed.mat']))
            % In principle, data from a single HS on a single region.
            disp('Re-referencing by Common Average Referencing (CARing).')
            tmp = ft_preproc_rereference(tmp, 'all', 'median');
    
            disp('Saving CARed file, will take a while.')
            save(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_CARed.mat']), 'tmp', '-v7.3');
        else
            disp('CARed data existing, loading.')
            load(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_CARed.mat']));
        end
    end

    %% Let's always filter
    % Now, channel by channel to keep memory usage low
    % TODO: paralellize?
    for b = 1:opt.numChannels
        fprintf('Channel %d of %d.\n', b, opt.numChannels);
        
        % Detrend channel (remove DC)
        disp('Detrending...')
        tmp(b,:) = ft_preproc_detrend(tmp(b,:));
    
        % Lowpass filter channel (Butterwort, 6th order, back&forth)
        disp('Lowpassing...')
        [tmp(b,:), ~, ~] = ft_preproc_lowpassfilter(tmp(b,:), opt.sampleRate, opt.lowpass, 6, 'but', 'twopass');
                    
        % FT's bandstop filter (btw 50 +-2 Hz, Butterwort, 2nd order, back&forth)
        if opt.linefilter > 0
            disp('Line denoising...')
            [tmp(b,:), ~, ~] = ft_preproc_bandstopfilter(tmp(b,:), opt.sampleRate, [opt.linefilter-2 opt.linefilter+2], 2, 'but', 'twopass', 'split');
        end
    end

    % Downsample
    disp('Downsampling.')
    [volt, ~, ~] = downsampleVolt(tmp, opt.sampleRate, opt.dwnsmplRate, 2);
    clear tmp

    % Save depending on previous treatment
    disp('Saving Filtered and downsampled data.')
    save(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_filt_dwn.mat']), 'volt', '-v7.3');
    
else
    disp('Filtered and downsampled data already existing, loading.')
    load(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_filt_dwn.mat']));
end

%% Get time series
% We simply create a time-vector from samples and divide it by the sampling rate.
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

disp('- Creating pseudo-FieldTrip structure...')
% Starting with labels as they have been extracted from Deuteron
for i=1:opt.numChannels
    nch = sprintf('%03d', i);
    data.label{i,1} = nch;
end

% The only trial contains all channels*time info                
data.trial{1}        = volt;

% The only trial is the whole time-series
data.time{1}         = time;

% The trial starts at t=0 and ends at t=t(end)
data.sampleinfo(1,:) = [1 length(time)];

disp('- Done.')
end