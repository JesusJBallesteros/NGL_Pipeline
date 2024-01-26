function [data] = Deuteron2Fieldtrip(opt)
% The function will extract data from Deuteron files, rearrange it as needed
% and filter one channel at a time. For FieldTrip, so far, we want the
% lowpass signal for LFP studies and we want to downsample to reduce
% the amount of data. We convert this into a flat .mat file that will be
% feeded into 'mat2FieldTrip'.

% Jesus 25.01.2024

%% Check existence of a FieldTrip file.
% If existing, load it instead and return to main script
if isfile(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_continous_FT.mat'))) || ...
   isfile(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_tparsed_FT.mat')))
    disp('A Fieldtrip-formatted file found in this directory, skipping.')
% 
%     try     data = load(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_continous_FT.mat')));
%     catch,  data = load(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_tparsed_FT.mat')));
%     end
% 
    data = [];
    return
end

%% Filter preparations
% Gather info to create and apply the lowpass filter
opt.dwnsmplRate = opt.sampleRate/32; % Matches INTAN's 32x downsample factor.

%% Proceed for DT2 Format. Deprecating.
if strcmp(opt.ext, 'DT2')
    disp('Converting DT2 files to pseudo-FieldTrip.')
    
    % Initiate matrix and sample index.
    data_tmp = int16([]);
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
        data_tmp(:,indexPos+1:indexPos+nSamples) = tempdata;
       
        % Get to next starting sample
        indexPos = indexPos+nSamples;
    end
end
%% Proceed for DF1 Format.
if strcmp(opt.ext, 'DF1') % opt.ext = DF1
    data_tmp = int16([]); % Create empty variable to store all data (do not pre-allocate the whole matrix)
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

        data_tmp = [data_tmp; tempdata];
    end
    clear tempdata fid

    % Reshape to sort as channels x samples.
    data_tmp = reshape(data_tmp, opt.numChannels, []);
end 

%% Let's always filter
data_mat = int16([]);
txt = sprintf('Filtering between %d and %d Hz. It may take a moment.\n', opt.lowpass(1), opt.lowpass(2));
fprintf(txt);

% To keep memory usage low, we proceed in a channel by channels basis
for i=1:opt.numChannels
    % Proceed with filter
    [tmp, ~, ~] = bandFilter(double(data_tmp(i,:)), [], opt.lowpass, opt.sampleRate);
    
    % Proceed with downsampling
    [data_mat(i,:), ~, ~] = downsampleVolt(tmp, opt.sampleRate, opt.dwnsmplRate);
end
clear data_temp

%% Get time series
% For LFP, having the voltage series already, we are going to use the number
% of samples there to create our own time-series. 

% Number of samples of resulting downsampled data
nSamples = length(data_mat);

% We create a time-vector in samples and divide it by the sampling rate.
time = (1:nSamples) / opt.dwnsmplRate; % results in sec.
        
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
    nch = sprintf('%03d', opt.channelOrder(i));
    data.label{i,1} = nch;
end

% The only trial contains all channels*time info                
data.trial{1}        = data_mat;

% The only trial is the whole time-series
data.time{1}         = time;

% The trial starts at 0 and ends at last sample
data.sampleinfo(1,:) = [1 nSamples];

disp('- Done.')
end