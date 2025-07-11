function [data] = Deuteron2mat(opt)
% The function will extract data from Deuteron files, rearrange it as needed
% and filter one channel at a time. For FieldTrip, so far, we want the
% lowpass signal for LFP studies and we want to downsample to reduce
% the amount of data. We convert this into a flat .mat file that will be
% feeded into 'mat2FieldTrip'.
%
% Version 31.10.2023 (Jesus)

%% Filter preparations
% Gather info to create and apply the lowpass filter
opt.dwnsmplRate = opt.sampleRate/32; % Matches INTAN's 32x downsample factor.

%% DT2 Format
if strcmp(opt.ext, 'DT2')
    disp('Format is FLAT. DEPRECATING')
    
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

%% DF1 Format
elseif strcmp(opt.ext, 'DF1') % opt.ext = DF1
    disp('Format is BLOCK.')

    stream   = 1;
    data_tmp = int16([]);

    for i = 1:length(opt.myFiles)
        fprintf('Obtaining data from Deuteron files %d out of %d. \n',i,length(opt.myFiles));
        
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name), 'r');
            tempdata = Deuteron_extractData(stream, fid, opt);
        fclose(fid);

        % Convert ADC steps into microvolts, so conversion to int16 is
        % possible without loss.
        tempdata = int16((opt.voltageResolution * (tempdata - opt.offset)) * 1000000);

        data_tmp = [data_tmp; tempdata];
    end

    % Reshape to sort as channels x samples.
    data_tmp = reshape(data_tmp, opt.numChannels, []);
end 
clear tempdata fid

if opt.set_filter
    % To keep memory usage low, we proceed in a channel by channels basis
    data_mat = int16([]);
    
    disp('Filtering and downsampling.')
        for i=1:opt.numChannels
        
            % Proceed with filter
            [tmp, ~, ~] = bandFilter(double(data_tmp(i,:)), [], opt.lowpass, opt.sampleRate);
            
            % Proceed with downsampling
            [data_mat(i,:), ~, ~] = downsampleVolt(tmp, opt.sampleRate, opt.dwnsmplRate);
        end
else
    data_mat = data_tmp;
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