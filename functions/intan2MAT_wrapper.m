function [data, sessions] = intan2MAT_wrapper(sessions, ss, opt)
%% If not using NWB we need a mix of INTAN file reading tools 
% to bring data into MATLAB. Detailed description here.
% It def needs to know if we have one file per channel or
% one file per type.
%
% Version 01.03.2023 Jesus

if ~isfield(opt,'PathRaw'),           opt.PathRaw           = pwd;                                 end
if ~isfield(opt,'FolderProcDataMat'), opt.FolderProcDataMat = sessions.info{ss}.savefolder;        end
if ~isfield(opt,'SavFileName'),       opt.SavFileName       = sessions.list(ss).name;              end

% Check for existing 'continous_FT.mat' files
saveFolder = opt.FolderProcDataMat;
cd(saveFolder);
files = dir([saveFolder + '\*continous_FT.mat']); 

% If there is none, proceed
if isempty(files)
    % Warn about file being process.
    disp('- Will convert session to pseudo-FT format.');
    cd(fullfile(sessions.folder,sessions.list(ss).name));
    nfiles = length(sessions.info{ss}.files);

    % Set default parameters for lowpass filter. Decide to leave here or input
    % it from main script? 0-250 should be good ennough.
    filt.freqBands = [0 250];
    Orig_sR   = sessions.info{ss}.amplifier_sample_rate;
    New_sR    = Orig_sR / 32;
    
    %% Only one file. The raw will be big. Will proceed by filtering and 
    % downsampling one channel at a time. Will remove the just processed channels 
    % as it progresses, to reduce memory usage.
    if nfiles == 1
        % Open file
        disp('- One file is being opened.')
        fid = fopen(sessions.info{ss}.files.name, 'r');
        
        % read voltage data according to INTAN
        disp('- All channels are being readed from single file.')
        tmp.volt = fread(fid, [sessions.info{ss}.nchannels inf], 'int16');
        fclose(fid);
        
        % 30KHz nSamples
        sessions.info{ss}.nSamples = length(tmp.volt);

        % If it is 'amp', needs filter and downsampling
        if strcmp(sessions.info{ss}.bandpass, 'amp')
            for ff = 1:sessions.info{ss}.nchannels
                fprintf('- Processing channel %d of %d.\n', ff, sessions.info{ss}.nchannels);
                
                % Proceed with filter
                [tmp.volt(1,:),filt.a,filt.b] = bandFilter(tmp.volt(1,:),[],filt.freqBands,Orig_sR);

                % Proceed with downsampling
                [tmp.v(ff,:),~,downsmpFactor] = downsampleVolt(tmp.volt(1,:),Orig_sR,New_sR);
                tmp.volt(1,:) = [];
            end

            % Keep record of used parameters in filtering and downsampling if any
            sessions.info{ss}.lowpass_filt = filt;
            sessions.info{ss}.lowpass_downsample = downsmpFactor;
            sessions.info{ss}.lowpass_sample_rate = New_sR;

        % If it is 'low', data is simply passed on
        elseif strcmp(sessions.info{ss}.bandpass, 'low')
            tmp.v = tmp.volt;
            tmp.volt = [];
        end

    elseif nfiles > 1
    %% Many files 
        disp('- Multiple files will be opened and readed, one by one.')
        for ff = 1:nfiles
            % Open file by file
            fid = fopen(sessions.info{ss}.files(ff).name, 'r');
    
            % Read voltage data according to INTAN and attach
            tmp.volt = fread(fid, [1 inf], 'int16');
            fclose(fid);
            
            % From wide data, we need to filter/downsample the data
            if strcmp(sessions.info{ss}.files(1).name(1:3), 'amp')
                fprintf('- Processing file %d of %d.\n', ff, nfiles);
                % Use filter and downsampling functions (designed for ETALO).
                % DETAILED explanation inside functions.
    
                % nSamples of 30KHz data, just once
                if ff==1
                    sessions.info{ss}.nSamples = length(tmp.volt);
                end

                % Proceed with filter
                [tmp.volt,filt.a,filt.b] = bandFilter(tmp.volt,[],filt.freqBands,Orig_sR);
    
                % Proceed with downsampling
                [tmp.v(ff,:),~,downsmpFactor] = downsampleVolt(tmp.volt,Orig_sR,New_sR);
    
                % Keep record of used parameters in filtering and
                % downsampling, just once
                if ff==1
                    sessions.info{ss}.lowpass_filt = filt;
                    sessions.info{ss}.lowpass_downsample = downsmpFactor;
                    sessions.info{ss}.lowpass_sample_rate = New_sR;
                end

            elseif strcmp(sessions.info{ss}.files(1).name(1:3), 'low')
                % If band is 'lowpass' already, directly pass it on
                tmp.v(ff,:) = tmp.volt;
            end

            % Clear original raw voltage data
            tmp.volt = [];
        end
    end
    
    % Convert to microvolts
    volt = tmp.v * 0.195;

    % clean up
    clear fid tmp
    
    %% Get time series
    % There is only one file for the time series: 'time.dat' at 30KHz. We don't want it.
    % For LFP, having the voltage series already, we are going to use the number
    % of samples there to create our own time-series. 
    
    % Number of samples of resulting downsampled data
    sessions.info{ss}.lowpass_nSamples = length(volt);

    % We create a time-vector in samples and divide it by the sampling rate.
    time = (1:sessions.info{ss}.lowpass_nSamples) / sessions.info{ss}.lowpass_sample_rate;
            
    % Easy access to total recording time.
    sessions.info{ss}.recording_time = sessions.info{ss}.lowpass_nSamples / sessions.info{ss}.lowpass_sample_rate;
    
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
    
    % Starting with labels as they have been extracted from the INTAN header
    for i=1:sessions.info{ss}.nchannels
        data.label{i,1} = convertStringsToChars(sessions.info{ss}.INTAN_hdr.amplifier_channels(i).native_channel_name);
    end
    
    % The only trial contains all channels*time info                
    data.trial{1}        = volt;
    
    % The only trial is the whole time-series
    data.time{1}         = time;
    
    % The trial starts at 0 and ends at last sample
    data.sampleinfo(1,:) = [1 sessions.info{ss}.lowpass_nSamples];
    
    disp('- Done.')

else

    % If there is any, exit the function and continue
    % Warn about existing pseudo-FT files
    disp('- A file in the FT format has been found for this session!');
    
    % Data output is empty
    data = [];
end   


end