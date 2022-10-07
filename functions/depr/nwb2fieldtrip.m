function [data, cored] = nwb2fieldtrip(sessions,cored)
% Create FieldTrip compatible files
%   Detailed explanation goes here

%% Start by creating the NWB core
disp('- Creating necessary dependencies:')

if cored==0
    cd functions\toolboxes\matnwb;
    addpath(genpath(pwd));
    generateCore()
    cored = 1;
end

% Navigate to session
cd(strcat(sessions.folder, '\', sessions.name));

% Filename and path
nwbFile = dir('*.nwb');

% Show schema version of the file. If this does not match the installed version,
% choose the schema closest to yours from NWB release site.
disp(util.getSchemaVersion(nwbFile.name))

%% Load data in nwb format
nwb = nwbRead(nwbFile.name);
disp(nwb)

%% Obtain header info and LFP data 
try
    disp('- Obtaining header info.')
    % Will try to find keys for electricalSeries. If there are more than
    % one will try to keep the lowpass. If there is only one, will go for
    % the 'aplifier' version of it. Therefore we have to address somehow
    % which one are we trying to get here.
    hdr = ft_read_header_INTAN(nwbFile.name);
    
    % For LFP files, it will find some 'lowpass'
    % electricalSeries, we should have account for this using the 'low'
    % keyword at the input?
    if strcmp(sessions.info.bandpass,'low')
        if isempty(hdr.Fs) 
            hdr.Fs = sessions.info.lowpass_sample_rate;
        end
    
        % Lowpass is probably downsampled. We adjust nSamples to convert
        % from samples to seconds
        hdr.nSamples = hdr.nSamples / sessions.info.lowpass_downsample;

    elseif strcmp(sessions.info.bandpass,'amp')
        if isempty(hdr.Fs) 
            hdr.Fs = sessions.info.amplifier_sample_rate;
        end

        % Amplifier is not downsampled. No need to modify nSamples
        % But we want to filter and downsample later, so we set it up.
        freqBands = [0 250];
        Orig_Fs = sessions.info.amplifier_sample_rate;
        smpRate = Orig_Fs / 32; % To equalize to any other INTAN lowpass
    end

catch ME
    disp('- Could not load in hdr information')
    rethrow(ME)
end

try
    disp('- Reading data. This will take a while.'); tic
    tmp.volt{1} = ft_read_data_INTAN(nwbFile.name, 'begsample', 1, 'endsample', hdr.nSamples);
    toc

    % From wide data, we need to filter/downsample the data
    if strcmp(sessions.info.bandpass{1}, 'amp')
        % Use filter and downsampling functions (designed for ETALO).
        % DETAILED explanation inside functions.
        for ch=1:hdr.nChans
            % Proceed with filter
            tofilt = double(tmp.volt{1}(ch,:));
            [tmp.volt{2},filt.a,filt.b] = bandFilter(tofilt,[],freqBands,smpRate);
            tofilt = [];
    
            % Proceed with downsampling
            [volt(ch,:),~,downsmpFactor] = downsampleVolt(tmp.volt{2},Orig_Fs,smpRate);
            tmp.volt{2} = [];    
        end
        
        % Remove raw data
        clear tmp tofilt

        % Keep record of used parameters
        sessions.info.lowpass_filt = filt;
        sessions.info.lowpass_downsample = downsmpFactor;
        sessions.info.lowpass_sample_rate = smpRate;

        % Convert to microvolts
        volt = volt * 0.195;

    else
        % If band is not 'amp', it must be 'lowpass' already (FOR NOW)
        % Convert to microvolts
        volt = tmp.volt{1} * 0.195;

        % Remove raw data
        clear tmp
    end    

catch ME
    disp('- Could not extract LFP data')
    rethrow(ME)
end

%% Circumvent native FT reading functions and create own FT-like data structure
% data.label      % cell-array containing strings, Nchan*1
% data.fsample    % sampling frequency in Hz, single number
% data.trial      % cell-array containing a data matrix for each
%                 % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
% data.time       % cell-array containing a time axis for each
%                 % trial (1*Ntrial), each time axis is a 1*Nsamples vector
% data.trialinfo  % this field is optional, but can be used to store
%                 % trial-specific information, such as condition numbers,
%                 % reaction times, correct responses etc. The dimensionality
%                 % is Ntrial*M, where M is an arbitrary number of columns.
% data.sampleinfo % optional array (Ntrial*2) containing the start and end
%                 % sample of each trial

disp('- Creating pseudo-FieldTrip structure...')
for i=1:hdr.nChans 
    data.label{i,1}      = convertStringsToChars(sprintf('A-%03d', str2double(hdr.label{i})));
end

% For 'single-trial' data
data.trial{1}      = volt;
data.time{1}       = (1:length(volt)) / sessions.info.lowpass_sample_rate;
data.sampleinfo    = [1 length(volt)];
disp('- Done.')
end