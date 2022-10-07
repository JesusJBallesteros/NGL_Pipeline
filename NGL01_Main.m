%% Jesus' Pipeline to read INTAN continous data
% Will read INTAN data, from selected sessions for a given animal, to MATLAB.
% For that, we can use two methods:
% First one involves a transformation INTAN-to-NWB format (HDF5 format) and 
%  then NWB-to-FieldTrip. This could, in priciple, keep metadata from INTAN 
%  files and access data in chunks. However the INTAN-to-NWB is kinda buggy,
%  and involves the installation and use of python (wrapped here).
% The second makes use of INTAN file reading into MATLAB and is much
%  straight forward. However all metadata would need to be added 'manually'
%  and the MAT files are supposed to be readed as a whole.
% Next, it should give format to fit in any FieldTrip pipeline.
% Should be able to recognize INTAN channel-per-file and type-per-file
% formats. For now, it will be limited to create continous data or
% fix-length trials, not using eventcodes.
% Eventually, merge with Sara's spike pipeline for Kilosort.

% TODO 
%       set several bandpass to read more than one sort of data.
%       There seems to be an ERROR on 2nd and following runs of the NWB functionalities.
%       
%       

% Last modified By Jesus J. Ballesteros 06.10.2022

%% 00. Needed input
% TES: {'20220609' '20220809' '20220922'}
% FRN: {'20181029' '20181105' '20181106'}
% FAT: {'20220909'}
% 427: {'20220929'}
input.mainfolder = pwd;                 % string. Main Code folder. Gets current by default
input.datafolder = "D:\Experiments\";   % chr array. Main data folder

input.animal   = '427';   % string 'TES', 'FAT' , 'FRN', '427' ...
input.dates    = 'all';   % cell array {'yyyymmdd' ...} or string 'all'
input.test_ch  = 2:2:12;  % int array If ~empty, plot snippet signal for channels
input.useNWB   = 0;       % int 1/0 for use/not use of NWB.

% TODO: 'high', 'spike' not yet available
input.bandpass = {'low' 'amp'};  % cell array up to {'low' 'amp' 'high' 'spike'}

%% 00. Dependencies and defaults
cd(input.mainfolder)
addpath functions\
addpath functions\toolboxes\fieldtrip_light
addpath(genpath('functions\toolboxes\chronux_2_12')) 

ft_defaults

% Adds to useNWB input to indicate that core has not been generated on first run. 
input.useNWB(2) = 0;

%% 01. Find sessions
% Get some/all sessions in animal folder
if iscell(input.dates) % input is cell array of dates
    input.dates = input.datafolder + input.animal + "\" + input.animal + "_" + input.dates(:);
    [sessions.folder,sessions.name,~] = fileparts(input.dates);
    sessions.folder = unique(sessions.folder);

elseif strcmp(input.dates, 'all') % input is 'all'
    sessions.folder = input.datafolder + input.animal;
    sessions.name   = dir(sessions.folder + '\' + input.animal + '*');
    sessions.name   = {sessions.name(:).name}';
end

% Count sessions
sessions.nSessions = length(sessions.name);

% Read into MATLAB
for ss=1:sessions.nSessions
    % Navigate to session folder.
    cd(strcat(sessions.folder,'\',sessions.name(ss)));
    
    %% 02. Find out INTAN settings and header file. Extract info
    % Only exists when 'filepertype' or 'fileperch'
    if isfile('info.rhd') 
        % Get more info from main INTAN header file. 
        %  Uses a modified Intan funtion to output info
        [sessions.info.INTAN_hdr] = mod_read_Intan_RHD2000_file('info.rhd');

        % Number of channels.
        sessions.info.nchannels = length(sessions.info.INTAN_hdr.amplifier_channels);

        % Here we look for files of each of the bandpass selected as input. 
        %  'amp' should always exist. Would be used as ultimate source of
        %  data if others do not. If others are requested AND exist, the 
        %  loop breaks and takes the last indexed file list.
        for i=1:numel(input.bandpass)
            files = dir([input.bandpass{i},'*.dat']);
             if ~isempty(files) && strcmp(input.bandpass{i},'low')
                sessions.info.files = dir('low*.dat');
                break
                % TODO: this is not ready yet
%              elseif ~isempty(files) && strcmp(input.bandpass{i},'high')
%                 sessions.info.files = dir('high*.dat');
%                 break
%              elseif ~isempty(files) && strcmp(input.bandpass{i},'spike')
%                 sessions.info.files = dir('spike*.dat');
%                 break
             else
                sessions.info.files = dir('amp*.dat');
             end
        end

        % The final bandpass to use.
        sessions.info.bandpass = input.bandpass{i};

        % How many files exist for this sessions.
        if length(sessions.info.files)>1
            % If there are many files, is 'fileperch'
            sessions.info.fileformat = 'fileperch';
        else
            % If there is only one, is 'filepertype'
            sessions.info.fileformat = 'filepertype';
        end
        % If we have several types ('low' & 'amp') with file per channel
        % format, 'fileperch' applies.

    else % TODO 'Intanformat'. Very low priority
        sessions.info.fileformat = 'Intanformat';
        disp('- Support for Intanformat files not yet ready.')
        return
    end
    
    % If available, parse info from 'settings.xml' (for easy access).
    %  Contains metadata that may be worth to keep. Some data is relocated
    %  to have easier access.
    if isfile('settings.xml')
        settingStruct = parseXML('settings.xml');
        sessions.info.amplifier_sample_rate   = str2double(settingStruct.Attributes(1).Value); 
        sessions.info.lowpass_downsample      = str2double(settingStruct.Children(2).Attributes(85).Value); 
        sessions.info.Version                 = settingStruct.Attributes(3).Value;
        sessions.info.Name                    = settingStruct.Name;
        sessions.info.Children                = settingStruct.Children;
        
        % We calculate lowpass sampling rate
        sessions.info.lowpass_sample_rate     = sessions.info.amplifier_sample_rate / sessions.info.lowpass_downsample;
    else
        % Old datasets do not necessarily have an associatted .xml file
        % We have to look for basic info in the header. Some cannot be
        % populated yet.
        sessions.info.amplifier_sample_rate   = sessions.info.INTAN_hdr.frequency_parameters.amplifier_sample_rate  ;
        sessions.info.lowpass_downsample      = []; 
        sessions.info.Version                 = str2double(sessions.info.INTAN_hdr.version);
        sessions.info.Name                    = [];
        sessions.info.Children                = [];
        
        % Lowpass sampling rate not available
        sessions.info.lowpass_sample_rate     = [];
    end
    clear settingStruct

    %% 03. Proceed with conversion, depending on desired path
    
    % Format IS 'fileperch' or 'filepertype'
    if ~strcmp(sessions.info.fileformat,'Intanformat')
         
         % We want a .NWB file
            % Run wrapper for the INTAN to NWB functionality:
            %   A good thing is that it does not care about the file format, 
            %   the NWB functionality will put any version together.
         if input.useNWB(1)
            % This NEEDS A PYTHON installation and the tooldbox inside!
            % Detailed explanation:
            % WHAT IT IS: function to convert data from INTAN to .NWB format.
            % WHAT IT DOES: Checks for Python engine in computer. Adds the necessary
            %  dependences. Locates input session, copies ALL files to the IntanToNWB
            %  folder and merges them into a new 'info.nwb' file. This file 
            %  is renamed to 'session_name.nwb'. Moves this new file back to 
            %  the original session folder. Removes the copied data from the 
            %  IntanToNWB folder.
            %
            % Requires Python installed in the machine. 
            %  To date, MATLAB 2021b accepts up to Python 3.9. Install the
            %  64 bits version:
            % (https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
            %  To check access to Python Modules from MATLAB, look that 'pe' is correctly populated when running the script.
            intan2NWB_wrapper(sessions)

            % TODO: There seems to be an error when using the NWB to Matlab
            % functionalities, where after a first run, the consecutive
            % ones will not find an .mex file. Not sure why.
        
            % Then, Convert the NWB into pseudo-FieldTrip
            % INPUT:    Sessions info
            %           Current session ordinal
            %           Input about 'nwbmat' core functions being built.
            % OUTPUT:   Data. Is a Pseudo-FieldTrip structure
            %           Modified input, if it was first run.
            cd(input.mainfolder) % Back to code folder
            [data, input.useNWB(2)] = nwb2fieldtrip(sessions, input.useNWB(2));
        
            
         % We DO NOT want NWB.        
         else 
            % Run wrapper for the INTAN to MATLAB.
            %  Includes a mix of INTAN funtions. Can be run for 
            %  both 'filepertype' or 'fileperch'
            [data, sessions] = intan2MAT_wrapper(sessions);
            
         end

    % Format is NOT 'fileperch' or 'filepertype'     
    else
        disp('- Support for Intanformat files not yet available.')
        return
    end

    % Check for created FT file.
    if isempty(data) 
        % If true, there was an existing FT file, and the step was
        % skipped. Load the pre-existing one.
        disp('Loading...');
        datafile = dir('*_FT.mat');
        load(datafile.name, 'FT_data');
        clear datafile data

        % A couple details could have been lost if not processing
        sessions.info.nSamples       = FT_data.sampleinfo(2);
        sessions.info.recording_time = sessions.info.nSamples / (sessions.info.amplifier_sample_rate/32);
    end

    %% 04. Give proper FieldTrip format
    
    % If the file is being created now, look for 'data'
    if exist('data','var')

        % Check that Fieldtrip likes what we have (it should).
        FT_data = ft_checkdata(data, 'feedback' ,'yes');
        clear data
    
        % Then give the FT_data a proper 'continous' state.
        cfg = [];
         cfg.continuous = 'yes';
    
        FT_data = ft_redefinetrial(cfg, FT_data);
        clear cfg
    
        % Save this session data. Generates a file with continous data for a
        %   SINGLE session only into the session folder.
        save(strcat(sessions.folder,'\',sessions.name,'\',sessions.name,'_continous_FT.mat'), ...
            'sessions', 'input', 'FT_data', '-v7.3')
    end

    % If the file comes from a loaded file, it is named 'FT_data'
    % And it should be on real FT format already.

    % This test plotting uses Chronux Multitaper approach to generate fast
    % single-tappered Spectrograms on a subset of channels for a small chunck
    % of time. Just to have a preview of how the signal looks like
    if ~isempty(input.test_ch)
        plot_testsignal(FT_data,input.test_ch)
    end

    clear data2save cfg FT_data
end

%% Plot funtion
function plot_testsignal(FT_data,ch)

for i=ch
    data = FT_data.trial{1,1}(i,:)';

    % Parameters for MT
    params.tapers = [1 1 1]; % [W T p] -> 2TW-p tapers are used 
    params.Fs = FT_data.fsample;
    params.fpass = [0 100];
    params.pad = 0;
    movingwin = [1 .25];  % movingwin(1) MUST == T
       
    % COMPUTE MULTI-TAPER SPECTROGRAM AND Z-SCORE-IT
    [spectral.data, spectral.t, spectral.f] = mtspecgramc(data, movingwin, params);
    cmax = max(max(max(mag2db(spectral.data))));
    cmin = min(min(min(mag2db(spectral.data))));
    
    % PLOT IT
    figure(1), pcolor(spectral.t, spectral.f, mag2db(spectral.data)'); hold on;
        shading interp; 
        colormap("turbo");
        set(gca, 'clim', [cmin*0.1 cmax]); 
        c = colorbar('location','eastoutside');
        c.Label.String = 'dB';
        ylabel('Frequency (Hz) or Voltage (uV/50)', 'fontsize', 12);
        xlabel('sec', 'fontsize', 12);
        plot(FT_data.time{1}, ((data)/50)+35, 'Color', 'w', 'LineWidth', 1)
        xlim([15 55]); 
        ylim([0 50]);
    title(sprintf('Spectrogram for Ch: %s',FT_data.label{i,1}));
%     sdf(1,'800x300_hdpi_spectrum'), box('off');
    saveas(gca,sprintf('testsignal_ch%s',FT_data.label{i,1}),'png')
    close all
end
end