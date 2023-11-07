%% Pipeline process INTAN and Deuteron continous data.
% Will read and process INTAN, DEUTERON (or ALLEGO) data, from selected sessions for a given animal.
% The main pipeline will be: INTAN/DEUTERON raw formats to be located, then
% converted to .bin files (spike sorting), and Fieldtrip .mat structures
% (for LFP). Once sorted, spike data will be attached to the FieldTrip
% structure. For Arena experiments, motion sensor data will be extracted and
% interpreted. Data will be trial-parsed using EventCodes.
%
% The hard disk data structure SHOULD fit the IKN standard published at:
% gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure
%
% DEPENDENCIES:
% Requires that all pipeline dependencies are properly located. 
% I suggest to include the 'mainfolder' in Matlab's permanent path system.
% The function 'set_default' will take care of the rest of folders on each run.
%
% INPUTS:
%       input.datadrive, char array with the drive where data is located. As 'D:\'
%       input.studyName, char array with the project name, matching the
%                          folder name where all data will be stored. As 'studyName'
%       input.subjects,  char array with either 'all' OR a single subject name e.g. 'DOE'
%       input.dates,     char array with either 'all' OR a cell array of dates 
%                          for a SINGLE subject e.g. {'YYYYMMDD' 'yyyymmdd' ...)
%
% OPTIONS: is a struct with many possible fields. All should have a
% corresponding default inside whatever function is being called. Main ones
% are:     
%     opt.bin,              Creation of .bin file, input to Kilosort 2/4.
%     opt.FTfile,           Creation of .mat file with FieldTrip format.
%     opt.RetrieveEvents,   Retrieve event log from Deuteron system.
%     opt.GetMotionSensors, Retrieve data from motion sensors in Deuteron.
%     opt.kilosort,         Asks to proceed with KS processing and waits to retrieve its results.
%     opt.set_filter,       If Deuteron data was adquired with a wideband.
%     opt.lowpass,          Lowpass band to extract LFP from wideband.
%     opt.highpass,         Highpass band to extract spike activity.
%
% OUTPUTS:
% For one single session or for a batch of sessions, from one single animal:
%       Fieldtrip (.mat), binary (.bin), HDF5 (.h5) and/or .nwb files from
%           1. Deuteron .DT2 or .DF1 data.
%           2. INTAN file-per-type and file-per-channel format data.
%           3. (ALLEGO data?)
%       EventRecord.mat file, from Deuteron session.
%       MotionData.mat file, From Deuteron sensors.
%       Plots snippets of time- and frequency-domain data, from FieldTrip
%       
% Last modified 05.04.2023 (Jesus)

% TODO LIST
% If Deuteron2Kilosort(opt) filter for DF1 format works, set filter out of format cases (generalize)
% Continue with 'Deuteron_GetDigInEvents' when we get a recording with EVENTS
% Check for Deuteron_GetDigInEvents(EventRecord) status.
% Check for FT trial-parsing using EventRecord with MAT2FieldTrip(data, opt, varargin)
%    Create a 'trial-parsed' stream in 'mat2FieldTrip' VS. add post-hoc parsing
% Extract nChannels from EventsRecord. Find first 'File started' then use
%    'strsplit(EventRecord(50).Details,{';','='})' and find the 6th cell
% There seems to be an ERROR on 2nd and following runs of the NWB functionalities.
%    Figure out what's going on with the NWB/H5 DLLs that block either when the other has been performed...

% Last updated 
% JESUS 31.10.2023

%% USER Inputs. Check A, B and C.
% A) CRITICAL
% Specify drive and folder where data is located AND this toolbox folder (If not already added to MATLAB folder system)
input.datadrive     = 'F:\';
input.studyName     = 'Pilot_SocialLearning'; 
input.toolbox       = 'C:\Code\ephys-data-pipeline'; % Default: 'C:\Code\ephys-data-pipeline'

% B) SUBJECTS AND SESSIONS
% To run the script on all subjects and sessions included in your project,
% use char array 'all'. For a session-to-session process, explicit the subject 
% and session/s to process using a cell array like {'001', ... , 'ETC'}.
input.subjects       = {'485' '257'}; % char array 'all', or a cell with a single subject denomination e.g. {'DOE'} or {'042'}
input.dates          = {'20231106'}; %'all'; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}

% C) GENERAL Options. 
% Those used for all sessions. Specific options can be set below or defaulted in the functions.
opt = struct(); % leave this, to empty possible residues from a previous run.

    opt.RetrieveEvents      = true;  % Retrieve event log from Deuteron system.
    opt.bin                 = true;   % Creation of .bin file, for Kilosort.

    opt.kilosort            = true;   % Call to kilosort processing. 
        % !! NEEDS configfile and chanmap saved under '...\analysisCode'
        opt.spkTh           = -2; % def: -4.5. It will override the KS configfile.
        opt.KSchanMapFile   = []; %'chanMapPoly3Deut.mat'; % 'chanMapPoly3Deut' 'chanMapPoly3' 'chanMapATLASTri'

    opt.FTfile              = true;   % Creation of .mat file, FieldTrip ready.
        % Only for FieldTrip .mat files. 
        opt.test_ch         = []; % Plots snippets of raw signals and spectrograms. An array of numerals for channels to plot.

    opt.GetMotionSensors    = false;  % JACOB gone MIA. Retrieve data from motion sensors in Deuteron.

%% 00. Check current inputs.
% Will set the rest of default inputs and dependencies.
input = set_default(input);

%% Run subjects
for s = 1:input.nsubjects
    %% 01. Find and list sessions, per animal
    % Read requested sessions from specified animal folder.
    sessions = findSessions(input);

    %% Loop subjects and sessions to process.
    for ss = 1:sessions(s).nsessions
        % Navigate to session's raw data folder.
        cd(fullfile(sessions(s).folder,sessions(s).list{ss}));
                
        % Progress report.
        txt = sprintf('\n --> Subject %s, session %d out of %d: %s \n', ...
                     input.subjects(s).name, ss, sessions(s).nsessions, sessions(s).list{ss});
        fprintf(txt);
    
        %% 02. Check Session type, version and folders.
        % Check System and version for current session.
        sessions(s).info = [];
        sessions(s).info = chckV();
    
        % Determine where processed data will be saved, for every session.
        opt.PathRaw           = pwd;
        opt.FolderProcDataMat = fullfile(input.processed, input.subjects(s).name, sessions(s).list{ss});
        opt.SavFileName       = sessions(s).list{ss}; 
        
        % Report and create folder.
        disp(strcat('Processed data will be saved to: >', opt.FolderProcDataMat));
        mkdir(opt.FolderProcDataMat);
    
        % Determine pipeline based on type of data.
        switch sessions(s).info.fileformat
            case {'DT2', 'DF1'} 
                %% 04.1 Deuteron Pipeline
                if input.ExtractData
                   disp('Deuteron data is NOT being filter, by default');
                   sessions(s) = Deuteron_PipelineWrapper(sessions(s), input, opt);
                end
    
            case {'fileperch', 'filepertype'}
                %% 04.2 INTAN Pipeline
                if input.ExtractData
                   sessions(s) = INTAN_PipelineWrapper(sessions(s), input, opt);
                end

            % Wrapped inside the above function, for cleaniness in this script
            %               % 01. Find out INTAN settings and header file. Extract info.
            %               %  Uses a modified Intan function, to make the basic information
            %               %  available at 'info{ss}' and a more detailed info at
            %               %  the '.INTAN_hdr' sub-structure.
            %               sessions(s) = findSetting(sessions(s));
            %     
            %               % 02. Create NWB file
            %               if input.useNWB % We want a .NWB file.
            %     
            %                   % Run wrapper for the INTAN to NWB functionality:               
            %                     % This NEEDS A PYTHON installation and the tooldbox inside!
            %                     % Detailed explanation:
            %                     % WHAT IT IS: function to convert data from INTAN to .NWB format.
            %                     % WHAT IT DOES: Checks for Python engine in computer. Adds the necessary
            %                     %  dependences. Locates input session, copies ALL files to the IntanToNWB
            %                     %  folder and merges them into a new 'info.nwb' file. This file 
            %                     %  is renamed to 'session_name.nwb'. Moves this new file back to 
            %                     %  the original session folder. Removes the copied data from the 
            %                     %  IntanToNWB folder.
            %                     %
            %                     % Requires Python installed in the machine. 
            %                     %  To date, MATLAB 2021b accepts up to Python 3.9. Install the
            %                     %  64 bits version:
            %                     % (https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
            %                     %  To check access to Python Modules from MATLAB, look that 'pe' is correctly populated when running the script.
            %                   intan2NWB_wrapper(input, opt);
            %               end 
            %     
            %               % 03. Run wrapper for the INTAN to Kilosort. Creates .bin and .h5 files
            %               if input.ExtractData 
            %                   if opt.bin
            %     
            %                     % Based on Sara, Aylin and Lukas' scripts.
            %                     % only if the .bin file does not exist yet.
            %                       if ~isfile(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']))
            %                         Intan2Kilosort_wrapper(sessions(s), opt);
            %                       end
            %                   end
            %               end
            %               
            %               if opt.FTfile
            %                   % 04. Run wrapper for the INTAN to FIELDTRIP.
            %                   % Includes a mix of INTAN funtions. CREATES and GIVES proper
            %                   % FieldTrip format without trial-parsing. 
            %                   intan2FieldTrip(sessions(s), opt)
            % 
            %                   % 04.1 Plotting. Uses Chronux Multitaper approach to generate fast
            %                   % single-tappered Spectrograms on a subset of channels for a small chunck
            %                   % of time. Just to have a preview of how the signal looks like in
            %                   % the LFP range.
            %                   if isfield(input, 'test_ch') && ~isempty(input.test_ch)
            %                       plot_testsignal(FT_data, input.test_ch, opt)
            %                   end
            %               end
    
            case 'Allego'
                warning('Allego format not implemented yet.');
    
            case 'Intanformat'
                warning('INTAN old format not implemented. Probably it wont be.');
    
            case 'NAN'
                warning('The format of this session could not be recognized. Skipping.');
                continue
    
            otherwise
                warning('Something went wrong during format verification. Skipping');
                sessions(s).info.fileformat = 'ERR'; % Flag for ERROR
                continue
        end 
                
        %% 04 Go Kilosorting
        if opt.kilosort
            % Kilosort Run without GUI.
            master_kilosort(sessions(s), input, opt) % 'opt' is a pipeline running variable.
        end
    
        % Clean up to move on to next session
        clear FT_data INTANdata txt
    end % sessions loop
end % subjects loop