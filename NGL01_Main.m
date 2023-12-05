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

% Last updated JESUS 31.10.2023

%% USER Inputs. Check A, B and C.
% A) CRITICAL
% Specify drive and folder where data is located AND this toolbox folder (If not already added to MATLAB folder system)
input.datadrive     = 'D:\';
input.studyName     = 'Pilot_SocialLearning'; 
input.toolbox       = 'C:\Code\ephys-data-pipeline'; % Default: 'C:\Code\ephys-data-pipeline'

% B) SUBJECTS AND SESSIONS
% To run the script on all subjects and sessions included in your project,
% use char array 'all'. For a session-to-session process, explicit the subject 
% and session/s to process using cell arrays.
input.subjects       = {'485'}; % char array 'all', or a cell with a single subject denomination e.g. {'DOE'} or {'042'}
input.dates          = {'20231110'}; %'all'; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}

% C) GENERAL Options. 
% Those used for all sessions. Specific options can be set below or defaulted in the functions.
opt = struct(); % leave this, to empty possible residues from a previous run.

    opt.RetrieveEvents      = true;  % Retrieve event log from Deuteron system.
        opt.useexe          = false; % Eventually, only option for Deuteron recordings (TODO)
        opt.usepar          = true;  % temporarily use of .par files from Juan's behavior paradigm
    
    opt.parsetrial          = true;      % Define trials based on retrieved events

    opt.bin                 = false;   % Creation of .bin file, for Kilosort.

    opt.kilosort            = false;   % Call to kilosort processing. 
        % NEEDS configfile saved under '...\analysisCode'
        opt.spkTh           = -2; % def: -4.5. It will override the KS configfile.
        opt.KSchanMapFile   = []; %'chanMapPoly3Deut.mat'; % 'chanMapPoly3Deut' 'chanMapPoly3' 'chanMapATLASTri'

    opt.FTfile              = true;   % Creation of .mat file, FieldTrip ready.

    opt.GetMotionSensors    = false;  % JACOB gone MIA. Retrieve data from motion sensors in Deuteron.

%% 00. Check current inputs.
% Will set the rest of default inputs and dependencies.
input = set_default(input);

% Loop subjects.
for s = 1:input.nsubjects
    %% 01. Find and list sessions. Determine the pipeline.
    % Read requested sessions from specified animal folder.
    sessions = findSessions(input);

    % Loop sessions.
    for ss = 1:sessions(s).nsessions
        % Navigate to session's raw data folder.
        cd(fullfile(sessions(s).folder,sessions(s).list{ss}));
                
        % Report.
        txt = sprintf('\n --> Subject %s, session %d out of %d: %s \n', ...
                     input.subjects(s).name, ss, sessions(s).nsessions, sessions(s).list{ss});
        fprintf(txt);
    
        % Check system and version.
        sessions(s).info = [];
        sessions(s).info = chckV();
    
        % Determine where processed session data will be saved.
        opt.PathRaw           = pwd;
        opt.FolderProcDataMat = fullfile(input.processed, input.subjects(s).name, sessions(s).list{ss});
        opt.behavFiles        = fullfile(input.bhvfolder, input.subjects(s).name, sessions(s).list{ss});
        opt.SavFileName       = sessions(s).list{ss}; 
        
        % Report and create folder.
        disp(strcat('Processed data will be saved to: >', opt.FolderProcDataMat));
        mkdir(opt.FolderProcDataMat);
    
        % Determine pipeline based on type of data.
        switch sessions(s).info.fileformat
            case {'DT2', 'DF1'} 
                %% 02.1 Deuteron Pipeline
                if input.ExtractData
                   disp('Deuteron data is NOT being filter, by default');
                   sessions(s) = Deuteron_PipelineWrapper(sessions(s), input, opt);
                end
    
            case {'fileperch', 'filepertype'}
                %% 02.2 INTAN Pipeline
                if input.ExtractData
                   sessions(s) = INTAN_PipelineWrapper(sessions(s), input, opt);
                end
    
            case 'Allego'
                warning('Allego format not implemented yet.');
        
            case 'NAN'
                warning('The format of this session could not be recognized. Skipping.');
                continue
    
            otherwise
                warning('Something went wrong during format verification. Skipping');
                sessions(s).info.fileformat = 'ERR'; % Flag for ERROR
                continue
        end 
                
        %% 03 Kilosort
        if opt.kilosort
            % Kilosort will run without GUI.
            master_kilosort(sessions(s), input, opt) % 'opt' is this pipeline running variable.
        end
    
        % Clean up to move on to next session
        clear FT_data INTANdata txt

    end % sessions loop
end % subjects loop