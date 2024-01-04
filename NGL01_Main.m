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

% Version 03.01.2024 (Jesus)

% %% USER Inputs. Check A, B and C.
% % A) CRITICAL
% % Specify drive and folder where data is located AND this toolbox folder (If not already added to MATLAB folder system)
% input.datadrive     = 'F:\';
% input.studyName     = 'Pilot_SocialLearning'; 
% input.toolbox       = 'C:\Code\ephys-data-pipeline'; % Default: 'C:\Code\ephys-data-pipeline'
% 
% % B) SUBJECTS AND SESSIONS
% % To run the script on all subjects and sessions included in your project,
% % use char array 'all'. For a session-to-session process, explicit the subject 
% % and session/s to process using cell arrays.
% input.subjects       = 'all'; %{'485'}; % char array 'all', or a cell with a single subject denomination e.g. {'DOE'} or {'042'}
% input.dates          = 'all'; %{'20231113'}; %'all'; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}
% 
% % C) GENERAL Options. 
% % Those used for all sessions. Specific options can be set below or defaulted in the functions.
% opt = struct(); % leave this, to empty possible residues from a previous run.
%     opt.kilosort            = true; % Call to kilosort processing.                      !! NEEDS configfile saved under 'studyName\analysisCode\'
%         opt.spkTh           = -2.5;     % Default: -4.5.
%         opt.KSchanMapFile   = '';       % e.g.'chanMapPoly3Deut', 'chanMapATLASTri'
%     opt.bombcell            = true;     % Run bombcell on the KS output, previously to manual curation
%         opt.rerun           = true;     % To overwrite previous results or not
%         opt.nRawSpikesToExtract = 1000; % Parameter for bombcell run
%     opt.phy                 = true;     % Calls phy for manual inspection or curation. !! It PUTS MATLAB on HOLD!
%     opt.FieldTrip           = true;     % Creation of .mat file, FieldTrip ready.
%     opt.RetrieveEvents      = true;     % Retrieve event log.
%         opt.useexe          = false;    % Eventually, only option for Deuteron recordings (TODO)
%         opt.usepar          = true;     % temporarily use of .par files from Juan's behavior paradigm
%     opt.parsetrial          = false;    % Define and split data into trials
%     opt.GetMotionSensors    = false;    % Retrieve data from motion sensors in Deuteron. (TODO)

%% 00. Check current inputs.
% Will set the rest of default inputs and dependencies.
cd(input.toolbox)
input = set_default(input);

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

% Loop subjects.
for x = 1:input.nsubjects
    % Loop sessions.
    for y = 1:input.sessions(x).nsessions
        input.run = [x y]; % Store current run as input to pass to functions

        % Navigate to session's raw data folder and report.
        % Check system and version. Determine where processed session data will be saved.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);
    
        % Determine pipeline based on type of data.
        switch input.sessions(input.run(1)).info.fileformat
            case {'DT2', 'DF1'} 
               %% 02.1 Deuteron Pipeline
               disp('Deuteron data is NOT being filter, by default');
               Deuteron_PipelineWrapper(input, opt);
    
            case {'fileperch', 'filepertype'}
               %% 02.2 INTAN Pipeline
               % input.sessions(input.run(1)) = INTAN_PipelineWrapper(sessions(input.run(1)), input, opt); %mod
               INTAN_PipelineWrapper(input, opt);
    
            otherwise
               warning('Something went wrong during format verification. Skipping');
               continue
        end 
                
        %% 03 Kilosort
        if opt.kilosort
            % Kilosort will run without GUI.
            master_kilosort(input, opt) %mod
        end
    
        %% 04 Bombcell
        if opt.bombcell
            % Kilosort will run without GUI.
            Bombcell_Main(opt) 
        end

        %% 05 Open Phy to manual curation or just inspection
        if opt.phy
            % Will change to current session directory and open phy.
            % ! Keeps MATLAB busy until interface is closed.
            cd(opt.FolderProcDataMat)
            system('phy template-gui params.py');
        end

        % Clean up to move on to next session
        clear FT_data INTANdata txt

    end % sessions loop
end % subjects loop