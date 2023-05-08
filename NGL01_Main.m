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
%    Prepare a small function that, if the IKN folder structure does not exists, it will be created and all data included there.
%    There seems to be an ERROR on 2nd and following runs of the NWB functionalities.
%    Figure out what's going on with the NWB/H5 DLLs that block either when the other has been performed...
%    Prepare to downsample highpass data to a half? For Data size reduction.
%    Continue with 'Deuteron_GetDigInEvents' when I get a recording with EVENTS
%    Create a 'trial-parsed' stream in 'mat2FieldTrip' VS. add post-hoc parsing
%    Figure out how to work with Allego files (most likely, after Allego's self preprocessing tool?)
%

%% Input storage drive, project name and toolbox folder:
input.datadrive     = 'D:\';
input.studyName     = 'ephysTest'; % For SPP people: 'Dorian\SPP'
input.toolbox       = 'C:\Code\Scripts\ephys-data-pipeline'; % Default: 'C:\Code\Scripts\ephys-data-pipeline'

% To run the script on all subjects and sessions included in your project,
% just leave both as 'all'. For a session-to-session process, explicit the
% subject and session/s to process. 
input.subjects       = {'478'};    % char array 'all', or a single subject denomination e.g. 'DOE'
input.dates          = {'20230505'};    % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}

%% General Options. What you want to obtain:
% Those used for all sessions. Specific options can be set below or defaulted in the functions.
opt = struct();
    % These are essential.
    opt.bin               = true;   % Creation of .bin file, for Kilosort.
    opt.kilosort          = true;   % Call to kilosort processing. NEEDS configfile and chanmap in \analysisCode 
    opt.FTfile            = false;   % Creation of .mat file, FieldTrip ready.
    opt.RetrieveEvents    = true;   % Retrieve event log from Deuteron system.
    opt.GetMotionSensors  = false;  % JACOB gone MIA. Retrieve data from motion sensors in Deuteron.

    % KS options. Leave commented to default, or explicit your own paths here (NOT Recommended). 
        % opt.KSConfigFile   = 'D:\ephysTest\analysisCode\'; % Location of configfile
        % opt.KSchanMapFile  = 'D:\ephysTest\analysisCode\chanMap*.mat'; % Location AND name of channel map file

    % Only for FieldTrip .mat files. Plots snippets of raw signals and spectrograms.
        % opt.test_ch    = []; % An array of numerals for channels to plot.
        
%% 00. Check inputs, set defaults and dependencies.
set_default(input);

%% 01. Find and list sessions, per animal
for s = 1:input.nsubjects
    % Read requested sessions from specified animal folder.
    sessions = findSessions(input);

    %% 02. Loop subjects and sessions to process.
    for ss = 1:sessions(s).nsessions
        % Navigate to session's raw data folder.
        cd(fullfile(sessions(s).folder,sessions(s).list{ss}));
                
        % Progress report.
        txt = sprintf('\n --> Subject %s, session %d out of %d: %s \n', ...
                     input.subjects(s).name, ss, sessions(s).nsessions, sessions(s).list{ss});
        fprintf(txt);
    
        %% 03. Check file type, version and folders.
        % Check System and version for current session.
        sessions(s).info = [];
        sessions(s).info = chckV();
    
        % Determine where processed data will be saved, done for every session.
        opt.PathRaw           = pwd;
        opt.FolderProcDataMat = fullfile(input.processed, input.subjects(s).name, sessions(s).list{ss});
        opt.SavFileName       = sessions(s).list{ss}; 
        
        % Report and create folder.
        disp(strcat('Processed data will be saved to: >', opt.FolderProcDataMat));
        mkdir(opt.FolderProcDataMat);
    
        %% 04. Determine pipeline based on type of data.
        switch sessions(s).info.fileformat
            case {'DT2', 'DF1'} % 'DT4', 'DT8', 'DAT', never seen.
                %% 04.1 Deuteron Pipeline. Neural Data
                if input.ExtractData
                   % So far, we are NOT applying any filters, bc we are only
                   % recording high pass data.
                   disp('Deuteron data is NOT being filter, by default');
                   Deuteron_PipelineWrapper(sessions(s), opt);
                end
    
            case {'fileperch', 'filepertype'}
              %% 04.2 INTAN Pipeline
              % 01. Find out INTAN settings and header file. Extract info.
              %  Uses a modified Intan function, to make the basic information
              %  available at 'info{ss}' and a more detailed info at
              %  the '.INTAN_hdr' sub-structure.
              sessions(s) = findSetting(sessions(s));
    
              % 02. Create NWB file
              if input.useNWB % We want a .NWB file.
    
                  % Run wrapper for the INTAN to NWB functionality:               
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
                  intan2NWB_wrapper(input, opt);
              end 
    
              % 03. Run wrapper for the INTAN to Kilosort. Creates .bin and .h5 files
              if input.ExtractData 
                  if opt.bin
    
                  % Based on Sara, Aylin and Lukas' scripts.
                  Intan2Kilosort_wrapper(sessions(s), opt);
                  end
              end
              
              if opt.FTfile
                  % 04. Run wrapper for the INTAN to FIELDTRIP.
                  % Includes a mix of INTAN funtions. CREATES and GIVES proper
                  % FieldTrip format without trial-parsing. 
                  intan2FieldTrip(sessions(s), opt)

                  % 04.1 Plotting. Uses Chronux Multitaper approach to generate fast
                  % single-tappered Spectrograms on a subset of channels for a small chunck
                  % of time. Just to have a preview of how the signal looks like in
                  % the LFP range.
                  if isfield(input, 'test_ch') && ~isempty(input.test_ch)
                      plot_testsignal(FT_data, input.test_ch, opt)
                  end
              end
    
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
                
        %% 04 Go Kilosorting!
        if opt.kilosort
            
            % Kilosort Run without GUI.
            master_kilosort(sessions, input) % ops (not opt) is created as kilosort settings
        
            % Manual curation of data in Phy would be done once all requested sessions are finished. Therefore, I think that
            % NGLXX_postPhy will be a new script, to run after all manual cuartion is done, allowing to recover the final
            % data from all sessions at once.
            
            % In principle, some basic curation could be automated (tossing very low firing rates, clusters with low
            % numbers of spikes, etc.)
        end
    
        % Clean up to move on to next session
        clear FT_data INTANdata txt

    end % sessions loop
end % subjects loop