%% Pipeline to start processing INTAN and Deuteron continous data.
% Will read and process INTAN, DEUTERON (or ALLEGO) data, from selected sessions for a given animal.
% Then, in general, will create single files in different formats for
% further processing. This includes merging many-files data into a
% single-continuous file, with filtering and downsampling as needed.
% Eventually will allow for creation of trial-parsed datafiles as well.
%
% The path to data will be: '...\datafolder\animal\dates'. This is 'rawfolder'
%     'datafolder', is any folder specified by user. 
%                   Different ones can be used for different projects.
%     'animal', is a folder with all sessions for an unique animal, 
%               identified with a 3 character code, i.e. '420' or 'FAT'.
%     'dates' are folders for single recordings named 'DOE_YYYYMMMDD'. 
%             Could have appends like '_01', '_Deut'... that need to be explicited.
% This pattern is constructed from the given inputs.
%
% DEPENDENCIES:
% Requires that all pipeline dependencies are properly set as matlab path. 
% I suggest to include the 'mainfolder' in Matlab's permanent path system.
% The function 'set_default' will take care of the rest of folders on each run.
%
% INPUTS:
%       mainfoldder:    chr array.  Full path to the code folder, as 'C:\...'. 
%                                   Where the toolbox lives.
%       datafolder:     string.     Full path to data folder, as "D:\...".
%                                   Where raw data will be searched for.
%       processed:      string.     A folder where newly created files will be saved.
%                                   Default would be '...\rawfolder\processed'
%       animal:         chr array.  A 3 character code as 'FAT', '427' ..., agreed upon.
%       dates:          cell of chr array. Specfic dates as {'yyyymmdd' 'yyyymmdd' ...} 
%                       or chr array. 'all'
%                                   'yyyymmdd_system' or other variations may exists.
%                                   Default: 'all'
%       ExtractData:    true/false  To create/skip binary and h5 files. Also extract motion sensor 
%                                   data if it comes from Deuteron.
%                                   Default: true.
%       useNWB:         true/false  To create/skip NWB file. Default: true.
%
% OUTPUTS:
% For one single session or for a batch of sessions, from one single animal:
%       Variables with inputs and paths.
%       Variable with a list of sessions and their associated info.
%       Fieldtrip (.mat), binary (.bin), HDF5 (.h5) and/or .nwb files from
%           1. Deuteron .DT2 or .DF1 data.
%           2. INTAN file-per-type and file-per-channel format data.
%           3. (ALLEGO data?)
%       EventRecord.mat file, from Deuteron session.
%       MotionData.mat file, From Deuteron sensors.
%       Plots snippets of time- and frequency-domain data, from FieldTrip
%   .mat, .h5 and .bin files will be saved under ../rawfolder/processed
%       
% Last modified 08.03.2023 (Jesus)

% TODO LIST 
%    There seems to be an ERROR on 2nd and following runs of the NWB functionalities.
%    Figure out what's going on with the NWB/H5 DLLs that block either when the other has been performed...
%    Prepare to downsample highpass data to a half? For Data size reduction.
%    Continue with 'Deuteron_GetDigInEvents' when I get a recording with EVENTS
%    Create a 'trial-parsed' stream in 'mat2FieldTrip'.
%    Consider allowing for more than one animal to be processed.
%    Figure out how to work with Allego files (most likely, after Allego's self preprocessing tool?)
%

%% Inputs. 
% The only input that cannot be defaulted is the ANIMAL to be used, for now:
input.animal     = '451'; % i.e '420' or FAT; 

% For the following, for any not provided a promp will pop-up. 
% If still empty, 'set_default' will use the defaults.
input.mainfolder = 'C:\Code\Scripts\ephys-data-pipeline'; % Default: 'C:\Code\Scripts\ephys-data-pipeline'
input.datafolder = 'D:\Experiments\';                     % Default: 'D:\Experiments\'
input.processed  = 'processed';                           % Default: 'processed'

% Make sure of your own file denominations. 
% Can be left empty, can be 'all', or can be a cell array like:
    % {'20230217_01' '20230217_02' '20230220_Deut' '20230221_Deut'...
    % '20230222_Deut' '20230220_Int' '20230221_Int' '20230222_Int'}; 
input.dates      = {'20230308_Int' '20230308_Deut'};   % Default: 'all'

% Due to a conflict at h5 python-matlab dlls, when the two following pipelines 
% are requested, the NWB will perform well but the data extraction will not. 
% It will crash for not completely known reason. It needs a Matlab restart between runs.
input.ExtractData = true;
input.useNWB      = false; % Meaning, do not run both 'true' (for now). 

%% Options.
% Here can go those used for all sessions. Variable ones (each session's path, or other) 
% can be set later on, normally automatized.
opt = struct();
    opt.h5                = false; % Creation of .h5 file (not really used, so far).
    opt.bin               = true;  % Creation of .bin file.
    opt.FTfile            = true;  % Create a FieldTrip-formatted .mat file.

    opt.RetrieveEvents    = false; % Retrieve event log from Deuteron system.
    opt.GetMotionSensors  = false; % Retrieve data from motion sensors in Deuteron.

%     opt.set_filter        = 0;     % Set to 1 when data is known to come as wideband 
                                   % i.e from a Deuteron recording with an open wideband.
                                   % In INTAN, this will be evaluated automatically, but with Deuteron is not, yet.
    opt.lowpass           = [  0  300]; % Lowpass band applied for FieldTrip pipeline.
    opt.highpass          = [300 7500]; % Highpass band applied for Kilosort pipeline.
    
    % This applies only to FieldTrip .mat files. Not really useful here other 
    % than for testing, or for checking that everything is running in a new 
    % dataset, to check for empty channels or other weird stuff. 
    % Can be used to plot snippets as example as well.
    opt.test_ch    = []; % An array of numerals for channels to plot.

%% 00. Check inputs, set defaults and dependencies.
set_default(input);

%% 01. Find and list sessions. 
% Read requested sessions from specified animal folder.
% This 'sessions' variable can be used for summary, book keeping and
% debugging at the end of the pipeline. But it will not be saved
% automatically. % smt TODO?
sessions = findSessions(input);

%% 02. Loop sessions to process.
for ss = 1:sessions.nSessions
    % Navigate to session raw data folder.
    cd(fullfile(sessions.folder,sessions.list(ss).name));
    
    % Progress report.
    txt = sprintf('\n --> Session %d out of %d: %s \n', ss, sessions.nSessions, sessions.list(ss).name);
    fprintf(txt);

    %% 03. Check file type, version and folders.
    % Check System and version, based on existing files. Get info. 
    sessions.info{ss} = chckV();

    % Determine where processed data will be saved, done for every session.
    opt.PathRaw           = pwd;
    opt.FolderProcDataMat = fullfile(pwd, input.processed);
    opt.SavFileName       = sessions.list(ss).name; 
    
    % Report and create folder.
    disp(strcat('Processed data will be saved to: >', opt.FolderProcDataMat));
    mkdir(opt.FolderProcDataMat);

    %% 04. Determine pipeline based on type of data.
    switch sessions.info{ss}.fileformat
        case {'DT2', 'DF1'} % 'DT4', 'DT8', 'DAT', never seen.
            %% 04.1 Deuteron Pipeline. Neural Data
            if input.ExtractData
               % So far, we are NOT applying any filters, bc we are only
               % recording high pass data.
               disp('Deuteron data is NOT being filter, by default');
               Deuteron_PipelineWrapper(sessions, ss, opt);
            end

        case {'fileperch', 'filepertype'}
          %% 04.2 INTAN Pipeline
          % 01. Find out INTAN settings and header file. Extract info.
          %  Uses a modified Intan function, to make the basic information
          %  available at 'sessions.info{ss}' and a more detailed info at
          %  the '.INTAN_hdr' sub-structure.
          sessions = findSetting(sessions, ss);

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
              if opt.h5 || opt.bin

              % Based on Sara, Aylin and Lukas' scripts.
              Intan2Kilosort_wrapperV2(sessions, ss, opt);
              end
          end
          
          if opt.FTfile
              % 04. Run wrapper for the INTAN to MATLAB.
              % Includes a mix of INTAN funtions. Outputs 'data' with plain
              % format. Can be feeded into next step for FT transformation.
              INTANdata = intan2mat_wrapper(sessions, ss, opt);
    
              % 05. CREATE and GIVE proper FieldTrip format. Give 'EventRecord'
              % variable as last input, if wanted to be trial-parsed. 
              % If the file comes from a loaded file, it will be named 'FT_data'
              % And it should be on real FT format already. otherwise, it
              % creates it.
              mat2FieldTrip(INTANdata, opt);
          end

        case 'Allego'
          %% 04.3 Allego Pipeline
            warning('Allego format not implemented yet.'); % TODO

        case 'Intanformat'
            warning('INTAN old format not implemented. Probably will not be.');

        case 'NAN'
            warning('The format of this session could not be recognized. Skipping.');
            continue

        otherwise
            warning('Something went wrong during format verification. Skipping');
            sessions.info{ss}.fileformat = 'ERR'; % Flag for ERROR
            continue
    end 

    %% 05 This test plotting uses Chronux Multitaper approach to generate fast
    % single-tappered Spectrograms on a subset of channels for a small chunck
    % of time. Just to have a preview of how the signal looks like in
    % the LFP range.
    if ~isempty(input.test_ch)
        plot_testsignal(FT_data, input.test_ch, opt)
    end
    
    %% 06 Clean up to move on to next session
    clear FT_data INTANdata txt

end 
