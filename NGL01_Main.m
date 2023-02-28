%% Jesus' Pipeline to read INTAN continous data
% Will read and process INTAN, DEUTERON or ALLEGO data, from selected sessions for a given animal.
%
% The path to data will be: '...\datafolder\animal\dates'. This is 'rawfolder'
%     'datafolder', is any folder specified by user. 
%                   Different ones can be used for different projects.
%     'animal', is a folder with all sessions for an unique animal, 
%               identified with a 3 character code, i.e. 'DOE'.
%     'dates' are folders for single recordings named 'DOE_YYYYMMMDD'. 
%             Could have appends like '_01', '_Deut'... that need to be explicited.
% This pattern is constructed from the given inputs.
%
% DEPENDENCIES:
% Requires that all pipeline dependencies are properly set as matlab path. 
% 'set_default' will take care of this when a proper 'mainfolder' is provided.
%
% INPUTS:
%       mainfoldder:    chr array.  Full path to pipeline Code folder, as 'C:\...'. 
%                                   To include all dependencies.
%       datafolder:     string.     Full path to data folder, as "D:\...".
%                                   Where raw data is got from.
%       processed:      string.     A folder where newly created files will be saved.
%                                   Default is '\processed'
%       animal:         chr array.  A 3 character code as 'FAT', '427' ..., agreed upon.
%       dates:          cell of chr array. Specfic dates as {'yyyymmdd' 'yyyymmdd' ...} 
%                       or chr array. 'all'
%                                   At testing stages, 'yyyymmdd_system' or variations may exists
%                                   Default: 'all'
%       useNWB:         true/false  To create/skip NWB file.
%                                   Default: true.
%       ExtractData:    true/false  To create/skip binary and h5 files. Also extract motion sensor 
%                                   data if it comes from Deuteron.
%                                   Default: true.
%       plots:          int array.  If not empty, to draw plots, as [1 0 0]
%                                   for [spectrograms , raster, raster&traces].
%                                   Default: empty [].
%       test_ch:        int array.  If not empty, to plot snippet of requested channels.
%                                   As i.e. [1, 2, 5:15, 32].
%                                   Default: empty [].
% OUTPUTS
% For one single session or for a batch of sessions, from one single animal:
%       A set of default inputs and paths.
%       A list of sessions, with their associated info.
%       Fieldtrip (.mat), binary (.bin), HDF5 (.h5) and .nwb files from
%           1. Deuteron .DT2 or .DF1 data.
%           2. INTAN file-per-type and file-per-channel format data.
%           3. (ALLEGO data?)
%       EventRecord.mat file, from Deuteron session.
%       MotionData.mat file, From Deuteron sensors.
%       Plots snippets of time- and frequency-domain data.
%   .mat, .h5 and .bin files will be saved under ../rawfolder/processed
%       
% Last modified 27.02.2023 (Jesus)

% TODO LIST 
%       There seems to be an ERROR on 2nd and following runs of the NWB functionalities.
%       Figure out what's going on with the NWB/H5 DLLs that block either when the other has been performed...
%       'high' bandpass not yet available
%       create the wrapper for an INTAN fileperchannel format to Kilosort
%       Create Fieldtrip files from Deuteron data
%       Figure out how to work with Allego files (most likely, after Allego's self preprocessing tool?)
%       Save the 'sessions' variable by default, at the end, with a date timestamp perhaps?
%       Continue with 'Deuteron_GetDigInEvents' when I get a recording with EVENTS
%       Create a 'trial-parsed' stream in MAT2FieldTrip.
%

%% Inputs. 
% If not provided, a promp will ask for them or 'set_default' will use the 
% defaults. It will also put them in the correct format if an incorrect one 
% was given.
input.mainfolder = []; % 'C:\Code\Scripts\ephys-data-pipeline';
input.datafolder = []; % 'D:\Experiments\';
input.animal     = []; % '420';
input.processed  = []; % ie: 'processed' (default). A subfolder will be created inside the session folder

% Dates will be set to 'all' if missing here. 
% while testing, 'yyyymmdd_system' or other variations may exists
input.dates       = [{'20230217_01' '20230217_02' '20230220_Int' '20230221_Int' '20230222_Int'}]; % can be left empty, 'all', or a list like:
                                    % {'20230217_01' '20230217_02'...
                                    % '20230220_Deut' '20230221_Deut'...
                                    % '20230222_Deut' '20230220_Int'...
                                    % '20230221_Int' '20230222_Int'}; 

% Due to a conflict at h5 python-matlab dlls, when the two following pipelines 
% are requested, the NWB will perform well but the data extraction will not. 
% It will crash for not completely known reason. It needs a Matlab restart between runs.
input.ExtractData = true;  
input.useNWB      = false; 
    input.pyfolder = []; % Only needed if useNWB = true. Recommended 'C:\Code\Python39\IntanToNWB'

% These apply to FieldTrip-ready .mat files, only.
input.plots      = []; % An logic array of 0/1s, to ask for specific plots. See details.
input.test_ch    = []; % An array of numerals for channels to plot.

%% 00. Check inputs, set defaults and dependencies.
set_default(input);

%% 01. Find and list sessions. 
% Read requested sessions from specified animal folder.
% This 'sessions' variable can be used for summary, book keeping and
% debugging at the end of the pipeline. But it will not be saved
% automatically. % smt TODO?
sessions = findSessions(input);

%% 02. Loop sessions to process
for ss = 1:sessions.nSessions
    % Progress report
    txt = sprintf('\n --> Session %d out of %d. Session name: %s \n', ss, sessions.nSessions, sessions.list(ss).name);
    fprintf(txt);

    %% 03. Check file type and versions
    % Navigate to session raw dat folder.
    cd(fullfile(sessions.folder,sessions.list(ss).name));
    
    % Check System and version, based on existing files. Get info. 
    sessions.info{ss} = chckV();
    disp(sessions.info{ss});

    % Determine where processed data will be saved
    sessions.info{ss}.savefolder = fullfile(pwd, input.processed);
    txt = strcat('Process data will be saved to: >', sessions.info{ss}.savefolder);
    disp(txt);
    clear txt

    %% 04. Determine pipeline based on type of data
    switch sessions.info{ss}.fileformat
        case {'DT2', 'DT4', 'DT8', 'DAT', 'DF1'}
            %% 04.1 Deuteron Pipeline. Neural Data
               
            % 01 TODO Create Fieldtrip files from Deuteron data
            % So far, Deuteron does not seem ideal for LFP, but it should be possible at some point.
             %%%
             % Then, here will go the LFP extraction and conversion to NWB? and FT.
             %%%   

            % 02 Create .bin (and .h5) files with spiking data from highpass data
            if input.ExtractData
               opt = struct();
               % To modify optional inputs:
                opt.h5                = false;
                opt.bin               = true;
                opt.RetrieveEvents    = false;
                opt.GetMotionSensors  = false;

               % TODO: implement the new format conversion
               Deuteron_PipelineWrapper(sessions, ss, opt);
            end

        case {'fileperch', 'filepertype'}
          %% 04.2 INTAN Pipeline
          % 01. Find out INTAN settings and header file. Extract info.
          %  Uses a modified Intan function, to make the basic information
          %  available at 'sessions.info{ss}' and a more detailed info at
          %  the '.INTAN_hdr' sub-structure.
          sessions = findIntanSetting(sessions, ss);

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
              intan2NWB_wrapper(input, sessions, ss)
          end 

          % 03. Run wrapper for the INTAN to Kilosort. Creates .bin and .h5 files
          if input.ExtractData
              % Based on Sara, Aylin and Lukas' scripts.
              % To modify optional inputs, can be done here:
              opt = struct();
                opt.h5                = true;
                opt.bin               = true;

%               Intan2Kilosort_wrapper(sessions, ss, opt);
              Intan2Kilosort_wrapperV2(sessions, ss, opt);
          end

          % 04. Run wrapper for the INTAN to MATLAB.
          % Includes a mix of INTAN funtions. Outputs 'data' with plain
          % format. Can be feeded into next step for FT transformation.
          [data, sessions] = intan2MAT_wrapper(input, sessions, ss);

          % 05. CREATE and GIVE proper FieldTrip format. Give 'EventRecord'
          % variable as last input, if wanted to be trial-parsed. 
          % If the file comes from a loaded file, it will be named 'FT_data'
          % And it should be on real FT format already. otherwise, it
          % creates it.
          MAT2FieldTrip(input, data, sessions, ss, []);

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
    if ~isempty(input.test_ch) & any(input.plots)
        plot_testsignal(FT_data,input.test_ch)
    end
    
    %% 06 Clean up to move on to next session
    clear data2save cfg FT_data

end 
