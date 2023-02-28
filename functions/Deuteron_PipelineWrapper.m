function Deuteron_PipelineWrapper(sessions, ss, varargin)
% Adaptation from the common pipeline for Deuteron. Wraps up the most common 
% processing lines necessary to get data from Deuteron raw files. This
% includes the Neural data and the motion sensors, so far. Could be
% expanded to extract audio as well.
%
% DEPENDENCIES
%   Deuteron_EventFileReaderDll: To extract Event Record from Deuteron Block format.
%   Deuteron2Kilosort: To compile recorded data in a single file per channel.
%                      Can also split the data based on event codes. 
%   Deuteron_GetMotionSensors: To extract data from motion sensors.
%   Deuteron_PlotMotionSensors: To process and visualize data from motion sensors.
%
% INPUTS:
%    sessions: struct. Variable containing info about sessions in process
%    ss:       int. Current session ordinal in the pipeline
%    in:        struct. optional inputs to override the defaults:
%                   RetrieveEvents: logic. possibility to load event codes to build a restriced matrix
%                                   (e.g., the matrix starts at the first 'itiON' and ends at 'end'
%                                    experiment, removing paradigm irrelevant periods) 
%                   StpSz: int. relative to HDF5file: chunks in which ...  
%
% GENERATES:
%    EventRecord.mat file and compressed events file.
%    .bin file, as channels x samples. If requested.
%    .h5 file, as channels x sample. If requested.
%       (both with channels in increasing order as required for Kilosort.)
%    MotionData.mat file, with [Accelerometer, Gyroscope, Magnetometer] variables
%       containing timeseries for each sensor readings, in physical units. Plus
%       a 'rotators' variable, containing the quaternions to create the
%       rotation matrices and other transformations.
%    Also, a file named 'D2K.mat' containing used inputs and outputs.
% 
% 24.02.2023 (Jesus)

if nargin < 3, opt = struct();
elseif nargin == 3, opt = varargin{1};
end

%% Defaults
if ~isfield(opt,'RetrieveEvents'),       opt.RetrieveEvents       = true;    end
if ~isfield(opt,'h5'),                   opt.h5                   = true;    end
if ~isfield(opt,'bin'),                  opt.bin                  = true;    end
if ~isfield(opt,'GetMotionSensors'),     opt.GetMotionSensors     = true;    end
if ~isfield(opt,'StpSz'),                opt.StpSz                = 1000000; end

% Paths and naming
if ~isfield(opt,'PathRaw'),              opt.PathRaw              = pwd;                                                  end
if ~isfield(opt,'FolderSingleChannels'), opt.FolderSingleChannels = fullfile(pwd,'oneFilePerChannel');                    end
if ~isfield(opt,'FolderProcDataMat'),    opt.FolderProcDataMat    = sessions.savefolder;                                  end
if ~isfield(opt,'SavFileName'),          opt.SavFileName          = sessions.list(ss).name;                               end

if ~isfield(opt,'DllFolder'),            opt.DllFolder            = 'C:\Code\Scripts\ephys-data-pipeline\functions\dlls'; end
if ~isfield(opt,'ReaderDll'),            opt.ReaderDll            = fullfile(opt.DllFolder, 'Event_File_Reader_8_3.dll'); end

%% Event data, using dll
if ~isfile('COMP_EVENTS.DF1')
    if opt.RetrieveEvents
        % Proceed to extract all events during session. Give some feedback.
        disp('Retrieving Events from Deuteron...')
        [EventRecord, sessions.info{ss}.numChannels] = ...
            Deuteron_EventFileReaderDll(opt, sessions, ss);
        disp(EventRecord);
        disp(['Found ', int2str(sessions.info{ss}.numChannels), ' channels']);

    else
       disp('Event extraction not requested, skipping...')
    end

else % Probably only useful while testing.
    disp('Found collected Events from Deuteron, skipping...')
    load("EventRecord.mat", "EventRecord");

    % Use event log to determine number of channels.
    modechange = find(strcmp({EventRecord.EventType}, 'Mode change')==1);
    geninfo = split(EventRecord(modechange(1)+1).Details, ";");
    geninfo = regexp(geninfo,'\d*','Match');
    sessions.info{ss}.numChannels = str2double(geninfo{3});
    clear geninfo modechange
end

%% Neural Data To .bin and .h5.
if ~isfile([sessions.list(ss).name '.h5'])
    if opt.h5 || opt.bin
        % Converts Deuteron DT2 and DF1 files into the 'oneFilePerChannel' format.
        % Creates full single files (.bin and .h5) to further use (i.e. with Kilosort)
        % Creates and saves D2K.mat file with few details (TODO, necessary?)
        disp('Generating single channel files from Deuteron...')
        Deuteron2Kilosort(opt, sessions, ss);
    else
        disp('No Neural data found in folder, but also not requested. Skipping...')
    end

else
    disp('Found .h5 file, skipping and recovering Header...')
    % TOD, case where data was already extracted. prob only useful while
    % testing.
end

%% Motion Data to Matlab
if opt.GetMotionSensors
    if ~isfile('MotionData.mat')
        disp('Extracting Motion Sensor data from Deuteron...')
        [Accelerometer, Gyroscope, Magnetometer] = ...
            Deuteron_GetMotionSensors(opt, sessions, ss);
    else
        disp('Motion Sensor data file found. Loading...')
        load("MotionData.mat","Accelerometer","Gyroscope","Magnetometer");
    end
    
    disp('Processing and Plotting Motion Sensor data.')
    % Add (..., 1, 1) to input, if visualization and video recording are wanted.
    Deuteron_PlotMotionSensors(Accelerometer, Gyroscope, Magnetometer, [], [])
end

end