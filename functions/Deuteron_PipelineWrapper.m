function Deuteron_PipelineWrapper(sessions, ss, varargin)
% Adaptation from the common pipeline for Deuteron. Prepares recorded data for spike sorting with Kilosort.
% For now, uses *.DT2 files and creates .h5 and .bin files.
%
% DEPENDENCIES
%    Deuteron2Kilosort: function to compile recorded data in a single file per channel.
%                      Can also filter the data and retrieve event codes. 
%    Deuteron_EventFileReaderDll: 
%
% INPUTS:
%    sessions: struct. Variable containing info about sessions in process
%    ss:       int. Current session ordinal in the pipeline
%    input:    struct. optional inputs to override the defaults:
%               createAvrgDatMat: logic. if true, another matrix (and respective binary file) are created 
%                                     with the average of all channels subtracted from every channel 
%               retrieveEvents: logic. possibility to load event codes to build a restriced matrix
%                                   (e.g., the matrix starts at the first 'itiON' and ends at 'end'
%                                    experiment, removing paradigm irrelevant periods) 
%               ApplyHighPassFilter: logic. Use High-pass filter
%               stpSz: int. relative to HDF5file: chunks in which ...  
%
% OUTPUT:
%    Binary file, channels(rows) per sample (columns), with channels
%    in increasing order as required for Kilosort.
%    h5 file, channels(rows) per sample (columns), with channels
%    in increasing order as required for Kilosort.
%    Also, a file named 'D2K.mat' containing used inputs and outputs.
% 
% 14.02.2023 (Jesus)

if nargin < 3, in = struct();
elseif nargin == 3, in = varargin{1};
end

%% Defaults
if ~isfield(in,'retrieveEvents'),       in.retrieveEvents       = true;     end
if ~isfield(in,'keeph5'),               in.keeph5               = true;     end
if ~isfield(in,'keepbin'),              in.keepbin              = true;     end
if ~isfield(in,'GetMotionSensors'),     in.GetMotionSensors     = true;     end

if ~isfield(in,'ApplyHighPassFilter'),  in.ApplyHighPassFilter  = false;    end
if ~isfield(in,'createAvrgDatMat'),     in.createAvrgDatMat     = false;    end
if ~isfield(in,'stpSz'),                in.stpSz                = 1000000;  end

% Paths and naming
if ~isfield(in,'pathRaw'),              in.pathRaw                  = pwd;                                                  end
if ~isfield(in,'folderSingleChannels'), in.folderSingleChannels     = fullfile(pwd,'oneFilePerChannel');                    end
if ~isfield(in,'folderProcDataMat'),    in.folderProcDataMat        = in.pathRaw;                                           end
if ~isfield(in,'dllFolder'),            in.dllFolder                = 'C:\Code\Scripts\ephys-data-pipeline\functions\dlls'; end
if ~isfield(in,'ReaderDll'),            in.ReaderDll                = fullfile(in.dllFolder, 'Event_File_Reader_8_3.dll');  end

if ~isfield(in,'savFileName'),          in.savFileName              = sessions.list(ss).name;                               end
if ~isfield(in,'savFileNameAvrg'),      in.savFileNameAvrg          = ['Averaged_',sessions.list(ss).name];                 end
if ~isfield(in,'folderProcDataMatAveraged'), in.folderProcDataMatAveraged = in.pathRaw;                                     end

%% Event data, using dll
if ~isfile('COMP_EVENTS.DF1')
    if in.retrieveEvents
        % Proceed to extract all events during session
        disp('Retrieving Events from Deuteron...')
        EventRecord = Deuteron_EventFileReaderDll(in, sessions, ss);
        
        % Use event log to determine number of channels.
        modechange = find(strcmp({EventRecord.EventType}, 'Mode change')==1);
        geninfo = split(EventRecord(modechange(1)+1).Details, ";");
        geninfo = regexp(geninfo,'\d*','Match');
        sessions.info{ss}.numChannels = str2double(geninfo{3});
        clear geninfo modechange

        % Proceed to extract DigIn events from full event record
        % TODO when I get a session with DIGIn events
    %     Events = Deuteron_GetDigInEvents(EventRecord);
        
        % TODO save event record and DigIn events at session folder
        save("EventRecord.mat","EventRecord","-mat");
    else
       disp('Event extraction not requested, skipping...')
    end
else
    disp('Found collected Events from Deuteronm, skipping...')
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
    % Converts Deuteron DT2 and DF1 files into the 'oneFilePerChannel' format.
    % Creates full single files (.bin and .h5) to further use (i.e. with Kilosort)
    % Creates and saves D2K.mat file with few details (TODO, necessary?)
    disp('Generating single channel files from Deuteron...')
    Deuteron2Kilosort(in, sessions, ss);
    
else
    disp('Found .h5 file, skipping and recovering Header...')
end

%% Motion Data to Matlab
if in.GetMotionSensors
    disp('Generating single channel files from Deuteron...')
    [Accelerometer, Gyroscope, Magnetometer, param] = Deuteron_GetMotionSensors(in, sessions, ss);

    % Plots 
    figure,
    if ~isempty(Accelerometer)
        MSData = Accelerometer;

        subplot(3,1,1)
        title('Accelerometer');
        % General plot
        plot(MSData.t, MSData.X); hold on
        plot(MSData.t, MSData.Y); hold on
        plot(MSData.t, MSData.Z);
        ylim([-param.acclMax param.acclMax]); ylabel('m/s^2');
        xlabel('ms');
        legend({'X' 'Y' 'Z'}, 'Box','off');
         box("off")
    end

    if ~isempty(Gyroscope)
        MSData = Gyroscope;

        subplot(3,1,2)
        title('Gyroscope');
        % General plot
        plot(MSData.t, MSData.X); hold on
        plot(MSData.t, MSData.Y); hold on
        plot(MSData.t, MSData.Z);
        ylim([-param.gyroMax param.gyroMax]); ylabel('deg/s')
        xlabel('ms');
        legend({'X' 'Y' 'Z'}, 'Box','off');
        box("off")
    end

    if ~isempty(Magnetometer)
        MSData = Magnetometer;

        subplot(3,1,3)
        title('Magnetometer');
        % General plot
        plot(MSData.t, MSData.X); hold on
        plot(MSData.t, MSData.Y); hold on
        plot(MSData.t, MSData.Z);
%         ylim([-param.magMax param.magMax]);
        ylim([-1e-4 1e-4]);  ylabel('Tesla')
        xlabel('ms');
        legend({'X' 'Y' 'Z'}, 'Box','off');
        box("off")
    end

end



end