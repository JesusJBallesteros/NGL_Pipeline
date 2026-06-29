function [EventRecord, opt] = Deuteron_ExtractEvents(input, opt)
% Deuteron_ExtractEvents  Extract and parse the event log from a Deuteron session.
%
% PURPOSE:
%   Produces an EventRecord struct in the same format as INTAN_ExtractEvents
%   for downstream compatibility. Dispatches on opt.uselog:
%     false (default) — invokes Event_File_Reader_9_0.exe on NEUR*.DF1 files
%                       to generate an EventRecord.CSV, then parses it.
%                       Also sets opt.channelOrder and opt.numChannels from
%                       the channel-mapping entry in the CSV.
%     true            — parses a logevents.txt text log directly; used when
%                       events were not transmitted to the system but were logged.
%   A session-start sanity check trims any events preceding the first
%   EventType 0 or 8 (recording start).
%
% USAGE:
%   [EventRecord, opt] = Deuteron_ExtractEvents(input, opt)
%   Called from EventProcess when input format is 'DF1' or 'DT2'.
%
% INPUTS:
%   input  - struct; must contain:
%              .exefile   full path to Event_File_Reader_9_0.exe
%   opt    - options struct; relevant fields:
%              .uselog        false = use EXE (default); true = use logevents.txt
%              .PathRaw       raw data folder (NEUR*.DF1 files location)
%              .FolderProcDataMat  output folder for EventRecord.CSV
%
% OUTPUTS:
%   EventRecord  - struct with fields:
%                    .EventNumber          (double)
%                    .EventType            (double)
%                    .TimeStamp            (string)
%                    .TimeMsFromMidnight   (double), relative to first ts
%                    .TimeSource           (string)
%                    .Details              (string)
%                    .TimeBreak            (Nx2 cell) — populated by check_timebreaks
%   opt          - updated with (EXE path only):
%                    .channelOrder   active channel IDs from the channel-map log entry
%                    .numChannels    numel(opt.channelOrder)
%
% Last modified 29.06.2026 (Jesus)

%% Case
if opt.uselog, EventRecord = extractFromLog(opt);
else,          [EventRecord, opt] = extractFromExe(input, opt);
end

% Expected session start check
if ~ismember(EventRecord.EventType(1),[0,8])
    firstEv = find(EventRecord.EventType==8,1,"first");
    EventRecord.EventType(1:firstEv-1)           = [];
    EventRecord.EventNumber(1:firstEv-1)         = [];
    EventRecord.TimeStamp(1:firstEv-1)           = [];
    EventRecord.TimeMsFromMidnight(1:firstEv-1)  = [];
    EventRecord.TimeSource(1:firstEv-1)          = [];
    EventRecord.Details(1:firstEv-1)             = [];
    warning('While extracting events, first event was neither 0 nor 8. Every event before the first 8 was removed.\n')
    fprintf('Successfully created ''EventRecord'' structure.\n');
end


end

%% Actual functions
function [EventRecord, opt] = extractFromExe(input, opt)
    % Use the Event_File_Reader_X_X or the .exe application without invoking the GUI
    % from a Deuteron recording with Block Format.

    %% Hardcoded variables
    maxFileIndex     = length(dir([opt.PathRaw '\NEUR*']));
    count            = 1; % just a counter
    
    %% Set up files to load 
    listOfFilesToLoad =  cell(maxFileIndex, 1);
    
    for fileIdx = 0:maxFileIndex-1
        indexStr = num2str(fileIdx,'%04.f');
        % if opt.useexe
            listOfFilesToLoad{count} = strcat('NEUR', indexStr, '.DF1');
        % else
        %     listOfFilesToLoad{count} = fullfile(opt.PathRaw, strcat('NEUR', indexStr, '.DF1'));
        % end
        count = count + 1;
    end
    
    % Create list of final files to process.
    % numberOfFiles = length(listOfFilesToLoad);
    
    % Executable requires a list as char array: 'NEUR0001 NEUR0001 ... NEURNNNN'
    listOfFilesToLoadchar = [];
    for i=1:maxFileIndex
        listOfFilesToLoadchar = [listOfFilesToLoadchar ' ' cell2mat(listOfFilesToLoad(i))];
    end
    
    %% Load events
    % For executable just command system('file.exe, [char array of files], output.csv').
    % Input to system is actually a single one of class char array. The spaces in between 'subinputs' need to be explicited.
    % if opt.useexe % Preferred way to go, due to simplicity.
    s = system([input.exefile, ...                              % use full path to executable
                listOfFilesToLoadchar, ' ', ...               % use char vector of full list of files
                opt.FolderProcDataMat, '\EventRecord.CSV']);  % export to .cvs 
        
    myRecord = readmatrix(fullfile(opt.FolderProcDataMat, '\EventRecord.csv'), 'OutputType', 'string'); % Read the output cvs
    
    % Remove headers if existing
    hasHeader = strcmp(myRecord(1,1),'Event number');
    if hasHeader, myRecord(1,:) = []; end

    % Find those logs with Digital-IN info.
    digCol = size(myRecord,2);
    if digCol>6
        warning('This log contains 8 columns instead of the regular 6.')
        myRecord(:,7:8) = [];
        digCol = size(myRecord,2);
    end
    edgeDect = contains(myRecord(:,digCol), 'edge on pin'); 
    bitRecord = myRecord(edgeDect, :); % Keep all fields

    % Number of remainer records
    numberOfRecords = length(bitRecord); 
    fprintf(['Events extracted. The number of records is: ' num2str(numberOfRecords) '\n']);

    % Check for time breaks in the session
    [timebreak] = check_timebreaks(myRecord);
    if ~isempty(timebreak{2})
        warning('A time break has been found.');
    end

    %% Translate edge detections into binary words. 
    % Every detected edge is a change of pin to either 1 (rising) or 0 (falling).
    % Prepare a variable with all 4 pins, all set to zero
    words = zeros(numberOfRecords+1, 4);
    
    % By default, recordings start as [1 1 0 0], but this is not recorded.
    words(1,:) = [1 1 0 0];

    % Get change direction from log description (8th/6th column) (raising == 1, falling == 0)
    edgeDirection = contains(bitRecord(:,digCol), 'rising'); % categorize rising and falling edges.
    
    % Get changed Pin from the same description. (8th/6th column)
    pin = regexp(bitRecord(:,digCol),'\d*','Match', 'once'); % Match the general expression '\d*', only once.
    pin = single(str2double(pin)); % make it single array

    % Place corresponding rising changes into corresponding pins
    for i = 1:numberOfRecords
        words(i+1,:) = words(i,:); % get bits current status
        words(i+1,pin(i)) = edgeDirection(i); % set to 1 or 0 as coded
        
        % translate the resulting word to decimal
        EventType(i) = binvec2dec(flip(words(i+1,:))); 
    end

    % Use EventRecord to determine number of channels.
    % As a final account for active channels, we use the explicit log about it
    % that Deuteron provides with every new file created while recording.
    mapDetc = find(contains(myRecord(:,digCol), 'Channel'), 1, "first"); % Find the log for a new file started. % Find the log for a new file started.
    geninfo = split(myRecord(mapDetc,digCol), "="); % Split the text contained in Details using '='.
    geninfo = regexp(geninfo,'\d*','Match'); % Match the general expression '\d*'.
    opt.channelOrder = str2double(geninfo{2}); % The Ch numbers should be the second part
    opt.numChannels = numel(opt.channelOrder);

    % Place extracted information into a proper EventRecord
    if size(char(myRecord(2,2)),2)==1, add = 1; 
    else, add=0; end

    EventRecord.EventNumber         = double(1:1:length(EventType))';
    EventRecord.EventType           = single(EventType)';
    EventRecord.TimeStamp           = string(bitRecord(:,2+add)); % Convert to string array
    EventRecord.TimeMsFromMidnight  = str2double(bitRecord(:,3+add)) - str2double(myRecord(1,3+add)); % in ms from start recording
    EventRecord.TimeSource          = nan(size(bitRecord,1),1);
    EventRecord.Details             = bitRecord(:,digCol);
    EventRecord.TimeBreak           = timebreak;

    if ~isempty(timebreak{2})
        %EventRecord.TimeBreak{1,2} = EventRecord.TimeBreak{1,2}-EventRecord.TimeMsFromMidnight(1);
        tbreakdur = EventRecord.TimeBreak{1,2}(2) - EventRecord.TimeBreak{1,2}(1);
        tbidx = EventRecord.TimeMsFromMidnight > EventRecord.TimeBreak{1,2}(1);
        EventRecord.TimeMsFromMidnight(tbidx) = EventRecord.TimeMsFromMidnight(tbidx) - tbreakdur;
        warning('The timebreak has been fixed and the EventRecord will be saved for check.')
    end    
end

function EventRecord = extractFromLog(opt)
    % Use this customized function to extract events from a text file
    % containing the log from a Deuteron recording with Block Format.
    
    if ~isfield(opt,'delimiters'),  opt.delimiters  = {',','='};    end
    if ~isfield(opt,'outputas'),    opt.outputas    = 'string';     end
    if ~isfield(opt,'iniPins'),     opt.iniPins     = [1 1 0 0];    end
    
    %% Read the text file containing the Deuteron log and output a matrix using
    % the given delimiters. By default it should output a matrix where columns are:
    % [local time, msec after midnight, SpikeLog SN, local HH:MM:SS.MSEC, InputCh, InputState, Port]
    logfile = "logevents.txt";
    logevents = readmatrix(logfile, 'OutputType', opt.outputas, 'Delimiter', opt.delimiters);
    
    %% Retrive all msec after midnight (column 2)
    tsmsec = str2double(logevents(:,2));
    
    % Retrieve all timestamps (column 4). They come as HH:MM:SS.mmmmmm
    ts = regexp(logevents(:,4),'(\d+):(\d+):(\d+).(\d+)','Match');
    
    % Retrieve pin number receiving status change (column 5)
    pinChange = regexp(logevents(:,5),'\d','Match'); % Find matching expressions to a single digit
    idx = cell2mat(cellfun(@length,pinChange,'UniformOutput', false)); % assess size of results
    pinChange = pinChange(idx==1); % Keep only those of length=1
    % pinChange = cellfun(@cell2mat,pinChange,'UniformOutput', false); % Convert each cell to matrix
    pinChange = cellfun(@str2num,pinChange,'UniformOutput', false); % Convert each cell to matrix
    pinChange = single(cell2mat(pinChange)); % Convert all values to single
    
    % Update also the valid timestamps
    tsmsec = tsmsec(idx==1); % Keep only those related to valid events
    ts = ts(idx==1); % Keep only those related to valid events
    
    % retrieve new status received by pin
    pinStatus = regexp(logevents(:,6),'\d','Match'); % Find expressions of input channel state and others
    pinStatus = pinStatus(idx==1); % Keep only those related to valid events
    % pinStatus = cellfun(@cell2mat,pinStatus,'UniformOutput',false); % Convert each cell to matrix
    pinStatus = cellfun(@str2num,pinStatus,'UniformOutput',false); % Convert each cell to matrix
    pinStatus = single(cell2mat(pinStatus)); % Convert all values to single
    
    %% Create a log of all pin states (including the initial one) and a vector
    % with the decimal values of such states
    stateLog = [opt.iniPins; zeros(size(pinChange,1), size(opt.iniPins,2))];
    newState = stateLog(1,:);
    for i=2:length(stateLog)
        newState(pinChange(i-1)) = pinStatus(i-1);
        stateLog(i,:) = newState;
    end
    stateLog(1,:) = []; % remove initial state, 
    % added artificially (as it IS the exisiting initial pinState but it IS NOT sent by the paradigm 
    % in the current session as part of it, but set by Deuteron as default when the system
    % boots up. Also we send it at the end of any previous session, to replicate this fact.
    
    stateLog = int2str(stateLog);
    
    %% Place extracted information into a proper EventRecord
    EventRecord.EventNumber         = double(1:1:length(stateLog))';
    EventRecord.EventType           = single(bin2dec(stateLog));
    EventRecord.TimeStamp           = string(ts); % Convert to string array. This keeps the real time
    EventRecord.TimeMsFromMidnight  = tsmsec - tsmsec(1); % Relativize this one to the beginning
    EventRecord.TimeSource          = nan(length(stateLog),1);
    EventRecord.Details             = nan(length(stateLog),1);
    EventRecord.TimeBreak           = {[] []};
end                                                                                                                                      