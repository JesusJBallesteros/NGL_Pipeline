function [EventRecord, opt] = Deuteron_ExtractEvents(opt)
% Use the Deuteron's application without invoking the GUI.
% This example creates a struct called EventRecords that has a length of the 
% number of records in the event log file with the following fields:
% The dll requires a set of dark commands to perform the same as exe, much clearer
% INPUT: 
%       opt: struct with relevant info about paths and requirements.
% OUTPUT:
%       EventRecords: struct with all events recorded during session.
%           EventNumber (double)
%           EventType (string)  
%           TimeStamp (string)
%           TimeMsFromMidnight (double)
%           TimeSource (string)
%           Details (string)
%           TimeBreak (Nx2 cell array)
%
% 23.08.2024, Jesus: Consolidation and integration in pipeline, with two
%                    possible ways to get the session information.
% 24.04.2025, Jesus: Added function to detect time breaks when getting 
%                    information from Event logs. Can detect more than one.
%                    Only implemented for 'extractFromExe' function.
%                    Added output field .TimeBreak to match the detection
%                    from Deuteron System. Unlikely that they will happen on
%                    Intan Systems, so it will just be an empty 1x2 cell array.
    
%% Case
if opt.useexe,  [EventRecord, opt] = extractFromExe(opt);
else,           EventRecord = extractFromLog(opt);
end

end

%% Actual functions
function [EventRecord, opt] = extractFromExe(opt)
    % Use the Event_File_Reader_X_X or the .exe application without invoking the GUI
    % from a Deuteron recording with Block Format.
    % This example creates a struct called EventRecords that has a length of the 
    % number of records in the event log file with the following fields:
    % The dll requires a set of dark commands to perform the same as exe, much clearer
    % INPUT: 
    %       opt: struct with relevant info about paths and requirements.
    % OUTPUT:
    %       EventRecords: struct with all events recorded during session.
    %           EventNumber (double)
    %           EventType (string)  
    %           TimeStamp (string)
    %           TimeMsFromMidnight (double)
    %           TimeSource (string)
    %           Details (string)
    % Jesus. 23.08.2024

    %% Hardcoded variables
    maxFileIndex     = length(dir([opt.PathRaw '\NEUR*'])) - 1; % cero indexed, so [0:Nfiles-1]
    count            = 1; % just a counter for processed files
    
    %% Set up files to load 
    listOfFilesToLoad =  cell(maxFileIndex + 1, 1); % cell(maxFileIndex - 1 + 2, 1);
    
    for fileIdx = 1:maxFileIndex
        indexStr = num2str(fileIdx,'%04.f');
        if opt.useexe
            listOfFilesToLoad{count} = strcat('NEUR', indexStr, '.DF1');
        else
            listOfFilesToLoad{count} = fullfile(opt.PathRaw, strcat('NEUR', indexStr, '.DF1'));
        end
        count = count + 1;
    end
    
    % Create list of final files to process.
    numberOfFiles = length(listOfFilesToLoad);
    
    % Executable requires a list as char array: 'NEUR0001 NEUR0001 ... NEURNNNN'
    listOfFilesToLoadchar = [];
    for i=1:numberOfFiles
        listOfFilesToLoadchar = [listOfFilesToLoadchar ' ' cell2mat(listOfFilesToLoad(i))];
    end
    
    %% Load events
    % For executable just command system('file.exe, [char array of files], output.csv').
    % Input to system is actually a single one of class char array. The spaces in between 'subinputs' need to be explicited.
    % if opt.useexe % Preferred way to go, due to simplicity.
    s = system([opt.exefile, ...                              % use full path to executable
                listOfFilesToLoadchar, ' ', ...               % use char vector of full list of files
                opt.FolderProcDataMat, '\EventRecord.CSV']);  % export to .cvs 
        
    myRecord = readmatrix(fullfile(opt.FolderProcDataMat, '\EventRecord.csv'), 'OutputType', 'string'); % Read the output cvs
    
    % Find those logs with Digital-IN info.
    edgeDect = contains(myRecord(:,8), 'Digital in'); % Could be found in (:,6) as well. Same tstamp
    bitRecord = myRecord(edgeDect, 1:8); % Keep fields 1:8

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
    
    % By default, recordings starts as [1 1 0 0], but this is not recorded.
    words(1,:) = [1 1 0 0];

    % Get change direction from log description (8th column) (raising == 1, falling == 0)
    edgeDirection = contains(bitRecord(:,8), 'rising'); % categorize rising and falling edges.
    
    % Get changed Pin from the same description. (8th column)
    pin = regexp(bitRecord(:,8),'\d*','Match', 'once'); % Match the general expression '\d*', only once.
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
    mapDetc = find(contains(myRecord(:,8), 'Channel'), 1, "first"); % Find the log for a new file started. % Find the log for a new file started.
    geninfo = split(myRecord(mapDetc,8), "="); % Split the text contained in Details using semicolons.
    geninfo = regexp(geninfo,'\d*','Match'); % Match the general expression '\d*'.
    opt.channelOrder = str2double(geninfo{2});
    opt.numChannels = numel(opt.channelOrder); % Transform the 3rd field (hardcoded) into double.

    % Place extracted information into a proper EventRecord
    EventRecord.EventNumber         = double(1:1:length(EventType))';
    EventRecord.EventType           = single(EventType)';
    EventRecord.TimeStamp           = string(bitRecord(:,3)); % Convert to string array
    EventRecord.TimeMsFromMidnight  = str2double(bitRecord(:,4));
    EventRecord.TimeSource          = nan(length(bitRecord(:,6)),1);
    EventRecord.Details             = nan(length(bitRecord(:,8)),1);
    EventRecord.TimeBreak           = timebreak;

    if ~isempty(timebreak{2})
        %EventRecord.TimeBreak{1,2} = EventRecord.TimeBreak{1,2}-EventRecord.TimeMsFromMidnight(1);
        tbreakdur = EventRecord.TimeBreak{1,2}(2) - EventRecord.TimeBreak{1,2}(1);
        tbidx = EventRecord.TimeMsFromMidnight > EventRecord.TimeBreak{1,2}(1);
        EventRecord.TimeMsFromMidnight(tbidx) = EventRecord.TimeMsFromMidnight(tbidx) - tbreakdur;
        warning('The timebreak has been fixed and the EventRecpord will be saved for check.')
    end

    fprintf('Successfully created ''EventRecord'' structure.\n');
    
    % %% Save event record and DigIn events (TODO) at session folder
    % save((opt.FolderProcDataMat + "\EventRecord.mat"),"EventRecord","-mat");
end

function EventRecord = extractFromLog(opt)
    % Use this customized function to extract events from a text file
    % containing the log from a Deuteron recording with Block Format.
    % This example creates a struct called EventRecords that has a length of the 
    % number of records in the event log file with the OUTPUT fields.
    % INPUT: 
    %       opt: struct with relevant info about paths and options.
    % OUTPUT:
    %       EventRecords: struct with all events recorded during session.
    %           EventNumber (double)
    %           EventType (single)  
    %           TimeStamp (string)
    %           TimeMsFromMidnight (double)
    %           TimeSource (NaN)
    %           Details (NaN)
    % Jesus. 30.05.2024
    
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
    pinChange = cellfun(@cell2mat,pinChange,'UniformOutput', false); % Convert each cell to matrix
    pinChange = single(str2double(pinChange)); % Convert all values to single
    
    % Update also the valid timestamps
    tsmsec = tsmsec(idx==1); % Keep only those related to valid events
    ts = ts(idx==1); % Keep only those related to valid events
    
    % retrieve new status received by pin
    pinStatus = regexp(logevents(:,6),'\d','Match'); % Find expressions of input channel state and others
    pinStatus = pinStatus(idx==1); % Keep only those related to valid events
    pinStatus = cellfun(@cell2mat,pinStatus,'UniformOutput',false); % Convert each cell to matrix
    pinStatus = single(str2double(pinStatus)); % Convert all values to single
    
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
    EventRecord.TimeStamp           = string(ts); % Convert to string array
    EventRecord.TimeMsFromMidnight  = tsmsec;
    EventRecord.TimeSource          = nan(length(stateLog),1);
    EventRecord.Details             = nan(length(stateLog),1);
    EventRecord.TimeBreak           = {[] []};

end