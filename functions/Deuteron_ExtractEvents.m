function [EventRecord, opt] = Deuteron_ExtractEvents(opt)
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

% Jesus. 04.01.2024

%% Hardcoded variables (TO Reduce)
%IncludeEventFile = 1; %  We always include EVENT000.DF1 bc does not matter.
%filePrefix       = 'NEUR'; % Need supressed by keeping it invariable at recording time
%minFileIndex     = 1; % ordinal of the first file to load (e.g. for NEUR0003, = 3); Need supressed by always having the SD card emptied.
%folderName       = opt.PathRaw; % bridged

maxFileIndex     = length(dir([opt.PathRaw '\NEUR*'])) - 1; % cero indexed, so [0:Nfiles-1]
count            = 1; % just a counter for processed files

%% Set up files to load 
% if IncludeEventFile % Always included
listOfFilesToLoad =  cell(maxFileIndex + 1, 1); % cell(maxFileIndex - 1 + 2, 1);
if opt.useexe,     listOfFilesToLoad{count} = 'EVENT000.DF1';
else,              listOfFilesToLoad{count} = fullfile(opt.PathRaw, 'EVENT000.DF1');
end
count = count + 1;
% end

for fileIdx = 1:maxFileIndex
    indexStr = num2str(fileIdx,'%04.f');
    if opt.useexe
        listOfFilesToLoad{count} = strcat('NEUR', indexStr, '.DF1');
    else
        listOfFilesToLoad{count} = fullfile(opt.PathRaw, strcat('NEUR', indexStr, '.DF1'));
    end
    count = count + 1;
end

%% Create list of final files to process.
numberOfFiles = length(listOfFilesToLoad);

% Executable requires a list as char array: 'NEUR0001 NEUR0001 ... NEURNNNN'
if opt.useexe
    listOfFilesToLoadchar = [];
    for i=1:numberOfFiles
        listOfFilesToLoadchar = [listOfFilesToLoadchar ' ' cell2mat(listOfFilesToLoad(i))];
    end

else % dll requires the use of .NET array, whatever that is...
    fileNames = NET.createArray('System.String',numberOfFiles);
    for i = 1:numberOfFiles
         fileNames.Set(i - 1, listOfFilesToLoad{i});
    end
end

%% Load events
% For executable just command system('file.exe, [char array of files], output.csv').
% Input to system is actually a single one of class char array. The spaces in between 'subinputs' need to be explicited.

if opt.useexe % Preferred way to go, due to simplicity.
    s = system([opt.exefile, ...                              % use full path to executable
                listOfFilesToLoadchar, ' ', ...               % use char vector of full list of files
                opt.FolderProcDataMat, '\EventRecord.CSV']);  % export to .cvs 

    myRecord = readmatrix(fullfile(opt.FolderProcDataMat, '\EventRecord.csv'), 'OutputType','string'); % Read the output cvs
    numberOfRecords = length(myRecord); % Extract number of records

% else  % Much less intuitive, with same results. DEPRECATING
%     % To cancel this while it is running, type c.Cancel()
%     
%     % Load in assembly
%     asminfo = NET.addAssembly(opt.ReaderDll);        % loads in .NET dll 
%     c = Event_File_Reader_9_0.EFRMatlabFunctions();  % an object containing the dll functions 
%     c.Initialize();                                  % initialize the dll
%     
%     lh = addlistener(c, 'WriteFileLoaded', @(o, e) fprintf('Event File loaded: %d\n', e.number)); % adds listener to read progress
%     
%     errorCode = c.LoadFiles(fileNames); % load eventlog file, returns an int as error code
%      
%     % Check if LoadFile executed successfully
%     if (errorCode ~= 0) % if LoadFile returned an error, display appropriate error message
%         errorStr = c.GetErrorAsString(errorCode);
%         error(char(errorStr));
%         return;
%     end
%         
%     % Get number of records.
%     pause(300) % Give some time to the process in the background (used to be necessary without listener)
%     numberOfRecords = c.GetNumberOfRecords(); % get number of records in event log
end

fprintf(['Events extracted. The number of records is: ' num2str(numberOfRecords) '\n']);

%% Loop through records and add them to an EventRecord structure.
if opt.useexe
    % Iterates backwards, preallocating array by assigning the final index first.
    for recIdx = numberOfRecords:-1:1 
        EventRecord(recIdx).EventNumber = str2double(char(myRecord(recIdx,1)));         % ?
        EventRecord(recIdx).TimeStamp = char(myRecord(recIdx,2));                       % Time stamp ?
        EventRecord(recIdx).TimeMsFromMidnight = str2double(char(myRecord(recIdx,3)));  % milisecs from midnight
        EventRecord(recIdx).TimeSource = char(myRecord(recIdx,4));                      % Source of time stamp.
        EventRecord(recIdx).EventType = char(myRecord(recIdx,5));                       % ?
        EventRecord(recIdx).Details = char(myRecord(recIdx,6));                         % Extra details
    end

% else % DEPRECATING
%     for recIdx = numberOfRecords:-1:1 % iterates backwards to preallocate array by assigning the final index first.
%         myRecord = c.GetIndexedRecord(recIdx - 1); 
%         EventRecord(recIdx).EventNumber = str2double(char(myRecord(1))); 
%         EventRecord(recIdx).EventType = char(myRecord(2));
%         EventRecord(recIdx).TimeStamp = char(myRecord(3));
%         EventRecord(recIdx).TimeMsFromMidnight = str2double(char(myRecord(4)));
%         EventRecord(recIdx).TimeSource = char(myRecord(5));
%         EventRecord(recIdx).Details = char(myRecord(6));
%     end
end
fprintf('Successfully created ''EventRecord'' structure.\n');

%% Use EventRecord to determine number of channels.
% As a final account for active channels, we use the explicit log about it
% that Deuteron provides with every new file created while recording.
filestarted = find(strcmp({EventRecord.EventType}, 'File started')==1); % Find the log for a new file started.
geninfo = split(EventRecord(filestarted(1)).Details, ";"); % Split the text contained in Details using semicolons.
geninfo = regexp(geninfo,'\d*','Match'); % Match the general expression '\d*'.
opt.numChannels = str2double(geninfo{3}); % Transform the 3rd field (hardcoded) into double.
if isempty(opt.channelOrder)
    opt.channelOrder = 1:1:opt.numChannels; % Order channels as incremental ordinals. (TODO: this? perhaps match the Deuteron map?)
end

%% Save event record and DigIn events (TODO) at session folder
save((opt.FolderProcDataMat + "\EventRecord.mat"),"EventRecord","-mat");

end