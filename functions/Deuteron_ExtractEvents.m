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
%
% Jesus. 21.12.2023

% Hardcoded variables (TO Reduce)
%IncludeEventFile = 1; %  We always include EVENT000.DF1 bc does not matter much.
%filePrefix       = 'NEUR'; % Need supressed by keeping it invariable at recording time
%minFileIndex     = 1; % ordinal of the first file to load (e.g. for NEUR0003, = 3); Need supressed by always having the SD card formatted.
%folderName       = opt.PathRaw; % bridged
maxFileIndex     = length(dir([opt.PathRaw '\NEUR*'])) - 1; % cero indexed, so [0:Nfiles-1]
count            = 1; % just a counter for processed files

%% Set up files to load 
% if IncludeEventFile % Always include
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

%% Final number of files to process. Create List
numberOfFiles = length(listOfFilesToLoad);

% executable likes a list of files as a character array: 'NEUR0001 NEUR0001 ... NEURNNNN'
if opt.useexe
    listOfFilesToLoadchar = [];
    for i=1:numberOfFiles
        listOfFilesToLoadchar = [listOfFilesToLoadchar ' ' cell2mat(listOfFilesToLoad(i))];
    end

else % dll likes a list of files as a .NET array
    fileNames = NET.createArray('System.String',numberOfFiles);
    for i = 1:numberOfFiles
         fileNames.Set(i - 1, listOfFilesToLoad{i});
    end
end

%% Load events
% For executable just command system('file.exe NEUR0001 NEUR0002 ... NEURNNNN recordfile.csv')
if opt.useexe 
    s = system([opt.exefile, ...                              % use full path to executable
                listOfFilesToLoadchar, ' ', ...               % use char vector of full list of files
                opt.FolderProcDataMat, '\EventRecord.CSV']);  % export to .cvs 

    myRecord = readmatrix(fullfile(opt.FolderProcDataMat, '\EventRecord.csv'), 'OutputType','string');
    numberOfRecords = length(myRecord);

else  % Much less intuitive to use with same final results. DEPRECATING
    % To cancel this while it is running, type c.Cancel()
    
    % Load in assembly
    asminfo = NET.addAssembly(opt.ReaderDll);        % loads in .NET dll 
    c = Event_File_Reader_9_0.EFRMatlabFunctions();  % an object containing the dll functions 
    c.Initialize();                                  % initialize the dll
    
    lh = addlistener(c, 'WriteFileLoaded', @(o, e) fprintf('Event File loaded: %d\n', e.number)); % adds listener to read progress
    
    errorCode = c.LoadFiles(fileNames); % load eventlog file, returns an int as error code
     
    % Check if LoadFile executed successfully
    if (errorCode ~= 0) % if LoadFile returned an error, display appropriate error message
        errorStr = c.GetErrorAsString(errorCode);
        error(char(errorStr));
        return;
    end
        
    % Get number of records.
    pause(300) % Give some time to the process in the background (used to be necessary without listener)
    numberOfRecords = c.GetNumberOfRecords(); % get number of records in event log
end

fprintf(['The number of records is: ' num2str(numberOfRecords) '\n']);

%% Loop through records and add them to EventRecords struct array.
if opt.useexe
    for recIdx = numberOfRecords:-1:1 % iterates backwards to preallocate array by assigning the final index first.
        EventRecord(recIdx).EventNumber = str2double(char(myRecord(recIdx,1))); 
        EventRecord(recIdx).TimeStamp = char(myRecord(recIdx,2));
        EventRecord(recIdx).TimeMsFromMidnight = str2double(char(myRecord(recIdx,3)));
        EventRecord(recIdx).TimeSource = char(myRecord(recIdx,4));
        EventRecord(recIdx).EventType = char(myRecord(recIdx,5));
        EventRecord(recIdx).Details = char(myRecord(recIdx,6));
    end

else % DEPRECATING
    for recIdx = numberOfRecords:-1:1 % iterates backwards to preallocate array by assigning the final index first.
        myRecord = c.GetIndexedRecord(recIdx - 1); 
        EventRecord(recIdx).EventNumber = str2double(char(myRecord(1))); 
        EventRecord(recIdx).EventType = char(myRecord(2));
        EventRecord(recIdx).TimeStamp = char(myRecord(3));
        EventRecord(recIdx).TimeMsFromMidnight = str2double(char(myRecord(4)));
        EventRecord(recIdx).TimeSource = char(myRecord(5));
        EventRecord(recIdx).Details = char(myRecord(6));
    end
end
fprintf('Successfully loaded file into EventRecords struct.\n');

%% Use eventlog to determine number of channels.
filestarted = find(strcmp({EventRecord.EventType}, 'File started')==1);
geninfo = split(EventRecord(filestarted(1)).Details, ";");
geninfo = regexp(geninfo,'\d*','Match');
opt.numChannels = str2double(geninfo{3});

if isempty(opt.channelOrder)
    opt.channelOrder = 1:1:opt.numChannels;
end

%% TODO save event record and DigIn events at session folder
save((opt.FolderProcDataMat + "\EventRecord.mat"),"EventRecord","-mat");

end