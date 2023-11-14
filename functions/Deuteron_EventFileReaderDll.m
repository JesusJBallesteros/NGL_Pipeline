function [EventRecord, nchan] = Deuteron_EventFileReaderDll(opt)
% Use the Event_File_Reader_X_X dll to extract events from a Deuteron
% recording with Block Format.
% The user must enter the full path of the dll and the file they wish to load. 
% This example creates a struct called EventRecords that has a length of the 
% number of records in the event log file with the following fields:
% INPUT: 
%       full path to the dll file
% OUTPUT:
%       EventRecords struct with all events recorded during session with:
%           EventNumber (double)
%           EventType (string)  
%           TimeStamp (string)
%           TimeMsFromMidnight (double)
%           TimeSource (string)
%           Details (string)
%
% From Deuteron software, modified by Jesus. 01.02.2023

%% Default settings
if ~isfield(opt,'useexe'),             opt.useexe                 = 0;         end
filePrefix = 'NEUR';

% EVENT000.DF1 may be included among the files containing events but it is
% not required. You can assemble your cell array of file names in listOfFilesToLoad
% with your own custom code, but you must afterwards run the following segment ("change this to a .NET array") 
IncludeEventFile = 1;
folderName = opt.PathRaw;

%% Set up files to load 
% Block file format
minFileIndex = 1; % the number of the first file to load (e.g. for NEUR0003, set minFileIndex = 3);
maxFileIndex = length(dir([folderName '\' filePrefix '*'])) - 1; % cero indexed, so [0:Nfiles-1]
count = 1;
if IncludeEventFile
    listOfFilesToLoad = cell(maxFileIndex - minFileIndex + 2, 1);
    if opt.useexe
        listOfFilesToLoad{1} = 'EVENT000.DF1';
    else
        listOfFilesToLoad{1} = fullfile(folderName, 'EVENT000.DF1');
    end
    count = count + 1;
end

for fileIdx = minFileIndex:maxFileIndex
    indexStr = num2str(fileIdx,'%04.f');
    if opt.useexe
        listOfFilesToLoad{count} = strcat(filePrefix, indexStr, '.DF1');
    else
        listOfFilesToLoad{count} = fullfile(folderName, strcat(filePrefix, indexStr, '.DF1'));
    end
    count = count + 1;
end
numberOfFiles = length(listOfFilesToLoad);

if opt.useexe
    listOfFilesToLoadchar = [];
    for i=1:numberOfFiles
        listOfFilesToLoadchar = [listOfFilesToLoadchar ' ' cell2mat(listOfFilesToLoad(i))];
    end
else
    %% change this to a .NET array
    fileNames = NET.createArray('System.String',numberOfFiles);
    for i = 1:numberOfFiles
         fileNames.Set(i - 1, listOfFilesToLoad{i});
    end
end

%% Load events
if opt.useexe 
    s = system([opt.exefile, ...                                        % use full path to executable
                listOfFilesToLoadchar, ...                              % use char vector of full list of files
                ' ', opt.FolderProcDataMat, '\EventRecord.CSV']);       % export to .cvs 
else    
    % to cancel this while it is running, type c.Cancel()
    % Load in assembly
    asminfo = NET.addAssembly(opt.ReaderDll);        % loads in .NET dll 
    c = Event_File_Reader_8_3.EFRMatlabFunctions();
    c.Initialize();
    
    lh = addlistener(c, 'WriteFileLoaded', @(o, e) fprintf('Event File loaded: %d\n', e.number));
    
    errorCode = c.LoadFiles(fileNames); % load eventlog file, returns an int as error code
     
    % Check if LoadFile executed successfully
    if (errorCode ~= 0) % if LoadFile returned an error, display appropriate error message
        errorStr = c.GetErrorAsString(errorCode);
        error(char(errorStr));
        return;
    end
    
    % To compress files instead:
    % to cancel this while it is running, type c.Cancel()
    errorCode = c.CompressFiles(fileNames, [folderName '\COMP_EVENTS.DF1']);
    
    % Check if executed successfully
    pause(30)
    if (errorCode ~= 0) % if LoadFile returned an error, display appropriate error message
        errorStr = c.GetErrorAsString(errorCode);
        error(char(errorStr));
        return;
    end
    
    % get number of records
    % Offer output about number of records
    numberOfRecords = c.GetNumberOfRecords(); % get number of records in event log
    fprintf(['The number of records is: ' num2str(numberOfRecords) '\n']);
end

%% Loop through records and add them to EventRecords struct array.
if opt.useexe
    myRecord = readmatrix(fullfile(opt.FolderProcDataMat, '\EventRecord.csv'), 'OutputType','string');
    numberOfRecords = length(myRecord);
    for recIdx = numberOfRecords:-1:1 % iterates backwards to preallocate array by assigning the final index first.
        EventRecord(recIdx).EventNumber = str2double(char(myRecord(recIdx,1))); 
        EventRecord(recIdx).TimeStamp = char(myRecord(recIdx,2));
        EventRecord(recIdx).TimeMsFromMidnight = str2double(char(myRecord(recIdx,3)));
        EventRecord(recIdx).TimeSource = char(myRecord(recIdx,4));
        EventRecord(recIdx).EventType = char(myRecord(recIdx,5));
        EventRecord(recIdx).Details = char(myRecord(recIdx,6));
    end

else
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

%% Use event log to determine number of channels.
filestarted = find(strcmp({EventRecord.EventType}, 'File started')==1);
geninfo = split(EventRecord(filestarted(1)).Details, ";");
geninfo = regexp(geninfo,'\d*','Match');
nchan = str2double(geninfo{3});

%% Proceed to extract DigIn events from full event record
% TODO when I get a session with DIGIn events
%     Events = Deuteron_GetDigInEvents(EventRecord);

%% TODO save event record and DigIn events at session folder
save((opt.FolderProcDataMat + "\EventRecord.mat"),"EventRecord","-mat");

end