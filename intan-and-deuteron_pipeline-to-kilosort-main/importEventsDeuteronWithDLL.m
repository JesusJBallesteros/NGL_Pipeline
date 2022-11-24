function [events] = importEventsDeuteronWithDLL(dllFolder,eventLogPath,sampleRate)

% INPUT:
%       path of 'Event_File_Reader_7_2.dll' (provided in the Deuteron's
%       website)
%       path of the EVENTLOG.NLE file for the recording session
%       sampling rate

% OUTPUT:
%       cell array containing time points (as minutes from midnight), in the
%       first column, the sample number,
%       in the the second column, and the pins that
%       were turned on (Deuteron detects only 'rising' edges) in the third column
% 
%       if no events were sent to the system, the cell array returns the a
%       2x2 matrix with the start and end of the recording,as minutes from 
%       midnight, in the first column, and the sample numbers of these timestamps 
%       in the second column

% VERSION HISTORY:
% Author:         Aylin Apostel
% Version:        1.1.0
% Last Change:    12.11.2021
%
% 12.11.2021, Aylin: first draft 
% 31.03.2022, Sara: Added conversion to sample number and option to get 
%                   event output when no transitions were logged in the Deuteron system
% 19.04.2022, Sara: Added the interaction between Matlab and Deuteron's .NET
%                   app for event reading, in order to automate this process.
%                   No further need to convert the event file "EVENTLOG.NLE" into a .CSV. 
%                   Currently, this code only works for recordings in flat format
%%
samplesMs = sampleRate/1000;

%% Read vent file into Matlab:

% User settings
% dllFolder = 'C:\Users\ACN\Desktop\DeuteronFileViewer\newestVesin'; %path of dll folder
eventFileReaderDll = fullfile(dllFolder, 'Event_File_Reader_7_2.dll'); 
% eventLogPath = 'C:\Users\ACN\Desktop\Sara\Deuteron_NeuralEventRecords';

% Set up files to load for flat file format

numberOfFiles = 1;
listOfFilesToLoad = cell(numberOfFiles, 1);
listOfFilesToLoad{numberOfFiles} = fullfile(eventLogPath,'EventLog.NLE'); % event log file name with full path
fileNames = NET.createArray('System.String',numberOfFiles);
for i = 1:numberOfFiles
    fileNames.Set(i - 1, listOfFilesToLoad{i});
end

%
% Load in assembly - as in Deuteron example code 
asminfo = NET.addAssembly(eventFileReaderDll); % loads in .NET dll 
c = Event_File_Reader_7_2.EFRMatlabFunctions();
c.Initialize();
lh = addlistener(c, 'WriteFileLoaded', @(o, e) fprintf('File loaded: %d\n', e.number));
Task = c.LoadFiles(fileNames);
c.Cancel(); % if you choose to cancel while it is running
errorCode = Task.Result;

% Check if LoadFile executed successfully
if (errorCode == 0) % if no error
    fprintf('EventLog file was loaded successfully.\n')
else % if LoadFile returned an error, display appropriate error message
    errorStr = c.GetErrorAsString(errorCode);
    error(char(errorStr));
    return;
end

% Compile event cell
numberOfRecords = c.GetNumberOfRecords(); % get number of records in event log
EventRecords = cell(numberOfRecords,6);
for recIdx = numberOfRecords:-1:1 % iterates backwards to preallocate array by assigning the final index first.
    myRecord = c.GetIndexedRecord(recIdx - 1);
    EventRecords{recIdx,1} = str2double(char(myRecord(1)));
    EventRecords{recIdx,2} = char(myRecord(3));
    EventRecords{recIdx,3} = str2double(strrep(char(myRecord(4)), ',', '.'));
    EventRecords{recIdx,4} = char(myRecord(5));
    EventRecords{recIdx,5} = char(myRecord(2));
    EventRecords{recIdx,6} = char(myRecord(6));
end
fprintf('\nEvent data successfully extracted\n');

%% take timestamp of first file opening
fileStart = EventRecords(contains(EventRecords(:,5),'File started'),:);
firstFileStartTime = cell2mat(fileStart(1,3));

% Get all lines with 'Digital in':
digIn     = contains(EventRecords(:,5),'Digital in');

if ~any(digIn) %no event codes were sent to the logger 
    events = cell(2,2); %first column: ms from midnight; second column: number of sample 
    events{1,1} = firstFileStartTime; %opening of 1st file (ms from midnight)
    events{2,1} = cell2mat(EventRecords(end,3));%end of recording session(ms from midnight)
    if events{2,1} ==0 %a restart of the logger (internal workings of Deuteron)
        events{2,1} = cell2mat(EventRecords(end-1,3));
    end 
    events{1,2} = 1; %sample nr at file opening 
    events{2,2} = (events{2,1} - firstFileStartTime(1,1))*samplesMs; %last sample recorded 
    
else % events were sent 
    digitalIn = EventRecords(any(digIn,2),:);
    events = cell(size(digitalIn,1),3); %###
    events(:,1) = digitalIn(:,3); % all timepoints
    
    % Which pin changed?
    for i = 1 : size(digitalIn,1)
        dat = digitalIn(i,:);
        if contains(dat{1,6},'1')
            pin = 1;
        elseif contains(dat{1,6},'2')
            pin = 2;
        elseif contains(dat{1,6},'3')
            pin = 3;
        elseif contains(dat{1,6},'4')
            pin = 4;
        end
        events(i,3) = {pin};
    end
    
    for g = 1:length(events)
        events(g,2) = num2cell((cell2mat(events(g,1))-firstFileStartTime(1,1))*samplesMs); %sample number of the 
    end                             %sample that was being recorded when the event occured 
end
end

