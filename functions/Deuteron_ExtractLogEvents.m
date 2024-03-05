function EventRecord = Deuteron_ExtractLogEvents(opt)
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

% Jesus. 24.01.2024

if ~isfield(opt,'delimiters'),  opt.delimiters  = {',','='};    end
if ~isfield(opt,'outputas'),    opt.outputas    = 'string';     end
if ~isfield(opt,'iniPins'),     opt.iniPins     = [1 1 0 0];    end

%% Read the text file containing the Deuteron log and output a matrix using
% the given delimiters. By default it should output a matrix where columns are:
% [local time, msec after midnight, SpikeLog SN, local HH:MM:SS.MSEC, InputCh, InputState, Port]
logevents = readmatrix('logevents.txt', 'OutputType', opt.outputas, 'Delimiter', opt.delimiters);

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
% stateLog(1,:) = []; % remove initial state
stateLog = int2str(stateLog);

% %% As a final account for active channels, we use the explicit log about it
% % that Deuteron provides with every new file created while recording.
% filestarted = find(strcmp({EventRecord.EventType}, 'File started')==1); % Find the log for a new file started.
% geninfo = split(EventRecord(filestarted(1)).Details, ";"); % Split the text contained in Details using semicolons.
% geninfo = regexp(geninfo,'\d*','Match'); % Match the general expression '\d*'.
% opt.numChannels = str2double(geninfo{3}); % Transform the 3rd field (hardcoded) into double.
% if isempty(opt.channelOrder)
%     opt.channelOrder = 1:1:opt.numChannels; % Order channels as incremental ordinals. (TODO: this? perhaps match the Deuteron map?)
% end

%% Place extracted information into a proper EventRecord
EventRecord.EventNumber = double(1:1:length(stateLog))';
EventRecord.EventType = single(bin2dec(stateLog));
EventRecord.TimeStamp = string(ts); % Convert to string array
EventRecord.TimeMsFromMidnight = tsmsec;
EventRecord.TimeSource = nan(length(stateLog),1);
EventRecord.Details = nan(length(stateLog),1);

end