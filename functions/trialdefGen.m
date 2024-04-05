function [events, trialdef, eventdef] = trialdefGen(EventRecord, opt)
% Testing in experiments with Deuteron block format with a text file
% generated from the software log. This log NEEDS to be saved and placed
% with the raw session data manually (for now).
%
% The Event system and descritipions are based on a probably-to-be standard, as
% Deuteron current capabilities include reading single pin changes, limited
% to four input pins only. Therefore we are restricted to a sucession of
% 4-pin states achieved by single-bit changes at a time. This makes for a
% total of 16 possible states (decimal values 0-15).
%
% INPUT: EventRecords: struct with all events recorded during session.
%           EventNumber (double)
%           EventType (string)  
%           TimeStamp (string)
%           TimeMsFromMidnight (double)
%           TimeSource (string)
%           Details (string)
%        opt: struct with optional field 'eventdef' and all necessary
%               subfields to define all possible eventcodes, as well as the event
%               code use to align to time zero.
%       
% OUTPUT: events: struct with fields
%           code {numtrials,1}, in decimal values as the standard from first event
%                 belonging to the trial t to the last one.
%           time {numtrials,1}, in seconds, aligned to a cero time fixed to an 
%                 specific event (normally, itiOn).
%         trialdef: array (ntrials,3) columns being [trial startTime, trial endTime, trial ZeroTime]
%         eventdef: the event definitions used to create trials, 
%                   either defaulted or the ones given by the user.

% Jesus. 05.03.2024

%% 01.1 Relativize timestamps to session start keeping it in msec
EventRecord.TimeMsFromMidnight = (EventRecord.TimeMsFromMidnight - EventRecord.TimeMsFromMidnight(1));

%% 01.2 Find trial start/end times using given definitions, count number of trials found.
% Index of events equal to the defined trial start event.
idx = find(EventRecord.EventType==opt.eventdef.itiOn); 
trialstarts = EventRecord.TimeMsFromMidnight(idx); % get corresponding timestamps.

% Now index the end of the trials
idx = find(EventRecord.EventType==opt.eventdef.end1 | ...
           EventRecord.EventType==opt.eventdef.end2 | ...
           EventRecord.EventType==opt.eventdef.end3);
trialends = EventRecord.TimeMsFromMidnight(idx); % get corresponding timestamps.

% Safety check, are trialstarts and trialends equal?
assert(length(trialstarts)==length(trialends),'Mismatch found between number of trialStatrt and trialEnd events!')

% I OK, use either as a reliable count for number of trials
ntrials = length(trialstarts); % count trial starts.

% Now index the events we want to align the trials as t0
idx = find(EventRecord.EventType==opt.eventdef.t0); 
trialt0 = EventRecord.TimeMsFromMidnight(idx); % get corresponding timestamps.

%% 01.3 Create a trial array  for Fieldtrip.
% Create a trial array 'trialdef' for FielTrip (Nx3, where columns are:
% 'trial start time', 'trial end time' and 'offset to zero').
trialdef      = nan(ntrials,3);
    trialdef(:,1) = trialstarts;
    trialdef(:,2) = trialends; 
    trialdef(:,3) = trialt0;

%% 02.1 Convert relativized timestamps to SECONDS
EventRecord.TimeSecFromMidnight = EventRecord.TimeMsFromMidnight/1000;

%% 02.2 Find trial start/end times using given definitions, count number of trials found.
% Index of events equal to the defined trial start event.
idx = find(EventRecord.EventType==opt.eventdef.itiOn); 
trialstarts = EventRecord.TimeSecFromMidnight(idx); % get corresponding timestamps IN SEC.

% use this as a reliable marker for number of trials
ntrials = length(trialstarts); % count trial starts.

% Now index the end of the trials
idx = find(EventRecord.EventType==opt.eventdef.end1 | ...
           EventRecord.EventType==opt.eventdef.end2 | ...
           EventRecord.EventType==opt.eventdef.end3);
trialends = EventRecord.TimeSecFromMidnight(idx); % get corresponding timestamps, IN SEC.

% Now index the events we want to align the trials as t0
idx = find(EventRecord.EventType==opt.eventdef.t0); 
trialt0 = EventRecord.TimeSecFromMidnight(idx); % get corresponding timestamps, in SEC.

%% 02.3 Create an NGl-standard event structure. With fields:
%   .code {numtrials,1}, in decimal values as the standard from first event belonging to the trial t to the last one.
%   .time {numtrials,1}, in SECONDS, aligned to a cero time fixed to an specific event (normally, itiOn).
events = struct('code',[],'time',[]);
for t = 1:ntrials
    % Grab all timestamps between time of start and time of end (both inclusive)
    trialstamps = EventRecord.TimeSecFromMidnight(EventRecord.TimeSecFromMidnight>=trialstarts(t) & ...
                                                 EventRecord.TimeSecFromMidnight<=trialends(t));
    % Relativize trial timestamps to t0
    trialstamps = trialstamps - trialt0(t); 

    % Grab all events ocurring between time of start and time of end (both inclusive)
    trialevents = EventRecord.EventType(EventRecord.TimeSecFromMidnight>=trialstarts(t) & ...
                                        EventRecord.TimeSecFromMidnight<=trialends(t));

    % Insert into the proper structure to be output.
    events.code{t,1} = trialevents; 
    events.time{t,1} = trialstamps; 
end

%% 03. Outputs
eventdef = opt.eventdef; % To keep track of definitions used
% trialdef % trial definition array for FieldTrip
% events % Event structure as NGL standard for spike processing

end