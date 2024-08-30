function [events, trialdef, eventdef] = trialdefGen(EventRecord, opt)
% The Event system and descritipions are based on a probably-to-be standard, as
% Deuteron current capabilities include reading single pin changes, limited
% to four input pins only. Therefore we are restricted to a sucession of
% 4-pin states achieved by single-bit changes at a time. This makes for a
% total of 16 possible states (decimal integers 0:15).
%
% INPUT: EventRecords: struct with all events recorded during session.
%           EventNumber (double)
%           EventType (string)  
%           TimeStamp (string)
%           TimeMsFromMidnight (double) --> Converted to SECONDS
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

% Jesus 28.08.2024
if ~isfield(opt,'trEvents'),        opt.trEvents            = [];                   end

events      = [];
trialdef    = [];
eventdef    = [];

%% 00 Sanity check for matching start/end events
% 0.1 Index of events equal to the defined trial start and trial end events.
idx.start   = EventRecord.EventType == opt.eventdef.itiOn; 
    starts  = find(EventRecord.EventType == opt.eventdef.itiOn);
idx.end     = EventRecord.EventType == opt.eventdef.end1 | ...
                    EventRecord.EventType == opt.eventdef.end2 | ...
                    EventRecord.EventType == opt.eventdef.end3;
    ends    = find(EventRecord.EventType == opt.eventdef.end1 | ...
                    EventRecord.EventType == opt.eventdef.end2 | ...
                    EventRecord.EventType == opt.eventdef.end3);
idx.abort   = EventRecord.EventType == opt.eventdef.oms1;

if ~(sum(idx.start)==sum(idx.end)) % matching start-end events
    warning('A mismatch between number of start/end trials found. Trying to fix it.')
    % Possible sources of start-end mismatch:
    if any(ends(ends<starts(1)))
        % trialend events BEFORE first trialstart. Possible error ending
        % a previous session, leaving the pins in a different state than 
        % the expected [1 1 0 0], generating succesive arbitrary events 
        % until a point where the preIni state is enforced. 
        % Solution, remove all events before first star trial event.
        EventRecord.EventNumber(1:starts(1)-1)   = [];
        EventRecord.EventType(1:starts(1)-1)     = [];
        EventRecord.TimeStamp(1:starts(1)-1)     = [];
        EventRecord.TimeMsFromMidnight(1:starts(1)-1) = [];
        EventRecord.TimeSource(1:starts(1)-1)    = [];
        EventRecord.Details(1:starts(1)-1)       = [];
        % Possible FIX to recover these initial trials? Assume firs sent event
        % is start trial. MANUAL CHECK!
        warning('Events before first start trial removed. Check if these trials are recoverable.')
    end
    
    if any(starts(starts>ends(end)))
        % This could be a lonely trial start with no apparent end. Error
        % at session level or at event reading? Get rid of this last orphan trial.
        EventRecord.EventNumber(starts(end):end)   = [];
        EventRecord.EventType(starts(end):end)     = [];
        EventRecord.TimeStamp(starts(end):end)     = [];
        EventRecord.TimeMsFromMidnight(starts(end):end) = [];
        EventRecord.TimeSource(starts(end):end)    = [];
        EventRecord.Details(starts(end):end)       = [];
    end

    % re-run idexing due to cover the changes
    idx.start   = EventRecord.EventType == opt.eventdef.itiOn; 
        starts  = find(EventRecord.EventType == opt.eventdef.itiOn);
    
    idx.end     = EventRecord.EventType == opt.eventdef.end1 | ...
                        EventRecord.EventType == opt.eventdef.end2 | ...
                        EventRecord.EventType == opt.eventdef.end3;
        ends    = find(EventRecord.EventType == opt.eventdef.end1 | ...
                        EventRecord.EventType == opt.eventdef.end2 | ...
                        EventRecord.EventType == opt.eventdef.end3);
    
    idx.abort   = EventRecord.EventType == opt.eventdef.oms1;
end

% 0.2 Relativize timestamps to session start keeping it in msec
EventRecord.TimeMsFromMidnight = EventRecord.TimeMsFromMidnight - EventRecord.TimeMsFromMidnight(1);

% 0.3 Convert relativized timestamps to SECONDS
if strcmp(opt.ext,'fileperch'), fs=30000;
elseif strcmp(opt.ext,'DF1'), fs=32000;
end
EventRecord.TimeSecFromMidnight = EventRecord.TimeMsFromMidnight/fs;

% 0.5 Safety check, in case of unsolved problem.
assert(length(starts)==length(ends),'Mismatch between number of start/end events unsolved!')

%% 01 If OK, use either as a reliable count for number of trials
ntrials = length(starts); % count trial starts.

%% 02 Create trialdef variables for FieldTrip. In MILISECONDS
% Check options and prepare given events to align trial times to.
opt.alignto = events2align(opt);

% the field 't0' is an cell array of decimal values and their char arrays.
% Then, create a 'trialdef' xxx array where 
% Nx3, where columns are 'trial start time', 'trial end time' and 'offset to zero'.
trialdef = cell(2,size(opt.alignto,1));

% Now take those time stamps according to the indexes
tstarts = EventRecord.TimeMsFromMidnight(idx.start); % get corresponding timestamps.
tsends  = EventRecord.TimeMsFromMidnight(idx.end); % get corresponding timestamps.

% Go over every event to align to and create the required trialdef with its proper alignment
for i=1:size(opt.alignto,1)
    trialdef{1,i} = opt.alignto{i,1}; % event to align to
    trialdef{2,i} = nan(ntrials,3); % pre allocate all trials
    
    trialdef{2,i}(:,1)  = tstarts; % set trial start times
    trialdef{2,i}(:,2)  = tsends; % set trial end times
        
    idx.align = EventRecord.EventType==opt.alignto{i,2}; % index those events to align to
    if any(idx.align) % if any is found
        taligned = EventRecord.TimeMsFromMidnight(idx.align); % get corresponding timestamps.
        if sum(idx.align)==ntrials % if found in all trials
            trialdef{2,i}(:,3) = taligned; % write in block
        else % not an event in every trial
            tt = 1; % count events to locate
            for t=1:ntrials % on each trial, if the event
                if tt<=length(taligned) && ... % until its last presence
                    (tstarts(t)<taligned(tt) && taligned(tt)<tsends(t)) % belongs to the trial
                    trialdef{2,i}(t,3) = taligned(tt); %  write it
                    tt = tt+1; % and move on to next event to align to
                end
            end
        end
    end
end
    
% Find special events (ITI events as treatments, tutors, etc), if any.
% Find the time of occurrence, as that should be enought to stablish them
% in any further analisys (trials < t=x vs. trials > t=x)
if ~isempty(opt.trEvents)
    for i=1:length(opt.trEvents)
        idx.(opt.trEvents{i}) = EventRecord.EventType==opt.eventdef.(opt.trEvents{i});
        idx.(opt.trEvents{i}) = EventRecord.TimeMsFromMidnight(idx.(opt.trEvents{i}));
        if isempty(idx.(opt.trEvents{i})), idx.(opt.trEvents{i}) = nan; end
    end
end

%% 06 Create an NGl-standard event structure. IN SECONDS
% For each requested time alignment, an 'events.(event_align)' structure with fields
%   .code {numtrials,1}, in decimal values as the standard from first event belonging to the trial t to the last one.
%   .time {numtrials,1}, in SECONDS, aligned to a cero time fixed to an specific event (normally, itiOn).
for i=1:size(opt.alignto,1)
    events.(opt.alignto{i,1}) = [];

    for t = 1:ntrials
        % Grab all events ocurring between time of start and time of end (inclusive)
        trialevents = EventRecord.EventType(EventRecord.TimeSecFromMidnight >= trialdef{2,i}(t,1)/fs & ...
                                            EventRecord.TimeSecFromMidnight <= trialdef{2,i}(t,2)/fs);
        if ismember(opt.eventdef.(opt.alignto{i,1}), trialevents)
            % Grab all timestamps between time of start and time of end (inclusive)
            trialstamps = EventRecord.TimeSecFromMidnight(EventRecord.TimeSecFromMidnight >= trialdef{2,i}(t,1)/fs & ...
                                                         EventRecord.TimeSecFromMidnight <= trialdef{2,i}(t,2)/fs);
            
            % Relativize trial timestamps to alignment offset
            trialstamps = trialstamps - trialdef{2,i}(t,3)/fs; 

        else
            trialevents = nan;
            trialstamps = nan;

        end
        % Insert into the proper structure to be output.
        events.(opt.alignto{i,1}).code{t,1} = trialevents; 
        events.(opt.alignto{i,1}).time{t,1} = trialstamps; 

    end
end

% We can add here the special events for ITI treatments, if any. This
% special case will be simply coded as a single event code and the 
% timestamps when it appears, in seconds.
if ~isempty(opt.trEvents)
    for i=1:length(opt.trEvents)
        events.(opt.trEvents{i}).code{1} = opt.eventdef.(opt.trEvents{i});
        events.(opt.trEvents{i}).time{1} = idx.(opt.trEvents{i})/fs;

        % add special .trial field for easy indexing at e.g. plotting
        for j=1:length(events.(opt.trEvents{i}).time{1})
            if ~isnan(idx.(opt.trEvents{i})(j))
                events.(opt.trEvents{i}).trial{1}(j) = sum(idx.(opt.trEvents{i})(j) > trialdef{2,1}(:,2))+1;
            else
                events.(opt.trEvents{i}).trial{1}(j) = nan;
            end
        end
    end
end

eventdef = opt.eventdef; % To keep track of definitions used

end