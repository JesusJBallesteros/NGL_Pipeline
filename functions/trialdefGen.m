function [events, trialdef, eventdef] = trialdefGen(EventRecord, opt, varargin)
% Testing in experiments with Deuteron block format with a text file
% generated from the software log. This log NEEDS to be saved and placed
% with the raw session data manually (for now).
%
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

% Jesus 11.07.2024
if ~isfield(opt,'trEvents'),        opt.trEvents            = [];                   end

%% 01 Check inputs
if nargin > 2,  useevents = varargin{1};
else,           useevents = 0;
end

events = [];
trialdef = [];
eventdef = [];

switch useevents
    case 0
        %% 00 Sanity check for matching start/end events
        % Index of events equal to the defined trial start and trial end events.
        idx.start   = find(EventRecord.EventType == opt.eventdef.itiOn); 
        idx.end     = find(EventRecord.EventType == opt.eventdef.end1 | ...
                           EventRecord.EventType == opt.eventdef.end2 | ...
                           EventRecord.EventType == opt.eventdef.end3);
        
        if ~(length(idx.start)==length(idx.end)) % matching start-end events
            warning('A mismatch between number of start/end trials found. Trying to fix it.')
            % Possible sources of start-end mismatch:
            if any(idx.end(idx.end<idx.start(1)))
                % trialend events BEFORE first trialstart. Possible error ending
                % a previous session, leaving the pins in a different state than 
                % the expected [1 1 0 0], generating succesive arbitrary events 
                % until a point where the preIni state is enforced. 
                % Solution, remove all events before first star trial event.
                EventRecord.EventNumber(1:idx.start(1)-1)   = [];
                EventRecord.EventType(1:idx.start(1)-1)     = [];
                EventRecord.TimeStamp(1:idx.start(1)-1)     = [];
                EventRecord.TimeMsFromMidnight(1:idx.start(1)-1) = [];
                EventRecord.TimeSource(1:idx.start(1)-1)    = [];
                EventRecord.Details(1:idx.start(1)-1)       = [];
                % Possible FIX to recover these initial trials? Assume firs sent event
                % is start trial. MANUAL CHECK!
    	        warning('Events before first start trial removed. Check if these trials are recoverable.')
        
            elseif any(idx.start(idx.start>idx.end(end)))
                % This is a lonely trial start with no apparent end. Error
                % at session level or at event reading? Get rid of this
                % lonely last trial.
                EventRecord.EventNumber(idx.start(end):end)   = [];
                EventRecord.EventType(idx.start(end):end)     = [];
                EventRecord.TimeStamp(idx.start(end):end)     = [];
                EventRecord.TimeMsFromMidnight(idx.start(end):end) = [];
                EventRecord.TimeSource(idx.start(end):end)    = [];
                EventRecord.Details(idx.start(end):end)       = [];
            end
            
            % re-run idexing due to cover the changes
            idx.start   = find(EventRecord.EventType==opt.eventdef.itiOn); 
            idx.end     = find(EventRecord.EventType==opt.eventdef.end1 | ...
                               EventRecord.EventType==opt.eventdef.end2 | ...
                               EventRecord.EventType==opt.eventdef.end3);
        end
        
        % 01 Relativize timestamps to session start keeping it in msec
        EventRecord.TimeMsFromMidnight = (EventRecord.TimeMsFromMidnight - EventRecord.TimeMsFromMidnight(1));
        
        % 02 Convert relativized timestamps to SECONDS
        EventRecord.TimeSecFromMidnight = EventRecord.TimeMsFromMidnight/1000;
        
        % 03 Find trial start/end times using given definitions
        % Index of events equal to the defined trial start and trial end events.
        idx.start   = find(EventRecord.EventType==opt.eventdef.itiOn); 
        idx.end     = find(EventRecord.EventType==opt.eventdef.end1 | ...
                            EventRecord.EventType==opt.eventdef.end2 | ...
                            EventRecord.EventType==opt.eventdef.end3);
        
        % Now take those index time values
        trialstarts = EventRecord.TimeMsFromMidnight(idx.start); % get corresponding timestamps.
        trialends = EventRecord.TimeMsFromMidnight(idx.end); % get corresponding timestamps.

    case 1
        % 00 Find trials with specific event combinations
        % Index of events equal to the defined trial start and trial end events.
        idx.start   = cellfun(@(x) any(x==opt.eventdef.startON), EventRecord.code, 'UniformOutput', 1);
        idx.end     = cellfun(@(x) any(x==opt.eventdef.chc),   EventRecord.code, 'UniformOutput', 1);
        
        % Keep only those that fulfill all conditions above
        EventRecord.code = EventRecord.code(idx.start & idx.end);
        EventRecord.time = EventRecord.time(idx.start & idx.end);

        % Re-index to keep track of only valid ones
        idx.start   = cellfun(@(x) any(x==opt.eventdef.startON), EventRecord.code, 'UniformOutput', 1);
        idx.end     = cellfun(@(x) any(x==opt.eventdef.chc),   EventRecord.code, 'UniformOutput', 1);

        % Now take those index time values
        trialstarts = cellfun(@(x) x(1), EventRecord.time, 'UniformOutput', 1); % get corresponding timestamps.
        trialends   = cellfun(@(x) x(end), EventRecord.time, 'UniformOutput', 1); % get corresponding timestamps.

end

%% 02 Safety check, in case of unsolved problem.
assert(length(trialstarts)==length(trialends),'Mismatch between number of start/end events unsolved!')

% If OK, use either as a reliable count for number of trials
ntrials = length(idx.start); % count trial starts.

%% 03 Create trialdef variables for FieldTrip. In MILISECONDS
% Check options and prepare given events to align trial times to.
opt.alignto = events2align(opt);

% the field 't0' is an cell array of decimal values and their char arrays.
% Then, create a 'trialdef' xxx array where 
% Nx3, where columns are 'trial start time', 'trial end time' and 'offset to zero'.
trialdef = cell(2,size(opt.alignto,1));

% Go over every event and create the required trialdef aligned for that event
switch useevents
    case 0
        for i=1:size(opt.alignto,1)
            trialdef{1,i} = opt.alignto{i,1};
            
            idx = find(EventRecord.EventType==opt.alignto{i,2}); 
            trialdef{2,i}(:,1) = trialstarts;
            trialdef{2,i}(:,2) = trialends; 
            trialdef{2,i}(:,3) = EventRecord.TimeMsFromMidnight(idx); % get corresponding timestamps.
        end
    
    case 1
        for p = 2:ntrials
            tsum = trialends(p-1);
                trialstarts(p) = trialstarts(p)+tsum;
                trialends(p) = trialends(p)+tsum;
        end

        for i=1:size(opt.alignto,1)
            trialdef{1,i} = opt.alignto{i,1};
            
            % Check all trials for alignment event, allow no uniform. 
            idx.align = cellfun(@(x) find(x==opt.eventdef.startON), EventRecord.code, 'UniformOutput', 0);
            
            % safe check that all trials contain alignment event
            chck = cellfun(@(x) isempty(x), idx.align, 'UniformOutput', 1);
            idx.align = cell2mat(idx.align); % cell 2 mat

            % Retrieve the times compatible with alignment event
            trialdef{2,i}(:,1) = trialstarts(~chck);
            trialdef{2,i}(:,2) = trialends(~chck);
            for p = 2:ntrials-(sum(chck))
                trialdef{2,i}(p,3) = EventRecord.time{p}(idx.align(p))+trialdef{2,i}(p-1,2);
            end
        end
end

% Moved down % Find special events (ITI events as treatments, tutors, etc), if any.
% % Find the time of occurrence, as that should be enought to stablish them
% % in any further analisys (trials < t=x vs. trials > t=x)
% if ~isempty(opt.trEvents)
%     for i=1:length(opt.trEvents)
%         tr_idx.(opt.trEvents{i}) = find(EventRecord.EventType==opt.eventdef.(opt.trEvents{i}));
%         tr_idx.(opt.trEvents{i}) = EventRecord.TimeMsFromMidnight(tr_idx.(opt.trEvents{i}));
%     end
% end

%% 04 Create an NGl-standard event structure. IN SECONDS
if ~useevents
    % For each requested time alignment, an 'events.(event_align)' structure with fields
    %   .code {numtrials,1}, in decimal values as the standard from first event belonging to the trial t to the last one.
    %   .time {numtrials,1}, in SECONDS, aligned to a cero time fixed to an specific event (normally, itiOn).
    for i=1:size(opt.alignto,1)
        events.(opt.alignto{i,1}) = [];
        for t = 1:ntrials
            % Grab all timestamps between time of start and time of end (inclusive)
            trialstamps = EventRecord.TimeSecFromMidnight(EventRecord.TimeSecFromMidnight >= trialdef{2,i}(t,1)/1000 & ...
                                                         EventRecord.TimeSecFromMidnight <= trialdef{2,i}(t,2)/1000);
            % Relativize trial timestamps to alignment offset
            trialstamps = trialstamps - trialdef{2,i}(t,3)/1000; 
        
            % Grab all events ocurring between time of start and time of end (inclusive)
            trialevents = EventRecord.EventType(EventRecord.TimeSecFromMidnight >= trialdef{2,i}(t,1)/1000 & ...
                                                EventRecord.TimeSecFromMidnight <= trialdef{2,i}(t,2)/1000);
    
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
            % Find special events (ITI events as treatments, tutors, etc), if any.
            % Find the time of occurrence, as that should be enought to stablish them
            % in any further analisys (trials < t=x vs. trials > t=x)
            tr_idx.(opt.trEvents{i}) = find(EventRecord.EventType==opt.eventdef.(opt.trEvents{i}));
            tr_idx.(opt.trEvents{i}) = EventRecord.TimeMsFromMidnight(tr_idx.(opt.trEvents{i}));

            % Regular treatments
            events.(opt.trEvents{i}).code{1} = opt.eventdef.(opt.trEvents{i});
            events.(opt.trEvents{i}).time{1} = tr_idx.(opt.trEvents{i})/1000;
    
            % add special .trial field for easy indexing at e.g. plotting
            for j=1:length(events.(opt.trEvents{i}).time{1})
                events.(opt.trEvents{i}).trial{1}(j) = sum(tr_idx.(opt.trEvents{i})(j) > trialdef{2,1}(:,2));
            end
        end
    end
end

%% 07 Outputs
opt.eventdef.t0 = opt.alignto;
eventdef = opt.eventdef; % To keep track of definitions used

% trialdef % trial definition array for FieldTrip
% events % Event structure as NGL standard for spike processing

end