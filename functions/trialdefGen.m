function [events, trialdef, eventdef, EventRecord] = trialdefGen(EventRecord, opt, varargin)
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
%       24.04.2025: timebreak fix implemented

if ~isfield(opt,'trEvents'),        opt.trEvents            = [];                   end

%% 01 Check inputs
if nargin > 2,  useevents = varargin{1};
                trialdef = varargin{2};
else,           useevents = {};
                trialdef = [];
end

events = [];
% eventdef = [];

if isempty(useevents)
    %% 00 Sanity check for matching start/end events
    % Index of events equal to the defined trial start and trial end events.
    idx.start   = find(EventRecord.EventType == opt.eventdef.itiOn); 
    idx.end     = find(EventRecord.EventType == opt.eventdef.end1 | ...
                       EventRecord.EventType == opt.eventdef.end2 | ...
                       EventRecord.EventType == opt.eventdef.end3);
    
    if ~(length(idx.start)==length(idx.end)) % matching start-end events
        warning('A mismatch between number of start/end trials found.')
        if exist(fullfile(opt.behavFiles,"EventRecord.mat"),"file") 
            load(fullfile(opt.behavFiles,"EventRecord.mat"), 'EventRecord');
            warning('A fixed EventRecord variable found.')
        else
            warning('Trying to fix it.')
            % Possible sources of start-end mismatch:
            if any(idx.end(idx.end<idx.start(1)))
                % trialend events BEFORE first trialstart. Possible error ending
                % a previous session, leaving the pins in a different state than 
                % the expected [1 1 0 0], generating succesive arbitrary events 
                % until a point where the preIni state is enforced. 
                % Solution, remove all events before first star trial event.
                EventRecord.EventNumber(4639:4650)   = []; % 2:5 and 4639:4650 % 2993:3008
                EventRecord.EventType(4639:4650)     = [];
                EventRecord.TimeStamp(4639:4650)     = [];
                EventRecord.TimeMsFromMidnight(4639:4650) = [];
                EventRecord.TimeSource(4639:4650)    = [];
                EventRecord.Details(4639:4650)       = [];
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
        end
        
        % re-run idexing to recover the changes
        idx.start   = find(EventRecord.EventType==opt.eventdef.itiOn); 
        idx.end     = find(EventRecord.EventType==opt.eventdef.end1 | ...
                           EventRecord.EventType==opt.eventdef.end2 | ...
                           EventRecord.EventType==opt.eventdef.end3);
    
        % 000. If there was a timebreak, relativize it to the first timestamp
        % Also, all times after the break need to be adjusted for the actual
        % time passed during the delay, erasing it in terms of recording time
        if isfield(EventRecord,'TimeBreak')
            if ~isempty(EventRecord.TimeBreak{1,2})
                EventRecord.TimeBreak{1,2} = EventRecord.TimeBreak{1,2}-EventRecord.TimeMsFromMidnight(1);
                
                tbreakdur = EventRecord.TimeBreak{1,2}(2) - EventRecord.TimeBreak{1,2}(1);
                tbidx = EventRecord.TimeMsFromMidnight > EventRecord.TimeBreak{1,2}(2);
                EventRecord.TimeMsFromMidnight(tbidx) = EventRecord.TimeMsFromMidnight(tbidx) - tbreakdur;
            end
        end
    end
    
    % 01 Relativize timestamps to session start keeping it in msec
    EventRecord.TimeMsFromMidnight = (EventRecord.TimeMsFromMidnight - EventRecord.TimeMsFromMidnight(1));
        
    % 02 Find trial start/end times using given definitions
    % Index of events equal to the defined trial start and trial end events.
    idx.start   = find(EventRecord.EventType==opt.eventdef.itiOn); 
    idx.end     = find(EventRecord.EventType==opt.eventdef.end1 | ...
                        EventRecord.EventType==opt.eventdef.end2 | ...
                        EventRecord.EventType==opt.eventdef.end3);
        
    % Now take those index time values
    trialstarts = EventRecord.TimeMsFromMidnight(idx.start); % get corresponding timestamps.
    trialends = EventRecord.TimeMsFromMidnight(idx.end); % get corresponding timestamps.

    % 03 Check for trial length consistency
    triallengths = trialends-trialstarts;
    Avtriallength = median(triallengths); 
    if sum(triallengths > Avtriallength*1.1)==1
        NotTrial = find(triallengths > Avtriallength*1.1);
        trialstarts(NotTrial) = [];
        trialends(NotTrial) = [];

        rmvtrial(1) = find(EventRecord.TimeMsFromMidnight==EventRecord.TimeMsFromMidnight(idx.start(NotTrial))==1);
        rmvtrial(2) = find(EventRecord.TimeMsFromMidnight==EventRecord.TimeMsFromMidnight(idx.end(NotTrial))==1);

        % Remove events in between to eliminate its trace
        EventRecord.EventNumber(rmvtrial(1):rmvtrial(2)) = [];
        EventRecord.EventType(rmvtrial(1):rmvtrial(2))  = [];
        EventRecord.TimeStamp(rmvtrial(1):rmvtrial(2))  = []; % Convert to string array
        EventRecord.TimeMsFromMidnight(rmvtrial(1):rmvtrial(2)) = [];
        EventRecord.TimeSource(rmvtrial(1):rmvtrial(2)) = [];
        EventRecord.Details(rmvtrial(1):rmvtrial(2))    = [];

        warning('Exactly one trial have been found unconsistently lenghty, and has been excluded.')
    elseif sum(triallengths > Avtriallength*1.1)>1
        warning('Several trials have unconsistent length!')
    else
        disp('All trials are consistent in duration.')
    end

    % 05 Convert relativized timestamps to SECONDS
    EventRecord.TimeSecFromMidnight = EventRecord.TimeMsFromMidnight/1000;

% else
    % %% 00 Find trials with specific event combinations
    % % Index of events equal to the defined trial start and trial end events.
    % idx.start   = cellfun(@(x) any(x==opt.eventdef.itiOn), EventRecord.code, 'UniformOutput', 1);
    % idx.end     = cellfun(@(x) any(x==opt.eventdef.end1 | ...
    %                                x==opt.eventdef.end2 | ...
    %                                x==opt.eventdef.end3 ),  EventRecord.code, 'UniformOutput', 1);
    % 
    % % Keep only those that fulfill all conditions above
    % EventRecord.code = EventRecord.code(idx.start & idx.end);
    % EventRecord.time = EventRecord.time(idx.start & idx.end);
    % 
    % % Re-index to keep track of only valid ones
    % idx.start   = cellfun(@(x) any(x==opt.eventdef.startON), EventRecord.code, 'UniformOutput', 1);
    % idx.end     = cellfun(@(x) any(x==opt.eventdef.chc),   EventRecord.code, 'UniformOutput', 1);
    % 
    % % Now take those index time values
    % trialstarts = cellfun(@(x) x(1), EventRecord.time, 'UniformOutput', 1); % get corresponding timestamps.
    % trialends   = cellfun(@(x) x(end), EventRecord.time, 'UniformOutput', 1); % get corresponding timestamps.

    %% 02 Safety check, in case of unsolved problem.
    assert(length(trialstarts)==length(trialends),'Mismatch between number of start/end events unsolved!')
    
    % If OK, use either as a reliable count for number of trials
    ntrials = length(trialstarts); % count trial starts.

    %% 03 Create trialdef variables. In MILISECONDS
    % Check options and prepare given events to align trial times to.
    opt.alignto = events2align(opt);
    
    % the field 't0' is an cell array of decimal values and their char arrays.
    % Then, create a 'trialdef' xxx array where 
    % Nx3, where columns are 'trial start time', 'trial end time' and 'offset to zero'.
    trialdef = cell(2,size(opt.alignto,1));
else
    % 01 Find trial start/end times using given definitions
    % Index of events equal to the defined trial start and trial end events.
    idx.start   = find(EventRecord.EventType==opt.eventdef.itiOn); 
    idx.end     = find(EventRecord.EventType==opt.eventdef.end1 | ...
                    EventRecord.EventType==opt.eventdef.end2 | ...
                    EventRecord.EventType==opt.eventdef.end3);
    
    % 01 Relativize timestamps to session start keeping it in msec
    EventRecord.TimeMsFromMidnight = (EventRecord.TimeMsFromMidnight - EventRecord.TimeMsFromMidnight(1));
    
    % 02 Convert relativized timestamps to SECONDS
    EventRecord.TimeSecFromMidnight = EventRecord.TimeMsFromMidnight/1000;

    % Now take those index time values
    trialstarts = EventRecord.TimeMsFromMidnight(idx.start); % get corresponding timestamps.
    trialends = EventRecord.TimeMsFromMidnight(idx.end); % get corresponding timestamps.

    ntrials = length(idx.start); % count trial starts.

    opt.alignto = events2align(opt);
end

if isempty(useevents)
    %% Go over every event and create the required trialdef aligned for that event
    for i=1:size(opt.alignto,1)
        correction = 0;
        trialdef{1,i} = opt.alignto{i,1};
        
        if i > 1
           trialdef{2,i} = nan(size(trialdef{2,1}));
        end

        idx = find(EventRecord.EventType==opt.alignto{i,2});
            % if strcmp(opt.alignto{i,1},'bhv') % TODO. Has to be used ONLY for bhv2 in S3-Extintion Arena
            %     idx(EventRecord.EventType(idx-1) ~= str2double(opt.alignto{i,3})) = [];
            %     if str2double(opt.alignto{i,3})==2, correction = 1000; end % Fix for bhv-rwd in S3-Extintion Arena
            % end

        % start times
        trialdef{2,i}(:,1) = trialstarts-opt.addtime;
            if trialdef{2,i}(1,1) < 0
                trialdef{2,i}(1,1) = trialdef{2,i}(1,1)+opt.addtime; 
            end
       
        % end times
        trialdef{2,i}(:,2) = trialends+opt.addtime;
        
        % zero times
        if size(trialdef{2,i},1) == size(idx,1)
            trialdef{2,i}(:,3) = EventRecord.TimeMsFromMidnight(idx)-correction;
        elseif size(trialdef{2,i},1) ~= size(idx,1)
            trl = 1;
            for td = 1:ntrials
                tmps = EventRecord.TimeMsFromMidnight(idx(trl))-correction;
                if tmps > trialdef{2,i}(td,1) && tmps < trialdef{2,i}(td,2)
                    trialdef{2,i}(td,3) = tmps;
                    trl = trl + 1;
                % else TODO
                % To find out if some trials are missing in trialdef, that
                % are in EvenrRecord (idx). This could mean the trialdef
                % generation has been modified somewhere here, due to
                % mismatch on start-end events.
                %     trialdef{2,i}(td,3) = nan;
                end
            end
        end
    end  
else
    %% Go over every New Event and create the required alignment
    oldN = length(trialdef);

    for i = 1:size(opt.newEvent,1)
        if isempty(opt.newEvent{i,1}), continue, end
        trialdef{1,oldN+i} = 'bhv';%opt.newEvent{i,1};
        trialdef{2,oldN+i} = nan(size(trialdef{2,1}));

        idx = find(EventRecord.EventType == opt.eventdef.bhv);%(opt.newEvent{i,1})); 
        
        % Commonly we use a single 'bhv' for any response, given Deuteron
        % limitations. This helps differentiating sequences of behavioral
        % responses.
        if size(opt.newEvent,2)<2, opt.newEvent{i,2} = ''; end 
        if strcmp(opt.newEvent{i,2},'1') 
            stim = opt.alignto(ismember(opt.alignto(:,1),'stimOn2'));
            idx(EventRecord.EventType(idx-1)==opt.eventdef.(stim{1})) = [];
        elseif strcmp(opt.newEvent{i,2},'2')
            stim = opt.alignto(ismember(opt.alignto(:,1),'stimOn1'));
            idx(EventRecord.EventType(idx-1)==opt.eventdef.(stim{1})) = [];
        % if there would be more than 2, keep adding accordingly
        % elseif strcmp(opt.newEvent{i,2},'3')
        % ...
        end
        
        % start times
        trialdef{2,oldN+i}(:,1) = trialstarts-opt.addtime;
       
        % end times
        trialdef{2,oldN+i}(:,2) = trialends+opt.addtime;
        
        % zero times
        if size(trialdef{2,oldN+i},1) == size(idx,1)
            trialdef{2,oldN+i}(:,3) = EventRecord.TimeMsFromMidnight(idx); 
        elseif size(trialdef{2,oldN+i},1) ~= size(idx,1)
            trl = 1;
            for td = 1:size(trialdef{2,1},1)
                if strcmp(opt.newEvent{i,2},'2')
                    tmps = EventRecord.TimeMsFromMidnight(idx(trl))-1000; % -1000 is an Extintion arena FIX for 079's #1-11)
                else
                    tmps = EventRecord.TimeMsFromMidnight(idx(trl));
                end
                if tmps > trialdef{2,oldN+i}(td,1) && tmps < trialdef{2,oldN+i}(td,2)
                    trialdef{2,oldN+i}(td,3) = tmps;
                    trl = trl + 1;
                % else
                %     trialdef{2,oldN+i}(td,3) = nan;
                end
            end
        end

        % % Check all trials for alignment event, allow no uniform. 
        % idx.align = cellfun(@(x) find(x==opt.eventdef.startON), EventRecord.code, 'UniformOutput', 0);
        % 
        % % safe check that all trials contain alignment event
        % chck = cellfun(@(x) isempty(x), idx.align, 'UniformOutput', 1);
        % idx.align = cell2mat(idx.align); % cell 2 mat
        % 
        % % Retrieve the times compatible with alignment event
        % trialdef{2,i}(:,1) = trialstarts(~chck)-opt.addtime;
        % trialdef{2,i}(:,2) = trialends(~chck)+opt.addtime;
        % for p = 2:ntrials-(sum(chck))
        %     trialdef{2,i}(p,3) = EventRecord.time{p}(idx.align(p))+trialdef{2,i}(p-1,2);
        % end
    end
end

% %% Trialdef double check for time jumps and its fixing
% timebreak = check_timebreaks(trialdef{2,1}); % {2,1}(:,3) should always be t0 for itiOn

%% 04 Create an NGl-standard event structure. IN SECONDS
if isempty(useevents)
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
else
    events = useevents;
    for i=1:size(opt.newEvent,1)
        % if isempty(opt.newEvent{i,1}), continue, end
        events.(opt.alignto{oldN+i,1}) = [];
        for t = 1:ntrials
            % Grab all timestamps between time of start and time of end (inclusive)
            trialstamps = EventRecord.TimeSecFromMidnight(EventRecord.TimeSecFromMidnight >= trialdef{2,oldN+i}(t,1)/1000 & ...
                                                         EventRecord.TimeSecFromMidnight <= trialdef{2,oldN+i}(t,2)/1000);
            
            % Grab all events ocurring between time of start and time of end (inclusive)
            trialevents = EventRecord.EventType(EventRecord.TimeSecFromMidnight >= trialdef{2,oldN+i}(t,1)/1000 & ...
                                                EventRecord.TimeSecFromMidnight <= trialdef{2,oldN+i}(t,2)/1000);
           
            % Relativize trial timestamps to alignment offset
            trialstamps = trialstamps - trialdef{2,oldN+i}(t,3)/1000;

            if strcmp(opt.newEvent{i,2},'2')
                trialstamps(trialstamps==1) = 0;
            end
    
            % Insert into the proper structure to be output.
            events.(opt.alignto{oldN+i,1}).code{t,1} = trialevents; 
            events.(opt.alignto{oldN+i,1}).time{t,1} = trialstamps; 
        end
    end
end

%% 07 Outputs
opt.eventdef.t0 = opt.alignto;
eventdef = opt.eventdef; % To keep track of definitions used

% trialdef % trial definition array for FieldTrip
% events % Event structure as NGL standard for spike processing

end