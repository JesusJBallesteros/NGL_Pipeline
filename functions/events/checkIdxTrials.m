function [idx, EventRecord] = checkIdxTrials(idx, EventRecord, opt)

    t = 1;    % try counter
    tr = 0;   % trial removed counter  
    % Evaluate number of starts and ends
    while numel(idx.start) ~= numel(idx.end)
        % With a mismatch,
        % Start warning the user and check if a finished file already
        % exists, probably already checked. Use it if so.
        warning('A mismatch between number of start/end trials found.')
        if exist(fullfile(opt.behavFiles,"EventRecord.mat"),"file") && t==1
            load(fullfile(opt.behavFiles,"EventRecord.mat"), 'EventRecord');
            warning('A fixed EventRecord variable found.\n')
            t = t+1; % first try
            continue
        end

        % Either the loaded file is not fixed, or there is no file, try to
        % fix it.

        % More starts than ends
        if numel(idx.end) < numel(idx.start)
            % Find an index of start events that are NOT preceded by an end event.
            chkList = find(~ismember(EventRecord.EventType(idx.start(2:end)-1), ...
                                [opt.eventdef.end1, opt.eventdef.end2, opt.eventdef.end3]))+1;
        
            % Point to those events in EventRecord.EventType to check
            chkList = idx.start(chkList);
            % We need to start from the end, to avoid general index changes. Make row
            chkList = flip(chkList');
        
            for i = chkList
                if ~ismember(EventRecord.EventType(i-1), ...
                            [opt.eventdef.end1, opt.eventdef.end2, opt.eventdef.end3])
                    if ~ismember(EventRecord.EventType(i-2), ...
                                [opt.eventdef.end1, opt.eventdef.end2, opt.eventdef.end3])

                        % if the i-2 event is indeed an end event, i-1 is probably a
                        % block change event. We need to keep it.
                        % We don't contemplate any case (yet) where we need to send 
                        % two intertrial event codes in INTAN, so it must be some
                        % sort of error, remove all until the next found end event.
                        lastend = find(ismember(EventRecord.EventType(1:i-1),...
                                        [opt.eventdef.end1, opt.eventdef.end2, opt.eventdef.end3]), ...
                                        1,"last");
                        tr = tr+1;
                        EventRecord.EventType(lastend+1:i-1) = [];
                        EventRecord.EventNumber(lastend+1:i-1) = [];
                        EventRecord.TimeStamp(lastend+1:i-1) = [];
                        EventRecord.TimeMsFromMidnight(lastend+1:i-1) = [];
                        EventRecord.TimeSource(lastend+1:i-1) = [];
                        EventRecord.Details(lastend+1:i-1) = [];
                    end
                end
            end

        % More ends than starts
        elseif numel(idx.end) > numel(idx.start)
            % Find an index of end events that are NOT followed by an start event.
            chkList = find(EventRecord.EventType(idx.end(1:end-1)+1)~=opt.eventdef.itiOn);

            % Point to those events in EventRecord.EventType to check
            chkList = idx.end(chkList);
            % We need to start from the end, to avoid general index changes. Make row
            chkList = flip(chkList');

            for i = chkList
                if EventRecord.EventType(i+1)~=opt.eventdef.itiOn
                    if EventRecord.EventType(i+2)~=opt.eventdef.itiOn
                        % if the i+2 event is indeed an start event, i+1 is probably a
                        % block change event. We need to keep it.
                        % We don't contemplate any case (yet) where we need to send 
                        % two intertrial event codes in INTAN, so it must be some
                        % sort of error, remove all until the next found end event.
                        nextstart = find(EventRecord.EventType(i+1:end)==opt.eventdef.itiOn, ...
                                       1,"first");
                        tr = tr+1;
                        EventRecord.EventType(i+1:nextstart) = [];
                        EventRecord.EventNumber(i+1:nextstart) = [];
                        EventRecord.TimeStamp(i+1:nextstart) = [];
                        EventRecord.TimeMsFromMidnight(i+1:nextstart) = [];
                        EventRecord.TimeSource(i+1:nextstart) = [];
                        EventRecord.Details(i+1:nextstart) = [];
                    end
                end
            end
        end

        % There is a chance that the problem is either an unproperly initiated
        % first trial or a unproperly finished last trial. Easy:
        if any(idx.end(idx.end<idx.start(1)))
            % trialend events BEFORE first trialstart. Possible error ending
            % a previous session, leaving the pins in a different state than 
            % the expected [1 1 0 0], generating succesive arbitrary events 
            % until a point where the preIni state is enforced. 
            % Solution, remove all events before first star trial event.
            EventRecord.EventNumber(1:idx.end(1))   = [];
            EventRecord.EventType(1:idx.end(1))     = [];
            EventRecord.TimeStamp(1:idx.end(1))     = [];
            EventRecord.TimeMsFromMidnight(1:idx.end(1)) = [];
            EventRecord.TimeSource(1:idx.end(1))    = [];
            EventRecord.Details(1:idx.end(1))       = [];
            % Possible FIX to recover these initial trials? Assume first sent event
            % is start trial. MANUAL CHECK!
            warning('Events before first trial start removed.')
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
            warning('Events after the last trial end removed.')
        end

        % re-evaluate numel(idxs) to escape while loop when fixed
        idx.start   = find(EventRecord.EventType == opt.eventdef.itiOn); 
        idx.end     = find(EventRecord.EventType == opt.eventdef.end1 | ...
                           EventRecord.EventType == opt.eventdef.end2 | ...
                           EventRecord.EventType == opt.eventdef.end3);
    end

    % No mismatch found or solved
    disp('Matching number of start/end events. Appears to be a good session.')
    if tr > 0
        fprintf('%d misshaped trials removed.', tr)
    end
end