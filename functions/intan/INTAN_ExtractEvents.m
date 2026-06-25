function EventRecord = INTAN_ExtractEvents(input, opt)
% INTAN_ExtractEvents  Read digital-input pin states and convert to EventRecord.
%
% PURPOSE:
%   Reads all board-DIGITAL-IN*.dat files from the current session folder,
%   detects rising and falling pin transitions (any pin changes), applies a
%   debounce window (28 samples, hardcoded) to suppress glitches, and 
%   converts each stable 16-bit pin state to a decimal event code using
%   binvec2dec (LSB-first). Produces an EventRecord struct in the same format
%   as Deuteron_ExtractEvents for downstream compatibility.
%   A security check at the end looks for the first 0/8 appearance. If
%   neither are the first event, all other events until the first 8
%   are removed.
%   A second repair walks every itiOn (0) and, when the immediately
%   preceding event is NOT preIni (8), inserts a synthetic 8 at
%   ts_0 - 30 samples — provided that 30-sample window is otherwise
%   empty. Covers the case where the behavioural task drops preIni
%   either entirely or on individual trials. The number of orphan
%   zeros, repaired zeros and unrepairable zeros is reported in a
%   single warning (NGL:INTAN_ExtractEvents:insertedPreIni).
%
% USAGE:
%   EventRecord = INTAN_ExtractEvents(input, opt)
%
% INPUTS:
%   input  - struct with input.run and
%              input.sessions(input.run(1)).info.amplifier_sample_rate
%   opt    - struct with opt.PathRaw (session raw-data folder)
%
% OUTPUT:
%   EventRecord - struct with fields:
%     .EventType           (nEvents × 1 double) decimal event codes (0–65535)
%     .EventNumber         (nEvents × 1 double) sequential event index
%     .TimeStamp           (nEvents × 1 double) sample indices of transitions
%     .TimeMsFromMidnight  (nEvents × 1 double) timestamps in ms
%     .TimeSource          (nEvents × 1 double) NaN (INTAN has no source field)
%     .Details             (nEvents × 1 double) NaN
%     .TimeBreak           (1 × 2 cell) {[] []} placeholder (no breaks on INTAN)
%
% EVENT CODING:
%   All 16 digital-input pins are read simultaneously. The 4 LSBs encode the
%   16 reserved events (codes 0–15); higher pins encode project-specific events
%   (codes ≥ 16). See eventDefinitions.m for the full vocabulary.
%   Decimal 0 (itiOn, all pins low) is a valid event — it marks trial start.
%
% Last modified 25.06.2026 (Jesus) - added missing-preIni repair: every
%                                     itiOn (0) without an immediately
%                                     preceding preIni (8) gets a
%                                     synthetic 8 inserted 30 samples
%                                     before it (if the window is empty),
%                                     with a clear warning.

%% Defaults
pth     = opt.PathRaw; % folder for reading events
smpDel  = 28;       % an event code is read after smpDel since first pin change, for additional smpDel since (to catch instabilities)

%% Initialize
INfiles = string(ls(fullfile(pth,"board-DIGITAL-IN*")));
npins = numel(INfiles);
pins = 1:1:npins;

%% Read all digital IN
for i = pins
    fid = fopen(fullfile(pth, ['board-DIGITAL-IN-' sprintf('%02d',i) '.dat'])); % Point to file
    tmp = fread(fid, inf, 'uint16'); % Get file data into tmp
    fclose(fid); % Close pointer
    if i == 1 % Using the first file, allocate memory
        nsampl = length(tmp);
        dIn = zeros(nsampl, npins, 'uint8');
    end
    % Add the data to dIn
    dIn(1:nsampl, i) = tmp; 
end
clear tmp i fid

%% Convert to sample # and event-code
% Keep samples during which pin 1 and/or 2 are up (per default, before first task event '0' is sent).
startState = dIn(1,:); % Read pins states as recording starts 

% % pinsOff = find(any(dIn(:,1:npins)~=startState,2), 1, "first");  % find first pin change
% % if pinsOff ~= 1 % Only if is not already the first sample
% %     dIn(1:pinsOff, :) = []; % remove all samples until then
% % end

%% Find samples at which any pin changes
checksum = diff([int8(startState); dIn],1,1); % fixed to admit negative values, added starting state to keep nsampl
ts = find(any(checksum,2)); % First ts found should be task start (8 or 0)
clear checksum 

% check for pin changes too close to each other (inconsistencies, or delayed
% change of a single event)
for i = 2:size(ts,1)
    if ts(i)-ts(i-1) < smpDel % if pin changes are too close
        ts(i-1) = nan;   % The previous event (i-1) becomes NAN because we
                         % can assume it was a transitory state towards the 
                         % desired event (current, i)
    end
end

% clean up NANs, to not be considered as new events
ts(isnan(ts)) = [];

%% CONVERT all events
% convert each binary word to its corresponding decimal using the npins bits
EventType = nan(size(ts,1),1);

for i = 1:size(ts,1)
    % Convert binary pins to decimal, as sum over smpDel forward to catch inconsitencies
    EventType(i) = binvec2dec(sum(dIn(ts(i):ts(i)+smpDel,:))); % binary vector to decimal integer
end
clear dIn    

%% Place extracted information into a proper EventRecord
EventRecord.EventType           = double(EventType);
EventRecord.EventNumber         = double(1:1:length(EventType))';
EventRecord.TimeStamp           = ts; % Updated 16.02.2026 to always keep original timeStamps (samples)
EventRecord.TimeMsFromMidnight  = ts/(input.sessions(input.run(1)).info.amplifier_sample_rate/1000);
EventRecord.TimeSource          = nan(length(EventType),1);
EventRecord.Details             = nan(length(EventType),1);
EventRecord.TimeBreak           = {[] []};

if ~ismember(EventType(1),[0,8])
    firstEv = find(EventRecord.EventType==8 | EventRecord.EventType==0,1,"first");
    EventRecord.EventType(1:firstEv-1)           = [];
    EventRecord.EventNumber(1:firstEv-1)         = [];
    EventRecord.TimeStamp(1:firstEv-1)           = [];
    EventRecord.TimeMsFromMidnight(1:firstEv-1)  = [];
    EventRecord.TimeSource(1:firstEv-1)          = [];
    EventRecord.Details(1:firstEv-1)             = [];
    warning('While extracting events, first event was neither 0 nor 8. Every event before the first encounter was removed.\n')
end

%% Repair missing preIni (8) events.
% Convention enforced by the trial-state machine: every itiOn (0) must be
% IMMEDIATELY preceded by a preIni (8). Some recordings drop preIni
% entirely (the behavioural script never sent 8) or lose it on
% individual trials. Detect either case and patch by inserting a
% synthetic 8 event 30 samples before each orphan 0 — only when that
% 30-sample window is otherwise empty (no other event there) and the
% inserted timestamp would not fall before sample 1. Warn loudly with
% counts so the experimenter sees there was a problem AND that the
% pipeline tried to recover.
zeroIdx     = find(EventRecord.EventType == 0);
orphanZeros = [];
for k = 1:numel(zeroIdx)
    z = zeroIdx(k);
    if z > 1 && EventRecord.EventType(z-1) == 8
        continue
    end
    orphanZeros(end+1, 1) = z;
end

if ~isempty(orphanZeros)
    insertTs  = zeros(0, 1);
    skippedTs = zeros(0, 1);
    for k = 1:numel(orphanZeros)
        z      = orphanZeros(k);
        tsZero = EventRecord.TimeStamp(z);
        tsIns  = tsZero - 30;
        if tsIns < 1
            skippedTs(end+1, 1) = tsZero;
            continue
        end
        % Window (tsIns, tsZero) must hold no existing event.
        inWindow = (EventRecord.TimeStamp >  tsIns) ...
                 & (EventRecord.TimeStamp <  tsZero);
        if any(inWindow)
            skippedTs(end+1, 1) = tsZero;
            continue
        end
        insertTs(end+1, 1) = tsIns;
    end

    % Apply insertions: append, then re-sort all parallel arrays by
    % TimeStamp and renumber EventNumber.
    if ~isempty(insertTs)
        nNew       = numel(insertTs);
        sampleRate = input.sessions(input.run(1)).info.amplifier_sample_rate;

        EventRecord.EventType          = [EventRecord.EventType; repmat(8, nNew, 1)];
        EventRecord.TimeStamp          = [EventRecord.TimeStamp; insertTs];
        EventRecord.TimeMsFromMidnight = [EventRecord.TimeMsFromMidnight; insertTs / (sampleRate/1000)];
        EventRecord.TimeSource         = [EventRecord.TimeSource; nan(nNew, 1)];
        EventRecord.Details            = [EventRecord.Details; nan(nNew, 1)];

        [EventRecord.TimeStamp, sortIdx] = sort(EventRecord.TimeStamp);
        EventRecord.EventType          = EventRecord.EventType(sortIdx);
        EventRecord.TimeMsFromMidnight = EventRecord.TimeMsFromMidnight(sortIdx);
        EventRecord.TimeSource         = EventRecord.TimeSource(sortIdx);
        EventRecord.Details            = EventRecord.Details(sortIdx);
        EventRecord.EventNumber        = (1:numel(EventRecord.EventType))';
    end

    warning('NGL:INTAN_ExtractEvents:insertedPreIni', ...
        ['%d itiOn (0) event(s) had no immediately preceding preIni (8). ', ...
         '%d repaired by inserting a synthetic 8 thirty samples before each 0. ', ...
         '%d could NOT be repaired (the 30-sample window already held another ', ...
         'event, or fell before sample 1); those zeros remain orphaned and ', ...
         'trialdefGen will see them out of sequence. ', ...
         'Most likely cause: the behavioural task was not sending preIni reliably ', ...
         '(check the task script and/or DIO wiring).'], ...
        numel(orphanZeros), numel(insertTs), numel(skippedTs));
end

fprintf('Successfully created ''EventRecord'' structure.\n');
end