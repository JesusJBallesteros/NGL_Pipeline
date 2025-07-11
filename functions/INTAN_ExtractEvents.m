function EventRecord = INTAN_ExtractEvents(opt)
% Based on original function readEvents()
% Use this function to read event-codes saved in Intan (one file per channel).
%   Current version looks for any change in a digital pin, using that time
%   as event start. It controls for inconsistencies by considering 28
%   samples after this time (0.9ms). Unit of time is sample index (not
%   seconds) and it is relativized to the first event (start of session).
%   This version is NOT backwards compatible.
%
% INPUTS-OPTIONAL
%  * opt           : struct with options and parameter for the current run
%
% OUTPUTS
%  * EventRecord   : Vector containing all event-codes.

% VERSION HISTORY:
% Author:        Jonas Rose
% Version:       2.0
% Last Change:
% 08.05.2016, Jonas: Release version
% 17.05.2016, Jonas: sampling rate is picked up from header file or input
% 17.05.2016, Jonas: bugfix, read events as uint16
% 05.05.2022, Aylin: corrected code for new Intan System
% 25.05.2022, Aylin: corrected code for standard and extra event codes
% 11.10.2022, Jesus: Testing old format compatibility (Does not affect new format)
% 23.08.2024, Jesus: V 2.0 Reduced input/output to basics. General modification
%                   of digital pins reading to accomodate the unification of standard and
%                   extra events into just events, as we will use whole 16 bit words. This
%                   unifies the standard coding between Deuteron and Intan. Implies that the
%                   decimal integer 0 now has a meaning, and it is not just a reset.
% 24.04.2025, Jesus: Added output field .TimeBreak to match the detection
%                   from Deuteron System. Unlikely that they will happen on
%                   Intan Systems, so it will just be an empty 1x2 cell array.

%% Defaults
pth     = opt.PathRaw; % folder for reading events
smpDel  = 14*2;       % an event code is read after smpDel since first pin change, for additional smpDel since (to catch instabilities)

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
% remove samples during which pin 1 and 2 are up (per default before reset of all pins,
% necessary in order to correctly extract all relevant events).
startState = dIn(1,:); % Read pins states as recording starts 
pinsOff = find(any(dIn(:,1:npins)~=startState,2), 1, "first");  % find first pin change
if pinsOff ~= 1 % Only if is not already the first sample
    dIn(1:pinsOff, :) = []; % remove all samples until then
end
clear pinsOff

%% Find samples at which any pin changes
checksum = diff([int8(zeros(1,npins)); dIn],1,1); % fixed to admit negative values by changing uint8 to int8
ts = find(any(checksum,2));
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

% !!! ***
% checksum = single(sum(dIn,2));
% ts = find(diff([0; checksum])~=0);
% clear checksum 

%% CONVERT all events
% convert each binary word to its corresponding decimal using the npins bits
EventType = nan(size(ts,1),1);

for i = 1:size(ts)
    % Convert binary pins to decimal, as sum over smpDel forward to catch inconsitencies
    EventType(i) = binvec2dec(sum(dIn(ts(i):ts(i)+smpDel,:))); % binary vector to decimal integer
end
clear dIn

%% Place extracted information into a proper EventRecord
EventRecord.EventType           = double(EventType);
EventRecord.EventNumber         = double(1:1:length(EventType))';
EventRecord.TimeStamp           = nan(length(EventType),1);
EventRecord.TimeMsFromMidnight  = ts/(input.sessions.info.amplifier_sample_rate/1000);
EventRecord.TimeSource          = nan(length(EventType),1);
EventRecord.Details             = nan(length(EventType),1);
EventRecord.TimeBreak           = {[] []};
end