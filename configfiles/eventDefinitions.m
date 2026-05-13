function [eventdef] = eventDefinitions(format)
% eventDefinitions  Define event code for general/specific projects.
%
% PURPOSE:
%   Returns the complete event definition struct for the current recording
%   format, combining the fixed reserved events (codes 0–15, hardware-locked)
%   with project-specific events (codes ≥ 16, user-configurable per format).
%   Called by EventProcess on first run of a session.
%
%   The function in the toolbox is meant to be copied and modified as needed,
%   then saved under analysisCode\
%
% USAGE:
%   eventdef = eventDefinitions(format)
%   where format is one of: 'fileperch', 'filepertype', 'DT2', 'DF1'
%
% INPUT:
%   format  - (char) recording format string as in info.fileformat
%
% OUTPUT:
%   eventdef - struct; field names are event names (char), values are decimal
%              integer codes. Reserved codes 0–15 are blocked.
%              Project-specific codes (≥ 16) are added per format block.
%
% CUSTOMISATION:
%   Copy this file from configfiles\ to your project's analysisCode\ folder.
%   Edit ONLY the two blocks labelled "ONLY MODIFY THIS TWO BLOCKS":
%     - Add INTAN project events under: elseif strcmpi(format,'fileperch')
%     - Add Deuteron project events under: if strcmpi(format,'DF1')
%   DO NOT modify reservedEvents() or any code above the modify blocks.
%   DO NOT reuse codes 0–15 for project-specific events.
%   Use decimal integers starting from 16 (included); each must be unique.
%
% RESERVED EVENTS (DO NOT CHANGE):
%   itiOn=0, stimOn1=1, stimOn2=2, bhv=3, end1=4, oms1=5, oms2=6, rwd=7,
%   preIni=8, tr1=9, end2=10, pun=11, na1=12, tr2=13, na2=14, end3=15
%
% Jesus 06.05.2026

%% INFO GENERAL. Evencodes are used to timestamp behavioral events within the ephys time series.
% This means, behavioral events (stimulus presentation, peaks) are
% generated or expected by Matlab, and either presented or captured by
% external hardware. These events NEED to be precisely represented as a
% time point, to which we can align the neural data afterwards.
%
% These events, to be aligned with neural data, need to me timely defined:
% If we are presenting a stimulus, we need to timestamp the exact moment when
% the stimulus is PRESENTED to the animal in the display, not when it is
% generated in matlab or sent from matlab to the display (there may be delays
% between the ORDER to display and the ACTUAL display). The same way, a
% peak response need to be represented as precise as possible to the ACTUAL
% peak of the animal.
%
% We follow the convention where each event is coded by a single sequence
% of 0s and 1s. This function defines the bit sequences reflecting this convention.

%% INFO. DEUTERON using the exact name that will be provided in the task running script,
% eventdef.eventName
% and assign it an unique sequence of 0/1s. The idea is that the event
% sequence within a trial only needs to change ONE pin at a time, to
% code for any specific event.
% After any end of trial code, the pins are needed to set back
% stepwise to one step away from [0 0 0 0], so when it changes
% to it stamps for the start of the trial.
% What we need to read is which pin changed (stamped time) and what was
% the previous state of all the other, to get the new word.
%
% Example sequence:
% ...
% end1      [0 1 0 0] (trial N-1 ends with 'omission' flag. Ready to send itiOn, so no transition events exist.)
% itiOn     [0 0 0 0] (trial N starts, baseline.)
% stimOn1   [0 0 0 1] (stimulus 1 is presented, keybuffer opens.)
% bhv       [0 0 1 1] (a response is registered within allowed time, keybuffer closes.)
% pun       [1 0 1 1] (the response was incorrect, punishment feedback is sent.)
% end3      [1 1 1 1] (trial ends with 'incorrect' flag. Not ready to send itiOn.)
% na1       [1 1 1 0] (setting back. Not ready to send itiOn.)
% na2       [1 1 0 0] (setting back. Not ready to send itiOn.)
% na3       [1 0 0 0] (setting back. Ready to send itiOn.)
% itiOn     [0 0 0 0] (trial N+1 starts, baseline.)
%...

%% INFO. INTAN.
% Evnts include always the 16 pins. As before, the 4 first pins are reserved for fixed
% Evnts (decimals 0-15), and "extra" events should be defined independently for each experiment.
%
% We are using the names from deuteron to keep the analysis consistent.
%
% Order of pins in Intan and Deuteron are inverted: 
%   decimal is always the same, how to read bin is changed.
%   the binVec is now calculated as binVec = single(dec2binvec(eventdef.(stdEvents), SETUP.events.pinsLen));
% By default Intan starts with [1 1 0 0 ...] so YOU MUST END each experiment with this event for consistency.

%% DO NOT MODIFY. Function
eventdef = reservedEvents(); % hic sunt dracones. DO NOT MODIFY!

% THis is a DESCRIPTION of the events. DO NOT uncomment/change anything without explicit consent.

% % % WITHIN TRIAL
% % % itiOn   = 0  [0 0 0 0]  Trial start.
% % % stimOn1 = 1  [0 0 0 1]  Stim1 presentation (INI, Sample, etc).
% % % stimOn2 = 2  [0 0 1 0]  Stim2 presentation (Match, choice, cue, etc) .
% % % oms1    = 5  [0 1 0 1]  Omission to Stim1.
% % % oms2    = 6  [0 1 1 0]  Omission to Stim2.
% % % bhv     = 3  [0 0 1 1]  A response, or behaviour of intertest, is detected.
% % % rwd     = 7  [0 1 1 1]  A reward is given.
% % % pun     = 11 [1 0 1 1]  A punishment is presented
% % % end1    = 4  [0 1 0 0]  End of trial after any omission.
% % % end2    = 10 [1 0 1 0]  End of trial after punishment.
% % % end3    = 15 [1 1 1 1]  End of trial after reward.

% % % OUT TRIAL
% % % tr1     = 9  [1 0 0 1]  Treatment/block/phase 1. Or odd blocks/phases/... Or block/phase/treatment start.
% % % tr2     = 13 [1 1 0 1]  Treatment/block/phase 2. Or even blocks/phases/... Or block/phase/treatment end.
% % % na2     = 14 [1 1 1 0]  Transition sequence. (meaningless to Intan).
% % % na1     = 12 [1 1 0 0]  Transition sequence. (meaningless to Intan).
% % % preIni  = 8  [1 0 0 0]  Transition sequence. (meaningless to Intan).

%% ONLY MODIFY THIS TWO BLOCKS.
if strcmpi(format,'DF1') 
    % Deuteron should accept specific events in a near future, stay tuned.
    % % Future use
    % % Future use
    % % Future use

elseif strcmpi(format,'fileperch')
    % PROJECT-SPECIFIC EVENTS. INTAN allows for it. 
    % With 16 bits you can represent 65536 different values, when taken 
    % as a group. It is convention to use 16 bits to represent the integers
    % from 0–65535. Since 0-15 are RESERVED, you can ONLY code additional 
    % 65519 decimal integers, from 16 to 65535. Use them wisely.
    
    % This EXAMPLE would work for an specific project
    % Descriptions      = Decimal;  % [binary];           % Comments
    eventdef.CtxtA      = 16;       % [0000100000000000]; %
    eventdef.CtxtB      = 17;       % [1000100000000000]; %
    eventdef.CtxtA2     = 18;       % [0100100000000000]; %
    eventdef.CtxtC      = 19;       % [1100100000000000]; %
    eventdef.CtxtC2     = 20;       % [0010100000000000]; %
    eventdef.CtxtD      = 21;       % [1010100000000000]; %
    eventdef.FS         = 22;       % [0110100000000000]; %
    eventdef.NS1        = 23;       % [1110100000000000]; %
    eventdef.NS2        = 24;       % [0001100000000000]; %
    eventdef.bhvChoice  = 25;       % [1001100000000000]; %
    eventdef.ACQphase   = 26;       % [0101100000000000]; %
    eventdef.EXTphase   = 27;       % [1101100000000000]; %
    eventdef.TESTphase  = 28;       % [0011100000000000]; %
end

end

%% hic sunt dracones. ALL RESERVED. DO NOT MODIFY
function [eventdef] = reservedEvents()
% Event struct        = Decimal;
eventdef.itiOn        = 0;
eventdef.stimOn1      = 1;
eventdef.stimOn2      = 2;
eventdef.bhv          = 3;
eventdef.end1         = 4;
eventdef.oms1         = 5;
eventdef.oms2         = 6;
eventdef.rwd          = 7;
eventdef.preIni       = 8;
eventdef.tr1          = 9;
eventdef.end2         = 10;
eventdef.pun          = 11;
eventdef.na1          = 12;
eventdef.tr2          = 13;
eventdef.na2          = 14;
eventdef.end3         = 15;
end
