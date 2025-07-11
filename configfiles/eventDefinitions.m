function [eventdef] = eventDefinitions(format)
%% Definition of events for your own paradigm and recording system.
% Save it inside your project file system, under 'analisysCode'
% Jesus 24.08.2024

%% INFO. Evencodes are used to timestamp behavioral events within the ephys time series.
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
% of 0s and 1s. For now, a single 1 within a rest of 0s.
%
% This function defines the bit sequences reflecting this convention.
%
% Copyright:	Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
% Contact:      OTBR-Toolbox@ruhr-uni-bochum.de
% Source Code:  https://gitlab.ruhr-uni-bochum.de/ikn/OTBR
%
% Please cite the OTBR Toolbox where this function is used.
% This function is called in myHardwareSetup
%
% Author: Sara
% Version: 1.1
% Date: 30.08.2021
%
% 15.10.2021, Sara : - added 'reset' to standard events.
% 25.03.2022, Tobias: removed computations that are done in initOTBR
% 24.08.2023, Jesus: Complete modification to fit Deuteron single-bit
%                    concept. Removed output
% 21.02.2024, Jesus: Final event code convention for Deuteron.
% 27.02.2024, Juan: added treatments [tr 1&2]

%% INFO. DEUTERON using the exact name that will be provided in the task running script,
    % eventdef.eventName
    % and assign it an unique sequence of 0/1s. The idea is that the event
    % sequence within a trial only needs to change ONE pin at a time, to
    % code for any specific event.
    % After any end of trial code, the pins are needed to set back
    % stepwise to one step away from [0 0 0 0], so when it changes
    % to that, it stamps for the start of the trial.
    % What we need to read is which pin changed (stamped time) and what's
    % the current state of all the other, to get the new word.
    %
    % Example sequence:
    % ...
    % end1      [0 1 0 0] (trial-1 ends with 'omission' flag. Ready to send itiOn.)
    % itiOn     [0 0 0 0] (the trial starts, baseline.)
    % stimOn1   [0 0 0 1] (stimulus is presented, keybuffer opens.)
    % bhv       [0 0 1 1] (a response is registered within allowed time, keybuffer closes.)
    % pun       [1 0 1 1] (the response was incorrect, punishment feedback is sent.)
    % end3      [1 1 1 1] (trial ends with 'incorrect' flag. Not ready to send itiOn.)
    %           [1 1 1 0] (setting back to itiOn. Not ready to send itiOn.)
    %           [1 1 0 0] (setting back to itiOn. Not ready to send itiOn.)
    %           [1 0 0 0] (setting back to itiOn. Ready to send itiOn.)
    % itiOn     [0 0 0 0] (trial+1 starts, baseline.)
    %...
%% INFO. INTAN. Old way
% | decimal | binary    | NAME      | meaning                                           |
% | ------- | ------    |------     |---------------------------------------------------|
% | 0       | 0000      | reset     | set all pins to zero                              |
% | 1       | 0001      | itiOn     | start of the trial                                |
% | 2       | 0010      | stimOn    | any stimulus on: auditory/ visual/ neural ...)    |
% | 3       | 0011      | off       | anything off: stimulus/ reward/ punishment ...)   |
% | 4       | 0100      | bhv       | any behavior: peck/ fixation/ location ...)       |
% | 5       | 0101      | rwd       | reward on                                         |
% | 6       | 0110      | pun       | punishiment on                                    |
% | 7       | 0111      | end       | trial end                                         |

%% INFO. INTAN. New way
% Evnts include always the 16 pins, there is no stdEvents and extraEvents, all Events are now
%   espciefied in the defineEvntCodes function. As before the 4 first pins are reserved for fixed
% Evnts (dec 0:15), extra Evnts should be defined independently for each experiment (dec 16 +)
% We are using the names from deuteron to keep the analysis consistent.
% Order of pins in Intan and Deuteron is inverted, dec is always the same, how to read bin is changed
%   the binVec is now calculated as binVec = single(dec2binvec(eventdef.(stdEvents), SETUP.events.pinsLen));
% By default Intan starts with [1 1 0 0 ...] so YOU MUST END each experiment with this Evnt for consistency.


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
    
    % THIS example would only work for Juan's S3 project in box
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
