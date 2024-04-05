function [EventRecord] = trialdefFromPar(opt)
% For now, as a safeguard againt failures in Deuteron reading 
% events from behavior: load '*res.mat' and read SaveEvnts variable.
% Loads the variable used to store experimental progress and control, which
% includes the eventcodes sent by Matlab to cameras and Deuteron, as
% decimal integers.

% Jesus. 26.03.2024

%% Read event definitions into new variable.
def = opt.eventdef;

%% Locate and load events.
% Find 'par' file for session and variable 'SaveEvnts'
bhvfile = ls('*res.mat');
evnt    = load(bhvfile, 'SaveEvnts'); % collect events in temp variable 

% Collect eventcodes and timestamps from 'SaveEvnts'
events = evnt.SaveEvnts; 
clear evnt

events = events(~isnan(events));
events = reshape(events,[],2);

%% Create a trial array (startTime, endTime, ZeroTime)
% Set all times relative to session start time
events(:,1) = events(:,1) - events(1,1);

% Find all timestaps for trial start and end events (trial boundries).
def.tstart = events(events(:,2) == def.itiOn); % Could be rather itiOn-t
def.tend   = events(events(:,2) == def.end1 | events(:,2) == def.end2 | events(:,2) == def.end3);
%     % rmv last 7 event (session end)
%     def.tend(end) = []; % Bc is duplicated by design. To be elimininated

% Find how many start trial events exist.
def.ntrials = length(def.tstart);

% Create a trial array 'trl' for FielTrip (Nx3, where columns are:
% 'trial start time', 'trial end time' and 'offset to zero').
trl      = zeros(def.ntrials,3);
trl(:,1) = def.tstart;
try trl(:,2) = def.tend; % tend exists
catch, trl(:,2) = def.tstart+6; % no tend found
end
% we leave offset to 0, so we can re-define where we want to zero time to
% be aligned flexibly in the future. Would rather be tstart+t

%% Place extracted information into a proper EventRecord
EventRecord.EventNumber = double(1:1:length(events))';
EventRecord.EventType = single(events(:,2));
EventRecord.TimeStamp = string(events(:,1)); % Convert to string array
EventRecord.TimeMsFromMidnight = round(events(:,1)*1000);
EventRecord.TimeSource = nan(length(events),1);
EventRecord.Details = nan(length(events),1);

end