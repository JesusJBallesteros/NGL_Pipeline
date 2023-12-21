function [trl, def] = trialdefFromPar(opt)
% For now, testing in Pilot_SocialLearning: load '*par.mat'
% Loads the variable used to store experimental progress and control, which
% includes the eventcodes sent by Matlab to cameras and Deuteron, as
% decimal integers.
%
% Jesus. 21.12.2023

%% Read event definitions and Output as new variable
def = opt.def;

%% Locate and load data
% find 'par' file for session and variable 'SaveEvnts'
bhvfile = ls([opt.behavFiles, '\*par.mat']);
evnt    = load([opt.behavFiles, '\', bhvfile], 'SaveEvnts'); % collect events in temp variable 

% Collect eventcodes and timestamps from 'SaveEvnts'
events = evnt.SaveEvnts; 
clear evnt

%% Create a trial array (startTime, endTime, ZeroTime)
% Set all times relative to session start time
events(:,1) = events(:,1) - events(1,1);

%Find all timestaps for trial start event.
def.tstart = events(events(:,2)==def.itiOn);
def.tend   = events(events(:,2)==def.trialEnd);
    def.tend(end) = []; % rmv last 7 event (session end)

% Find how many start trial events exist
def.ntrials = length(def.tstart);

% Create a trial array 'trl' for FielTrip (nx3 where columns are 'trial start time',
% 'trial end time' and 'offset to zero').
trl      = zeros(def.ntrials,3);
trl(:,1) = def.tstart;
try trl(:,2) = def.tend;
catch, trl(:,2) = def.tstart+6;
end
% we leave offset to 0, so we can re-define where we want to zero time to
% be aligned flexibly in the future.

end