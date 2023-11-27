function [trl, def] = trialdefFromPar(opt)
% For now, testing in Pilot_SocialLearning: load '*par.mat'

%% Get event definitions
def = opt.def;

%% Locate and load data
% find par file for session and variable 'SaveEvnts'
bhvfile = ls([opt.behavFiles, '\*par.mat']);
evnt = load([opt.behavFiles, '\', bhvfile], 'SaveEvnts'); % collect events in temp variable 

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

% Find how many trial strat events are there
def.ntrials = length(def.tstart);

% Create a trial array 'trl' for FT (nx3 where columns are 'trial start time',
% 'trial end time' and 'offset to zero').
trl      = zeros(def.ntrials,3);
trl(:,1) = def.tstart;
trl(:,2) = def.tend;
% we leave offset to 0, so we can re-define where we want to zero time to
% be aligned flexibly in the future.

end