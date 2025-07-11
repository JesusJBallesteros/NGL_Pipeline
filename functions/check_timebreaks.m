function [timebreak] = check_timebreaks(data)
% By analyzing the trial definition timings or a event log from Deuteron, 
% we can check for consistency or continuity in the recorded timing 
% across the session. 
% For trialdef, regular ITIs should stay close to a median value,
% while deviations could occur at block or phase changes. If so, those
% would probably happen multiple times, again in a consistent fashion. If
% yet any other deviation happens, it could mean a break in the recording
% session, i.e. a a time where perhaps there is not actual data
% acquisition, for a number of reasons. In that case, it would be at least
% advisable to check if the spike timing has become de-synchronised due to
% this time break.
% For Event records, the detection of a 'Stopped recording' should suffice
% to locate those when unespected, i.e. not at the end.
% INPUT: a Nx3 numerical array (trialdef)
%        or
%        a Nx10 string array (myRecord)
% OUTPUT: 
% 'timebreak': cell array with {1} the trial number after the break
%                          and {2} the trial definition for that trial
% Jesus 24.04.2025

%% A variable will determine if there is need for correction. Create empty.
timebreak = {[] []};

%% Check what INPUT type we have: trialdef (1) or Deuteron eventlog (2)
if isstring(data), datatype = 2;
else,              datatype = 1; 
end

%% Proceed with case
switch datatype
    case 1
        % Get the time difference between every itiOn
        trialtimes = diff(data(:,3));
        mediantrialtime = median(trialtimes,"omitnan");
        
        longtrialtimes = trialtimes(trialtimes > mediantrialtime*1.2);
        nlong = numel(longtrialtimes);
        if nlong > 1
            medianlongtrialtimes = median(longtrialtimes,"omitnan");
            toolongtrialtimes = longtrialtimes(longtrialtimes > medianlongtrialtimes*1.1);
            nlong = numel(toolongtrialtimes);
            if nlong > 0
                trial = find(trialtimes==toolongtrialtimes)+1;      
                timebreak = {trial, data(trial,:)-data(trial-1,:)};
            end
        end

    case 2
        % Find those logs with Stopped recording info.
        breakDect = contains(data(1:end-1,6), 'Stopped recording'); % Actual end of recording not accountable
        breakDect = find(breakDect==1);
        if ~isempty(breakDect)
            breakDect = [breakDect breakDect+1]; % Get that stamp and the following, restart
            % There should not be more than one, but just in case
            for b = 1:size(breakDect,1)
                % Get the times for the last and first, to obtain real time
                % passed in between them.
                timebreak{b,1} = []; %
                timebreak{b,2} = [str2double(data(breakDect(b,1), 3)) str2double(data(breakDect(b,2), 3))]; % Keep 3rd field
            end
        end
end