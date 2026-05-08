function [timebreak] = check_timebreaks(data)
% check_timebreaks  Detect recording breaks from trial definitions or Deuteron event logs.
%
% PURPOSE:
%   Analyses recording continuity using one of two input types:
%     trialdef (numeric) — checks inter-trial intervals against the session
%       median; statistical outliers that stand out above all regular long ITIs
%       are flagged as potential recording breaks.
%     Deuteron event log (string matrix) — scans column 6 for 'Stopped recording'
%       entries that appear before the final record; each such entry indicates
%       a gap during which acquisition was paused.
%
% USAGE:
%   timebreak = check_timebreaks(data)
%   Called from Deuteron_ExtractEvents (event-log path) and trialdefGen (trialdef path).
%
% INPUT:
%   data  - [N × 3 double]  trial definition array [start end t0] in ms
%           OR
%           [N × 10 string] Deuteron event log CSV matrix (column 6 = Details)
%
% OUTPUT:
%   timebreak  - {1 × 2} cell array:
%                  {1} trial index (or event row) immediately after the break
%                  {2} trial-definition row (or timing info) at that point
%                Returns {[] []} when no break is detected.
%
% Last modified 08.05.2026 (Jesus)

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
                timebr