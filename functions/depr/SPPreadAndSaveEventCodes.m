function readAndSaveEventCodes(pathRawData,savePath,subject,session,eventLib)
%%
%
% Use this function to read out the events of a session from the raw
% INTAN-files and save them in MATLAB-Files in raw and trial-by-trial
% format. Event codes are saved in individual folders within
% 'data/sorted/...' per subject and named according to the corresponding 
% session. 
%
%INPUTS
%  * 'pathRawData'      : string of path to raw data (INTAN-files)
%  * 'savePath'         : string of path for saving
%  * 'subject'          : strings of subject
%  * 'session'          : string of session date (e.g. '220528')
%  * 'eventLib'         : event-code library containing fields
%        .start         : start code of all trials (= itiOn)
%        .end           : end code of all trials (= end)

%OUTPUTS
%  * 'events'           : struct labelled 'events###_YYMMDD' (### = bird)
%        .trialEvents   : matrix containing smpInd (index of sample 
%                         corresponding to each event-code), standardEvent, 
%                         and extraEvent (0 if not defined) sorted into
%                         trials 
%        .rawEvents     : matrix containing smpInd, standard event code, 
%                         extra event code of entire recording

% VERSION HISTORY:
% Author:         Lukas Hahn + Aylin Apostel
% Version:        2.0.0
% Last Change:    30.05.2022
%
% 02.08.2019, Lukas: v1.0.0 release version
% 25.05.2022, Aylin: updated code for new sendEvent function 
% 30.05.2022, Aylin: updated code for correct folder structure

% BUGS/ TODO:
% 

%% UPDATE - MAJOR CHANGES
% - only one session of one subject

%%% Select raw data files: 
% even for recording with multiple electrodes only one DigitalIn file per pin
eventPath = fullfile(pathRawData,subject,session);

%%% Read event codes out of respective intan files:
[ev, ~, ~] = readEvents('path',eventPath,'smprate',30000,'pins',2:15);
rawEvents  = [ev.smpInd ev.eventStandard ev.eventExtra]; % [sample index, standard event, extra event]

% Sort event codes into trials:
spp = 1;
if spp==1 %only for SPP
    trialStarts = [1:size(rawEvents,1)-1]';
    trialEnds = [[1:size(rawEvents,1)-1]+1]';
else
    trialStarts = find(rawEvents(:,2) == eventLib.start); % itiOn (standard event)
    trialEnds   = find(rawEvents(:,2) == eventLib.end);   % end (standard event)
end
trialEvents = cell(size(trialStarts,1),1);
for j = 1:size(trialStarts,1)
    trialEvents{j,1} = rawEvents(trialStarts(j,1):trialEnds(j,1),:);
end
events.trialEvents = trialEvents;
events.rawEvents   = rawEvents;
fileName           = ['events' subject '_' session];

% Create new directory to save struct events
mkdir(fullfile(savePath,subject),session);          % create new directory for current session
savDir = [fullfile(savePath,subject,session),'\'];	% where to save the data
save([savDir,fileName],'events');

end