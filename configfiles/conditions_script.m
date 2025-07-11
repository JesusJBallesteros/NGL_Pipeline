%% Personalize this script to extract conditions relevant or interesting to your Project
% Commonly, the top ones (Response, correct, incorrect, omissions) could be generalized 
% to any task, meanwhile the folloeing ones would be rather specific.
% Save this SCRIPT under your 'analysisCode' folder.
% Note that this is NOT a FUNCTION.

%% Create the conditions variable for trial indexing
% Initialize all fields with zeros
zerovec = zeros(length(events.itiOn.code),1);
conditions = struct( 'correct', zerovec, 'incorrect', zerovec, 'omission', zerovec, ...
                    'response', zerovec, 'aborted', zerovec, 'stimulus', zerovec);

% Trial by trial, classify them. Using 'itiOn' bc it should be the minimal
% choice of alignment, and events are the same for any alignment.
for i=1:length(events.itiOn.code)
    trialvect = events.itiOn.code{i,1};

    % Response / Correct
    if sum(ismember(trialvect, [opt.eventdef.bhv opt.eventdef.rwd]))==2
        conditions.response(i) = 1;
        conditions.correct(i) = 1;
    end

    % Response / Incorrect
    if sum(ismember(trialvect, [opt.eventdef.bhv opt.eventdef.pun]))==2
        conditions.response(i) = 1;
        conditions.incorrect(i) = 1;
    end

    % Omissions
    if any(ismember(trialvect, [opt.eventdef.oms1 opt.eventdef.oms2]))
        conditions.omission(i) = 1;
    end

    % Stimulus used 
    % 'tr2' defines the Novel Stimuli in Extinction_Arena project
    if ismember(trialvect, [opt.eventdef.tr2])
        conditions.stimulus(i) = 1;
    end

    % FOR MORE add as (with corresponding logical index): 
    % if ismember(trialvect, [opt.eventdef.xxx opt.eventdef.yyy])
    %     condition.xxxyyy(i) = 1;
    % end
    % if ismember(trialvect, [opt.eventdef.xxx opt.eventdef.yyy])
    %     condition.xxxyyy(i) = 1;
    % end
end