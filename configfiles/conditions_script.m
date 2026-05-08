%% Personalize this script to extract relevant conditions to your Project
%
% PURPOSE
%  To create a set of 'conditions' based on events found per trial, for an
%  easy indexing of behavioral or experimental classification. E.g. a
%  condition could be 'correct', therefore indexing all those trials ending
%  with a reward flag. Afterwards, analysis could be directed to only
%  'correct' condition, automatically limiting it to those trials, already
%  indexed. Other cases could be 'Stimulus A', 'Stimulus A|B', 
%  'Correct Stimulus A'...
%
% USAGE
%  Save this script under your 'analysisCode' folder. It runs as script (not 
%  function) called at EventProcess(), once events and trialdef have been
%  generated. 
%  1. MODIFY the 'conditions' struct construction below, adding "'condition', zerovec" pairs as needed.
%  2. CREATE the necessary logical indexing in the TRIAL LOOP: 
%     Example: ismember(trialvect, [opt.eventdef.tr2]) -> if tr2 appears in this trial, classify as...
% 
%  Conditions could be determine positively (is C) or negatively (not C),
%  so later on an indexig would be 'correct' vs '~correct' but it is 
%  recommended to stay in the positive side and explicit them. In this
%  example the indexing would result in 'correct' vs 'incorrect'.
%
% BASIC CONDITIONS
%  Commonly used ones (response, correct, incorrect, omissions...) could be generalized 
%  to any task, meanwhile others (Novel Stimuli) would be rather specific.
%  While the event name exists, even if the code never appears, the condition would be just empty. 

%% Create the conditions variable for trial indexing
% Initialize all fields with zeros
zerovec = zeros(length(events.itiOn.code),1);
conditions = struct( 'correct', zerovec, 'incorrect', zerovec, 'omission', zerovec, ...
                    'response', zerovec, 'aborted', zerovec, 'stimulus', zerovec); % Add "'condition', zerovec" pairs as needed

%% TRIAL LOOP
% Create and describe the logical indexing that would allow for each trial
% to be clasiffied as 'condition' or not. There cannot be ambiguities, each
% trial is C or NOT C. Each trial can be part of many conditions.
for i=1:length(events.itiOn.code)
    trialvect = events.itiOn.code{i,1};

    % Response & Correct
    if sum(ismember(trialvect, [opt.eventdef.bhv opt.eventdef.rwd]))==2
        conditions.response(i) = 1;
        conditions.correct(i) = 1;
    end

    % Response & Incorrect
    if sum(ismember(trialvect, [opt.eventdef.bhv opt.eventdef.pun]))==2
        conditions.response(i) = 1;
        conditions.incorrect(i) = 1;
    end

    % Omissions
    if any(ismember(trialvect, [opt.eventdef.oms1 opt.eventdef.oms2]))
        conditions.omission(i) = 1;
    end

    % Novel Stimuli
    % In this example, 'tr2' defines the Novel Stimuli in Extinction_Arena
    % but this could be your own stimuli code.
    % 
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