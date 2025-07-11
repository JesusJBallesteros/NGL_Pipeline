function conditions = get_trialConditions(opt, events)
% Funtion for indexing trials based on a single characteristic
% matching an specific event. Useful to extract trials of interest,
% grouped by e.g. stimulus presented, context, animal response, etc.
% Later analyses can be applied directly to each of the indexes only,
% excluding the trials of no interest.
%
% Jesus 29.08.2024

    % Initialize all fields with zeros
    zerovec = zeros(length(events.itiOn.code),1);
    % Basic ones
    conditions = struct( 'aborted', zerovec, 'correct', zerovec, 'incorrect', zerovec, ...
                         'omission2', zerovec, 'responsive', zerovec); 
    
    % Trial by trial, classify them. Using 'itiOn' bc it should be the minimal
    % choice of alignment, and events are the same for any alignment.
    for i=1:length(events.itiOn.code)
        trialvect = events.itiOn.code{i,1};

        % Get if the trial was not initiated
        if any(ismember(trialvect, [opt.eventdef.oms1]))
            conditions.aborted(i) = 1;
        end
        % Get if the trial was correct
        if sum(ismember(trialvect, [opt.eventdef.bhv opt.eventdef.rwd]))==2
            conditions.responsive(i) = 1;
            conditions.correct(i) = 1;
        end
        % Get if the trial was incorrect
        if sum(ismember(trialvect, [opt.eventdef.bhv opt.eventdef.pun]))==2
            conditions.responsive(i) = 1;
            conditions.incorrect(i) = 1;
        end
        % Get if the stimulus2 was omitted
        if any(ismember(trialvect, [opt.eventdef.oms2]))
            conditions.omission2(i) = 1;
        end
        
%% Examples and extensions
%         % example, to obtain wich context was used at trial (by code)
%         if any(ismember(trialvect, [opt.eventdef.CtxtA opt.eventdef.CtxtA2 ...
%                                     opt.eventdef.CtxtB ...
%                                     opt.eventdef.CtxtC opt.eventdef.CtxtC2]))
%             condition.ctxt(i) = trialvect(ismember(trialvect, [opt.eventdef.CtxtA opt.eventdef.CtxtA2 ...
%                                     opt.eventdef.CtxtB ...
%                                     opt.eventdef.CtxtC opt.eventdef.CtxtC2]));
%         end
%
%         % example, wich stimulus was shown at trial (by code)
%         if any(ismember(trialvect, [opt.eventdef.FS opt.eventdef.NS1 opt.eventdef.NS2]))
%             condition.ctxt(i) = trialvect(ismember(trialvect, [opt.eventdef.FS opt.eventdef.NS1 opt.eventdef.NS2]));
%         end
%
%         % example, to know the response made by the animal, if an event code 
%         % is sent for THAT specific answer (ie, peak at location 4, sendevent X;
%         % but if peaks at location 9, sendevent Y).
%         if any(ismember(trialvect, [opt.eventdef.EVENTX/Y/Z]))
%             condition.response(i) = 4/5/6;
%         end
%
%         % FOR MORE add as (with corresponding logical index): 
%         % if ismember(trialvect, [opt.eventdef.xxx opt.eventdef.yyy])
%         %     condition.xxxyyy(i) = 1;
%         % end
%         % if ismember(trialvect, [opt.eventdef.xxx opt.eventdef.yyy])
%         %     condition.xxxyyy(i) = 1;
%         % end
    end
end