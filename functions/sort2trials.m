function [neurons, neurons_FT] = sort2trials(spike, trialdef, opt)
% This function takes the spike data and the trial definition to sort spike
% timestamps into the different trials where they belong, and reliativizes
% this time to the current trial alignment 

% Do as many rounds as existing trial alignments
nclus = length(spike.label);
neurons = struct();
neurons_FT = struct();

% For each cluster
for c = 1:nclus 
    % for each requested alignment
    for a = 1:size(trialdef,2)
        ntrial = size(trialdef{2,a},1);
        neurons.(opt.alignto{1,a}){c,1} = cell(ntrial,1);
        % for each trial
        for i=1:ntrial 
            st = spike.timestamp{1,c}*1000; % seconds to msec
            
            idx = st >= trialdef{2,a}(i,1) & ... % index for spiketimes btw trial strat
                  st <  trialdef{2,a}(i,2);      % and trial end
            neurons.(opt.alignto{1,a}){c,1}{i,1} = st(idx) - trialdef{2,a}(i,3); % and relativize time to the given alignment point
        end
    end
end

%% TODO make FT spike trial parsing work
% %% With purpouse of putting them togehter with the LFP for i.e.
% % spike-field analisys. Do as many rounds as existing trial alignments
% % for each requested alignment
% for a=1:size(trialdef,2)
%     cfg = [];
%     cfg.trl = trialdef{2,a}; % Trial times come in msec
%     cfg.timestampspersecond = 3200; % timestamps per second to synchr spike timestamps to trial samples.
%     % For some reason this has to be divided? (vs what I expected, 32000, as sampled from Deuteron)
% 
%     % Transform manual 'spike' into a FT ready structure. 
%     neurons_FT.(opt.alignto{1,a}) = ft_spike_maketrials(cfg, spike);
% end
% 
% % Save output
% save(fullfile(opt.analysis, "neurons_FT.mat"), 'neurons_FT', '-mat')
% 
end