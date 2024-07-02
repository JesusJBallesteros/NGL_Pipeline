function [neurons, neurons_FT] = sort2trials(spike, trialdef, opt)
% This function takes the spike data and the trial definition to sort spike
% timestamps into the different trials where they belong, and reliativizes
% this time to the current trial alignment 
nclus = length(spike.label);
neurons = struct();
neurons_FT = struct();

% For each cluster
for c = 1:nclus 
    % for each requested alignment
    for a = 1:size(opt.alignto,2)
        ntrial = size(trialdef{2,a},1);
        neurons.(opt.alignto{1,a}){c,1} = cell(ntrial,1);
        % for each trial
        for i=1:ntrial 
            st = spike.timestamp{1,c}*1000; % convert spike times to msec
            
            % index for spiketimes ...
            idx = st >= trialdef{2,a}(i,1) & ... % btw trial start
                  st <  trialdef{2,a}(i,2);      % and trial end
            
            % relativize times to the given alignment point
            neurons.(opt.alignto{1,a}){c,1}{i,1} = st(idx) - trialdef{2,a}(i,3); 
        end
    end
end

%% TODO make FT spike trial parsing work
% %% With purpouse of putting them togehter with the LFP for i.e.
% % spike-field analisys. Do as many rounds as existing trial alignments
% % for each requested alignment
% for a = 1:size(opt.alignto,2)
%     cfg = [];
%     cfg.trl = trialdef{2,a}; % Trial times come in msec
%     cfg.timestampspersecond = 3200; % timestamps per second to synchr spike timestamps to trial samples.
%     % For some reason this has to be divided? (vs what I expected, 32000, as sampled from Deuteron)
% 
%     % Transform manual 'spike' into most basic FT structure. 
%     % Keep only fields recognized by FT. The output spike structure usually contains
%     %   spike.label     = 1xNchans cell-array, with channel labels
%     %   spike.waveform  = 1xNchans cell-array, each element contains a matrix (Nleads x Nsamples X Nspikes)
%     %   spike.waveformdimord = '{chan}_lead_time_spike'
%     %   spike.timestamp = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
%     %   spike.unit      = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
% 
%     spikeFT.label     = spike.label;
%     spikeFT.timestamp = cellfun(@transpose, spike.timestamp, 'UniformOutput', false);
%     neurons_FT.(opt.alignto{1,a}) = ft_spike_maketrials(cfg, spikeFT);
% end
% 
% % Save output
% save(fullfile(opt.analysis, "neurons_FT.mat"), 'neurons_FT', '-mat')

end