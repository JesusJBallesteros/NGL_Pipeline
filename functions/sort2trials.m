function [neurons, neurons_FT] = sort2trials(spike, trialdef, opt)
% This function takes the spike data and the trial definition to sort spike
% timestamps into the different trials where they belong, and reliativizes
% this time to the current trial alignmen. As a special case, a non-cell
% 'trialdef' input triggers the spike timestamps relative to the
% tutor-tutee interaction times in a Social paradigms.

neurons = [];
neurons_FT = struct();

nclus = length(spike.label);

% For each cluster
for c = 1:nclus
    if iscell(trialdef)
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

    else
        % Spiking indexing for Social interactions. Checks blob interaction 
        % times and extract spiking activity around them.
        ntrial = size(trialdef,1);
        neurons{c,1} = cell(ntrial,1);
        % for each trial
        for i=1:ntrial
            st = spike.timestamp{1,c}*1000; % convert spike times to msec

            % index for spiketimes at each interaction, including 5 seconds
            % before and after it happens
            idx = st >= (trialdef(i,1)-5)*1000 & ... % btw interaction start
                  st <  (trialdef(i,2)+5)*1000;      % and interaction end
            
            % relativize times to the given alignment point
            neurons{c,1}{i,1} = st(idx) - trialdef(i,1)*1000;
        end

        %% IN DEVELOPMENT
%         % In addition, find event times and and codes at video-asessment
%         % file, for social interaction cues extracted by students at Juan's
%         % Social paradigm
%         for i=1:size(opt.alignto,1)
%             events.(opt.alignto{i,1}) = [];
%         
%             for t = 1:ntrial
%                 % Grab all timestamps between time of start and time of end (inclusive)
%                 trialstamps = EventRecord.TimeSecFromMidnight(EventRecord.TimeSecFromMidnight >= trialdef{2,i}(t,1)/1000 & ...
%                                                              EventRecord.TimeSecFromMidnight <= trialdef{2,i}(t,2)/1000);
%                 % Relativize trial timestamps to alignment offset
%                 trialstamps = trialstamps - trialdef{2,i}(t,3)/1000; 
%             
%                 % Grab all events ocurring between time of start and time of end (inclusive)
%                 trialevents = EventRecord.EventType(EventRecord.TimeSecFromMidnight >= trialdef{2,i}(t,1)/1000 & ...
%                                                     EventRecord.TimeSecFromMidnight <= trialdef{2,i}(t,2)/1000);
%         
%                 % Insert into the proper structure to be output.
%                 events.(opt.alignto{i,1}).code{t,1} = trialevents; 
%                 events.(opt.alignto{i,1}).time{t,1} = trialstamps; 
%             end
%         end

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