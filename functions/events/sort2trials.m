function [neurons, neurons_FT] = sort2trials(spike, trialdef, opt)
% sort2trials  Bin per-cluster spike timestamps into per-trial cell arrays,
%              relativised to a trial-anchor alignment event.
%
% PURPOSE:
%   Walks every (cluster, trial) pair and keeps only the spikes whose
%   timestamp falls inside the trial window defined by trialdef. Output
%   times are expressed in ms relative to the alignment event of that
%   trial. Two operating modes are selected by the type of 'trialdef':
%
%     1) STANDARD (trialdef is a CELL)
%        Used after trialdefGen. trialdef is a (2 x nAlignments) cell with
%        rows {1,a}=alignment name (char), {2,a}=[nTrials x 3] start/end/t0
%        in ms. Produces a struct keyed by alignment name:
%            neurons.<align>{c,1}{i,1} = vector of ms timestamps for
%                                       cluster c, trial i, relative to
%                                       trialdef{2,a}(i,3).
%        Also fills neurons.ROI(c) = spike.roi(c).
%
%     2) SOCIAL  (trialdef is NUMERIC, e.g. blob.Merges)
%        Used for social-arena interaction indexing. trialdef is a
%        [nInteractions x >=2] numeric matrix of [tStart tEnd] in SECONDS.
%        Each interaction is expanded by HARDCODED ±5 s padding; spikes
%        falling inside the padded window are kept and relativised to
%        tStart (ms). Output is a flat cell-of-cells:
%            neurons{c,1}{i,1} = vector of ms timestamps for cluster c,
%                                interaction i.
%        Note the shape difference vs mode 1 (no alignment keying, no
%        .ROI field).
%
% USAGE:
%   [neurons, neurons_FT] = sort2trials(spike, trialdef, opt)
%
% INPUTS:
%   spike    - NGL spike struct from loadSpikes. Required: .label,
%              .timestamp (cell-per-cluster, SECONDS), .roi (for mode 1).
%   trialdef - cell (mode 1) or numeric matrix (mode 2); see above.
%   opt      - options struct. Required (mode 1): opt.alignto (cell of
%              alignment names). Mode 2 does not use opt.alignto today.
%
% OUTPUTS:
%   neurons    - struct (mode 1) or cell (mode 2); shape per the section
%                above. Spike times in MILLISECONDS, relative to alignment.
%   neurons_FT - reserved for a FieldTrip-style spike struct. Currently
%                always returned as an empty struct(); the implementation
%                lives commented at the bottom of this file (see audit
%                item R / "TODO fix Fieldtrip extraction" in NGL02).
%
% KNOWN ISSUES:
%   - The ±5 s window in mode 2 is a magic number. Should become opt or
%     param-driven (audit item R).
%   - Mode 2 ("social") is project-specific and is wedged into a
%     general-purpose function. The intended refactor (audit item R)
%     splits it into a separate sort2trials_blob.m or gates it by an
%     explicit mode argument so call sites can be unambiguous.
%   - The %% IN DEVELOPMENT block inside mode 2 (event extraction from
%     video assessment) and the %% TODO FT spike trial parsing block at
%     the bottom are dead code waiting for a decision.
%
% Jesus. 29.05.2026

neurons = [];
neurons_FT = struct();

nclus = length(spike.label);

% For each cluster
for c = 1:nclus
    if iscell(trialdef)
        %% for each requested alignment
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
        % Forward per-cluster metadata onto neurons, so downstream
        % functions (calculate_fireRate_general, plotting) have access to
        % them without needing the original spike struct.
        neurons.ROI(c)          = spike.roi(c);
        if isfield(spike,'KSLabel'),     neurons.KSLabel(c)     = spike.KSLabel(c);     end
        if isfield(spike,'HumanLabel'),  neurons.HumanLabel(c)  = spike.HumanLabel(c);  end
        if isfield(spike,'bc_unitType'), neurons.bc_unitType(c) = spike.bc_unitType(c); end
        if isfield(spike,'phyLabel'),    neurons.phyLabel(c)    = spike.phyLabel(c);    end
    else
        %% Spiking indexing for Social interactions. Checks blob interaction 
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
        % % In addition, find event times and and codes at video-asessment
        % % file, for social interaction cues extracted by students at Juan's
        % % Social paradigm
        % for i=1:size(opt.alignto,1)
        %     events.(opt.alignto{i,1}) = [];
        % 
        %     for t = 1:ntrial
        %         % Grab all timestamps between time of start and time of end (inclusive)
        %         trialstamps = EventRecord.TimeSecFromMidnight(EventRecord.TimeSecFromMidnight >= trialdef{2,i}(t,1)/1000 & ...
        %                                                      EventRecord.TimeSecFromMidnight <= trialdef{2,i}(t,2)/1000);
        %         % Relativize trial timestamps to alignment offset
        %         trialstamps = trialstamps - trialdef{2,i}(t,3)/1000; 
        % 
        %         % Grab all events ocurring between time of start and time of end (inclusive)
        %         trialevents = EventRecord.EventType(EventRecord.TimeSecFromMidnight >= trialdef{2,i}(t,1)/1000 & ...
        %                                             EventRecord.TimeSecFromMidnight <= trialdef{2,i}(t,2)/1000);
        % 
        %         % Insert into the proper structure to be output.
        %         events.(opt.alignto{i,1}).code{t,1} = trialevents; 
        %         events.(opt.alignto{i,1}).time{t,1} = trialstamps; 
        %     end
        % end

    end
end
end