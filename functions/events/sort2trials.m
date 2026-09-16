function [neurons, neurons_FT] = sort2trials(spike, trialdef, opt)
% sort2trials  Bin per-cluster spike timestamps into per-trial cell arrays,
%              relativised to a trial-anchor alignment event.
%
% PURPOSE:
%   Walks every (cluster, trial) pair and keeps only the spikes whose
%   timestamp falls inside the trial window defined by trialdef. Output
%   times are expressed in ms relative to the alignment event of that
%   trial.
%
%   trialdef is the (2 x nAlignments) cell produced by trialdefGen, with
%   rows {1,a}=alignment name (char), {2,a}=[nTrials x 3] start/end/t0
%   in ms. The output struct is keyed by alignment name:
%       neurons.<align>{c,1}{i,1} = vector of ms timestamps for cluster c,
%                                   trial i, relative to trialdef{2,a}(i,3).
%   Per-cluster metadata copied from spike: ROI, KSLabel, HumanLabel,
%   bc_unitType, phyLabel.
%
%   For social-arena indexing against a numeric blob-interaction matrix,
%   use the sibling sort2trials_blob.m (audit item R, split out 02.06.2026).
%
% USAGE:
%   [neurons, neurons_FT] = sort2trials(spike, trialdef, opt)
%
% INPUTS:
%   spike    - NGL spike struct from loadSpikes. Required: .label,
%              .timestamp (cell-per-cluster, SECONDS), .roi.
%   trialdef - (2 x nAlignments) cell from trialdefGen.
%   opt      - options struct. Required: opt.alignto (cell of alignment names).
%
% OUTPUTS:
%   neurons    - struct keyed by alignment name; see PURPOSE for shape.
%                Spike times in MILLISECONDS, relative to alignment.
%   neurons_FT - reserved for a FieldTrip-style spike struct. Currently
%                always returned as an empty struct(); the implementation
%                lives commented at the bottom of this file.
%
% KNOWN ISSUES:
%   - The %% TODO FT spike trial parsing block at the bottom is dead code
%     waiting for a decision.
%
% SEE ALSO:
%   sort2trials_blob (social-arena blob indexing path).
%
% Last modified 02.06.2026 (Jesus) - dropped social branch into sibling (#10 R)

% This function now handles ONLY the standard alignment-keyed case.
% Social-arena blob indexing was split out into sort2trials_blob.m as
% audit item R. A non-cell trialdef here is therefore a caller error.
assert(iscell(trialdef), 'NGL:sort2trials:badTrialdef', ...
    ['sort2trials expects trialdef to be a (2 x nAlignments) cell from ', ...
     'trialdefGen. For social-arena blob indexing (numeric trialdef matrix) ', ...
     'call sort2trials_blob.m instead.']);

neurons = [];
neurons_FT = struct();

nclus = length(spike.label);

% For each cluster
for c = 1:nclus
    %% for each requested alignment
    for a = 1:size(opt.alignto,2)
        % Match the alignment by NAME, not by position. trialdef carries its
        % own names in row 1, and a caller that asks for a subset of them (or
        % for them in another order) would otherwise get another alignment's
        % windows returned under the name it asked for - silently, and looking
        % entirely plausible. Position is kept only as the fallback for a
        % trialdef whose row 1 does not name this alignment.
        col = find(strcmp(trialdef(1,:), opt.alignto{1,a}), 1);
        if isempty(col)
            assert(a <= size(trialdef,2), 'NGL:sort2trials:noAlignment', ...
                ['alignment ''%s'' is not in trialdef, and there is no ', ...
                 'column %d to fall back to.'], opt.alignto{1,a}, a);
            col = a;
        end
        ntrial = size(trialdef{2,col},1);
        neurons.(opt.alignto{1,a}){c,1} = cell(ntrial,1);
        % for each trial
        for i=1:ntrial
            st = spike.timestamp{1,c}*1000; % convert spike times to msec

            % index for spiketimes ...
            idx = st >= trialdef{2,col}(i,1) & ... % btw trial start
                  st <  trialdef{2,col}(i,2);      % and trial end

            % relativize times to the given alignment point
            neurons.(opt.alignto{1,a}){c,1}{i,1} = st(idx) - trialdef{2,col}(i,3);
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
end
end