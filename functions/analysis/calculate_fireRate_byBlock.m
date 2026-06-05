function fireRate = calculate_fireRate_byBlock(neurons, events, condition, opt, param)
% calculate_fireRate_byBlock  Per-cluster firing rate split across one
%                             treatment/block partition of the session.
%
% PURPOSE:
%   Same firing-rate machinery as calculate_fireRate_general, but the
%   session is divided into N blocks (typically corresponding to a
%   pharmacological / behavioural treatment schedule: basal / treatment /
%   post, or similar) and firing rate is computed independently within
%   each block. Use this when you need to compare FR across phases of a
%   session.
%
%   Block boundaries are passed in via param.blockBounds and are the
%   caller's responsibility. The function does not try to infer blocks
%   from a project-specific condition field, on purpose: that detection
%   logic differs per paradigm and belongs in the user's NGL_SetAndRunMe
%   (or in a project-specific helper they call before this function).
%
% USAGE:
%   % 1) Decide block partition (project-specific). Example using a
%   %    'treatment' code in condition that bumps from 0 -> 1 -> 2 ...
%   bumps = find(diff(condition.v) > 0) + 1;
%   nTrials = numel(condition.correct);
%   param.blockBounds = [1, bumps, nTrials + 1];   % N+1 indices, N blocks
%
%   % 2) Call.
%   fireRate = calculate_fireRate_byBlock(neurons, events, condition, opt, param);
%
% INPUTS:
%   neurons   - struct from sort2trials (standard mode); see
%               calculate_fireRate_general header.
%   events    - per-trial event struct from trialdefGen (unused by this
%               function, kept for signature parity).
%   condition - per-trial condition struct. Used fields:
%                 .correct       1xNtrials  (sets the total trial count)
%                 .aborted       1xNtrials logical (aborted-trial mask)
%                 .<trial2plot>  optional, see param.trial2plot below.
%   opt       - resolved options struct. Required: opt.alignto.
%   param     - analysis/plot tuning struct. REQUIRED field:
%                 .blockBounds  1xN+1 numeric, strictly increasing trial
%                               indices delimiting N blocks. Block b runs
%                               over trials param.blockBounds(b) ..
%                               param.blockBounds(b+1)-1 (inclusive of
%                               the first index, exclusive of the second).
%               Other fields default exactly as in
%               calculate_fireRate_general (see PARAM DEFAULTS).
%
% OUTPUT:
%   fireRate  - struct with cell arrays sized {Nclust, Nblocks}:
%                 .sps{c,b}      raw spikes/s for cluster c, block b
%                 .Norm{c,b}     baseline-normalised
%                 .meanNorm{c,b} mean of .Norm across the block's trials
%
% PARAM DEFAULTS (same as calculate_fireRate_general):
%   .trial2plot 'allInitiated'
%   .binSize    200  ms
%   .stepSz      20  ms
%   .interval   [-2000 10000]
%   .smpRate    1000 Hz
%   .baseline   -interval(1)  ms
%   .plot       true   (plotting is per-cluster across blocks, then
%                       session-level multi-cluster like _general)
%
% NOTES:
%   - The previous in-line "treatment / level" logic that lived commented
%     in calculate_fireRate_general has been folded into this function.
%   - Block-internal trial selection uses the same param.trial2plot
%     filter as _general; the block partition is applied on top.
%   - Hierarchical block-on-block analyses (e.g. treatment x sub-block)
%     are not supported in a single call. For now, call this function
%     twice with different param.blockBounds or write a paradigm-specific
%     wrapper.
%
% CALLS:
%   calcFireRate, plot_single_fireRate (if param.plot),
%   plot_multi_fireRate (if param.plot).
%
% SEE ALSO:
%   calculate_fireRate_general    no-blocks version.
%   calculate_fireRate_extintion  legacy project-specific implementation;
%                                 a candidate for being rewritten as a
%                                 thin wrapper around this function plus
%                                 Extintion-paradigm block detection.
%
% Last modified 29.05.2026 (Jesus) - extracted from _general (#13)

%% Validate the required field.
assert(isfield(param,'blockBounds') && isnumeric(param.blockBounds) && ...
       isvector(param.blockBounds) && numel(param.blockBounds) >= 2 && ...
       all(diff(param.blockBounds) > 0), ...
    'NGL:calculate_fireRate_byBlock:badBlockBounds', ...
    ['param.blockBounds must be a strictly increasing vector of trial ', ...
     'indices of length nBlocks+1.']);

nBlocks = numel(param.blockBounds) - 1;

%% Default options (same as _general).
if ~isfield(param,'trial2plot'), param.trial2plot = 'allInitiated'; end
if ~isfield(param,'binSize'),    param.binSize    = 200;            end
if ~isfield(param,'stepSz'),     param.stepSz     = 20;             end
if ~isfield(param,'interval'),   param.interval   = [-2000 10000];  end
if ~isfield(param,'smpRate'),    param.smpRate    = 1000;           end
if ~isfield(param,'baseline'),   param.baseline   = -param.interval(1); end
if ~isfield(param,'plot'),       param.plot       = true;           end

%% Initialize
toalignto = opt.alignto;
param.ROI = neurons.ROI';

% Forward per-cluster curation labels (same as calculate_fireRate_general).
if isfield(neurons,'KSLabel'),     param.KSLabel     = neurons.KSLabel';     end
if isfield(neurons,'HumanLabel'),  param.HumanLabel  = neurons.HumanLabel';  end
if isfield(neurons,'bc_unitType'), param.bc_unitType = neurons.bc_unitType'; end
if isfield(neurons,'phyLabel'),    param.phyLabel    = neurons.phyLabel';    end

% Expose blockBounds to calcFireRate's normalisation logic, which expects
% it under the name 'blockchange' (legacy). We pass interior boundaries
% only (drop the first/last sentinels), matching the historical convention.
param.blockchange = param.blockBounds(2:end-1);

%% Per-alignment, per-cluster, per-block firing rate
for a = 1:length(toalignto)
    neuronSet = neurons.(toalignto{a});

    for c = 1:length(neuronSet)
        %% Build the per-trial spike set for this cluster (full session).
        % We make one full-session selector then mask down per block.
        clusterTrials = neuronSet{c};

        % Trials with no spikes at all -> NaN sentinels.
        emptytrials = cellfun(@isempty, clusterTrials);
        clusterTrials(emptytrials) = {NaN};

        % Behavioural filter (applied to all blocks identically) via the
        % shared applyTrialFilter helper (#19). NOTE: unlike
        % calculate_fireRate_general, _byBlock's output IS inherently
        % per-block-shaped, so we keep the per-block row reduction here
        % via the block-membership mask + the filter mask.
        validMask = applyTrialFilter(condition, param.trial2plot);
        clusterTrials(~validMask) = {[]};

        %% Per-block FR
        for b = 1:nBlocks
            % Trials belonging to this block.
            trialIdx = false(size(clusterTrials));
            blockTrials = param.blockBounds(b):(param.blockBounds(b+1) - 1);
            blockTrials = blockTrials(blockTrials <= numel(clusterTrials));
            trialIdx(blockTrials) = true;

            % Combine with non-empty mask (skip filtered-out trials).
            spikes2use = trialIdx(:)' & ~cellfun(@isempty, clusterTrials);

            param.cl = [a c b];  % alignment, cluster, block

            % Output indexed by [cluster, alignment, block] so multiple
            % alignments don't overwrite each other (#26). Pre-#26 the
            % alignment loop was ignored in indexing.
            [fireRate.sps{c,a,b}, fireRate.Norm{c,a,b}, fireRate.meanNorm{c,a,b}] = ...
                calcFireRate(clusterTrials(spikes2use), param, opt);

            if iscell(fireRate.sps{c,a,b}),      fireRate.sps{c,a,b}      = cell2mat(fireRate.sps{c,a,b});      end
            if iscell(fireRate.Norm{c,a,b}),     fireRate.Norm{c,a,b}     = cell2mat(fireRate.Norm{c,a,b});     end
            if iscell(fireRate.meanNorm{c,a,b}), fireRate.meanNorm{c,a,b} = cell2mat(fireRate.meanNorm{c,a,b}); end

            if param.plot
                plot_single_fireRate(fireRate.sps{c,a,b}, fireRate.Norm{c,a,b}, param, opt)
            end
        end
    end
end

%% Session-level multi-cluster plot (across all blocks).
if param.plot
    plot_multi_fireRate(fireRate.meanNorm, param, opt)
end

end % function end
