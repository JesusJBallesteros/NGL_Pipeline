function fireRate = calculate_fireRate_general(neurons, events, condition, opt, param)
% calculate_fireRate_general  Per-cluster firing rate for a full session,
%                             treating every trial uniformly (no blocks,
%                             no treatments).
%
% PURPOSE:
%   Computes per-trial, per-cluster firing rate (and optional baseline-
%   normalised firing rate) across an entire session, for every alignment
%   in opt.alignto. Trials are filtered by behavioural condition
%   (param.trial2plot), then handed to calcFireRate.
%
%   This is the simple, no-blocks path. If a project needs to split a
%   session into treatment blocks (e.g. basal / treatment / post), use
%   calculate_fireRate_byBlock.m instead.
%
% USAGE:
%   fireRate = calculate_fireRate_general(neurons, events, condition, opt, param)
%
% INPUTS:
%   neurons   - struct from sort2trials (standard mode). Required:
%                 .<alignment>{c,1}{i,1}  per-trial spike vectors (ms,
%                                         relative to alignment)
%                 .ROI                    1xNclust cell of area labels
%   events    - per-trial event struct from trialdefGen (not used by the
%               no-blocks path, kept in signature for parity with the
%               block-aware sibling).
%   condition - per-trial condition struct (from conditions_script). Used
%               fields:
%                 .correct    1xNtrials logical / numeric (all trials)
%                 .aborted    1xNtrials logical (aborted-trial mask)
%                 .<trial2plot>  optional, if param.trial2plot is something
%                                other than 'allInitiated' the field of
%                                that name is read as a trial mask.
%   opt       - resolved options struct (post-set_default). Required:
%                 .alignto       cell of alignment names (e.g. {'itiOn'}).
%   param     - analysis/plot tuning struct. All fields are filled with
%               safe inline defaults if absent (the function is
%               "semi-independent" by design). See PARAM DEFAULTS below.
%
% OUTPUT:
%   fireRate  - struct with cell arrays sized {Nclust, Nalign}:
%                 .sps       [Ntotal x Nbins] raw spikes/s per cell.
%                 .Norm      [Ntotal x Nbins] baseline-normalised.
%                 .meanNorm  [1 x Nbins] mean of .Norm across ALL trials.
%               *** SHAPE CONTRACT (#19 + #26, 02.06.2026) ***
%               Second cell dimension is the ALIGNMENT index, matching
%               opt.alignto. Pre-#26 code used {c,1} for every alignment
%               iteration, so only the LAST alignment's results survived.
%               Now consumers index as fireRate.sps{c, a} for cluster c,
%               alignment a in opt.alignto.
%               .sps and .Norm always have Ntotal rows, where Ntotal is
%               the per-trial dimension of `neurons.<align>`. Aborted /
%               filtered-out trials are NOT removed by this function:
%               consumers apply their own trial2plot selection (via
%               applyTrialFilter) so the row index stays aligned with
%               the per-trial condition vectors. .meanNorm is now the
%               mean across all valid (non-NaN) trials including
%               aborted ones; if you want the legacy "aborted-removed"
%               mean, compute it via applyTrialFilter + nanmean in your
%               consumer.
%
% PARAM DEFAULTS (inline; override by setting before the call):
%   .trial2plot 'allInitiated'  no longer filters rows here (#19), but
%                               downstream consumers read this field.
%   .binSize    200  ms         FR sliding-bin width
%   .stepSz      20  ms         sliding-bin step
%   .interval   [-2000 10000]   window around alignment, ms
%   .smpRate    1000  Hz        used by calcFireRate to convert counts/s
%   .baseline   -interval(1)    ms before alignment used for normalisation
%   .plot       true            also draw single- and multi-cluster plots
%
% CALLS:
%   calcFireRate
%   plot_single_fireRate (if param.plot)
%   plot_multi_fireRate (if param.plot)
%
% SEE ALSO:
%   calculate_fireRate_byBlock  iterates per pre-defined treatment/block (param.blockBounds).
%   calculate_fireRate_extintion  project-specific (Extintion paradigm).
%
% Last modified 29.05.2026 (Jesus)

%% Default options.
% this function preserves the FULL trial axis on every cluster,
% consumers (plot_fireRate_session, popDyn methods) apply their own
% selection via applyTrialFilter. We still default trial2plot here for
% downstream consumers that read the same param struct.
if ~isfield(param,'trial2plot'), param.trial2plot = 'allInitiated'; end % which trials (downstream use)
if ~isfield(param,'binSize'),    param.binSize    = 200;            end % ms
if ~isfield(param,'stepSz'),     param.stepSz     = 20;             end % ms
if ~isfield(param,'interval'),   param.interval   = [-2000 10000];  end % ms
if ~isfield(param,'smpRate'),    param.smpRate    = 1000;           end % smp/s
if ~isfield(param,'baseline'),   param.baseline   = -param.interval(1); end % ms
if ~isfield(param,'plot'),       param.plot       = true;           end % logic

%% Initialize
toalignto = opt.alignto;
param.ROI = neurons.ROI';

% Forward per-cluster curation labels to the plotting helpers as a
% subtitle (alongside ROI). Empty cells fall back to '' so the subtitle
% builder can skip missing pieces gracefully.
if isfield(neurons,'KSLabel'),     param.KSLabel     = neurons.KSLabel';     end
if isfield(neurons,'HumanLabel'),  param.HumanLabel  = neurons.HumanLabel';  end
if isfield(neurons,'bc_unitType'), param.bc_unitType = neurons.bc_unitType'; end
if isfield(neurons,'phyLabel'),    param.phyLabel    = neurons.phyLabel';    end

%% Per-alignment, per-cluster firing rate
for a = 1:length(toalignto)
    neuronSet = neurons.(toalignto{a});

    for c = 1:length(neuronSet)
        %% Build the per-trial spike set for this cluster.
        toCalculate = neuronSet{c}; % cell of vectors, one per trial

        % Trials with no spikes at all become NaN sentinels (calcFireRate
        % treats NaN as "skip this row" so they appear as gaps, not zeros).
        emptytrials = cellfun(@isempty, toCalculate);
        toCalculate(emptytrials) = {NaN};

        %% FireRate calculation. ALL trials are passed through; no
        %  row-dropping for trial2plot here. Downstream consumers apply
        %  their own selection via applyTrialFilter so the row index of
        %  fireRate.sps{c} stays aligned with the per-trial condition
        %  vectors at all times.
        param.cl = [a c 1];  % alignment, cluster, level (always 1 here)

        % Output indexed by [cluster, alignment] so multiple alignments
        % don't overwrite each other (#26). Previous code used {c,1}.
        [fireRate.sps{c,a}, fireRate.Norm{c,a}, fireRate.meanNorm{c,a}] = ...
            calcFireRate(toCalculate, param, opt);

        % Normalize cell-of-vectors -> matrix where appropriate.
        if iscell(fireRate.sps{c,a}),      fireRate.sps{c,a}      = cell2mat(fireRate.sps{c,a});      end
        if iscell(fireRate.Norm{c,a}),     fireRate.Norm{c,a}     = cell2mat(fireRate.Norm{c,a});     end
        if iscell(fireRate.meanNorm{c,a}), fireRate.meanNorm{c,a} = cell2mat(fireRate.meanNorm{c,a}); end
    end
end

end % function end
