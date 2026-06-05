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
%   fireRate  - struct with cell arrays sized {Nclust, 1}:
%                 .sps       raw spikes/s per trial-bin (matrix per cell)
%                 .Norm      baseline-normalised, per trial-bin
%                 .meanNorm  mean of .Norm across trials, per cell
%               Indexing is kept at {c,1} (column singleton) so downstream
%               code that iterates `for p = 1:size(fireRate.sps,2)`
%               continues to work even though there is only one level.
%
% PARAM DEFAULTS (inline; override by setting before the call):
%   .trial2plot 'allInitiated'  which trials enter the FR calculation
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
if ~isfield(param,'trial2plot'), param.trial2plot = 'allInitiated'; end % which trials
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

        % Behavioural filter. All trials are valid by default; trial2plot
        % can narrow to a named condition field, except 'allInitiated'
        % which means "drop aborted trials only".
        if strcmp(param.trial2plot, 'allInitiated')
            toCalculate(logical(condition.aborted)) = {[]};
        else
            assert(isfield(condition, param.trial2plot), ...
                'NGL:calculate_fireRate_general:unknownTrialFilter', ...
                'param.trial2plot = ''%s'' but condition has no such field.', ...
                param.trial2plot);
            toCalculate(logical(~condition.(param.trial2plot))) = {[]};
        end

        %% FireRate calculation. Empty cells are skipped.
        param.cl  = [a c 1];  % alignment, cluster, level (always 1 here)
        spikes2use = ~cellfun(@isempty, toCalculate);

        [fireRate.sps{c,1}, fireRate.Norm{c,1}, fireRate.meanNorm{c,1}] = ...
            calcFireRate(toCalculate(spikes2use), param, opt);

        % Normalize cell-of-vectors -> matrix where appropriate.
        if iscell(fireRate.sps{c,1}),      fireRate.sps{c,1}      = cell2mat(fireRate.sps{c,1});      end
        if iscell(fireRate.Norm{c,1}),     fireRate.Norm{c,1}     = cell2mat(fireRate.Norm{c,1});     end
        if iscell(fireRate.meanNorm{c,1}), fireRate.meanNorm{c,1} = cell2mat(fireRate.meanNorm{c,1}); end
    end
end
% Plotting decoupled: the per-cluster and session-level plots that used
% to be drawn inside this function are now in plot_fireRate_session.m
% (audit item S). NGL02_postPhy calls that helper after this function.

end % function end
