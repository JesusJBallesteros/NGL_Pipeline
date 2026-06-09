function result = calculate_pca_from_pool(pools, condLabels, params)
% calculate_pca_from_pool  Population PCA from buildFireRatePool output.
%
% PURPOSE:
%   Cross-subject / cross-session version of calculate_neural_pca for
%   NGL04_PCA. The single-session function consumed a fireRate struct
%   from one session; here we consume a CELL of pool structs (one per
%   condition level) that already pool clusters across (subject, session).
%   We bin each cluster's trial spike vectors into firing rates, build
%   the [Nclust x Nbins x Ncond] condition-mean tensor, fit PCA, and
%   return:
%     - the condition-mean trajectories in PC space (the solid lines)
%     - per-session mean trajectories projected into the same PC space
%       (the grey overlay - "how does the average generalize?")
%     - a trial-bootstrap 95% CI tube around each condition mean (the
%       alternative overlay).
%
% USAGE:
%   result = calculate_pca_from_pool(pools, condLabels, params);
%
% INPUTS:
%   pools       {1 x Ncond} cell of pool structs from buildFireRatePool.
%               All pools must share the same alignment, labelField and
%               labelValue; only the condition filter differs. Each pool
%               brings its .byCluster array; clusters are matched across
%               pools by (subjIdx, sessIdx, clusterIdx).
%   condLabels  {1 x Ncond} cell of condition names (char). Used as
%               legend labels; their order matches `pools`.
%   params      struct:
%                 .intervalMs       [pre post] window, ms relative to alignment
%                 .binSize_ms       bin width, ms (e.g. 200)
%                 .stepSz_ms        bin step, ms (e.g. 20)
%                 .smoothSigma_s    Gaussian sigma in SECONDS (e.g. 0.050)
%                 .nComponents      number of PCs to keep (>=2)
%                 .nBootstrap       bootstrap reps for the CI tube
%                                   (0 disables CI). Default 100.
%                 .smpRate          spike-time sample rate (1000 for ms)
%                 .rngSeed          (optional) integer for reproducibility
%
% OUTPUT:
%   result      struct:
%                 .method        'PCA-pool'
%                 .coeff         [Nclust x K] loadings
%                 .mu            [1 x Nclust] per-neuron mean removed before PCA
%                 .explained     [K x 1] percent variance
%                 .timeAxis      [1 x Nbins] seconds, t=0 at alignment
%                 .conditions    {Ncond x 1} cell of condLabels
%                 .clusters      Nclust struct array of cluster identity
%                                (subjIdx, sessIdx, clusterIdx, sessionKey, area)
%                 .meanFR        [Nclust x Nbins x Ncond] condition-mean rates
%                                (smoothed; raw inputs to PCA after reshape)
%                 .traj_mean     [Nbins x K x Ncond] condition-mean
%                                trajectories in PC space (solid lines)
%                 .traj_session  {Ncond x 1} cell, each [Nbins x K x Nsess]
%                                per-session contribution (faint grey overlay).
%                                Per-session traces sum to the condition mean.
%                 .sessionKeys   {Nsess x 1} cell of session-key strings,
%                                indexed identically across all .traj_session
%                                slices.
%                 .nClustPerSess [Nsess x 1] number of clusters in each session
%                 .ciLo          [Nbins x K x Ncond] 5th percentile of bootstrap
%                                projected condition mean (NaN if nBootstrap==0)
%                 .ciHi          [Nbins x K x Ncond] 95th percentile (NaN if 0)
%                 .params        echoed params
%
% NOTES ON CROSS-SESSION HANDLING:
%   The "per-session" trajectory is the MARGINAL contribution of that
%   session's clusters to the condition-mean projection. Because the PCA
%   was fit on the full pool, each session's projection (computed using
%   only its rows of coeff against the centered rates) sums across
%   sessions to the condition-mean trajectory. Session traces will have
%   amplitudes proportional to how many clusters that session contributed
%   - large sessions push the mean more than small ones. This is the
%   most rigorous "session generalization" view; if amplitude-matched
%   comparison is preferred later, the orchestrator can rescale.
%
%   The bootstrap CI tube is a TRIAL bootstrap (resample trials per
%   cluster, recompute condition mean, project through the original
%   coeff). This captures within-cluster trial-level uncertainty in the
%   mean trajectory. Session-level variability is shown by the per-
%   session traces instead.
%
% SEE ALSO:
%   buildFireRatePool, calculate_neural_pca, smooth_spikes,
%   calcFireRate (toolboxes/BDPAT_NGL/), plot_pca_state_space.
%
% Last modified 09.06.2026 (Jesus)

    %% Defaults
    if ~isfield(params,'nBootstrap'),    params.nBootstrap    = 100;   end
    if ~isfield(params,'smpRate'),       params.smpRate       = 1000;  end
    if ~isfield(params,'smoothSigma_s'), params.smoothSigma_s = 0.050; end
    if ~isfield(params,'rngSeed'),       params.rngSeed       = [];    end

    intervalMs    = params.intervalMs;
    binSize_ms    = params.binSize_ms;
    stepSz_ms     = params.stepSz_ms;
    smoothSigma_s = params.smoothSigma_s;
    K             = params.nComponents;
    nBoot         = params.nBootstrap;

    Ncond = numel(pools);
    assert(Ncond >= 1, 'NGL:calculate_pca_from_pool:noPools', ...
        'pools must be a non-empty cell of pool structs.');
    if numel(condLabels) ~= Ncond
        error('NGL:calculate_pca_from_pool:condLabelMismatch', ...
            'condLabels (n=%d) must match pools (n=%d).', numel(condLabels), Ncond);
    end

    %% Build the union cluster set across pools.
    % Each cluster is keyed by (subjIdx, sessIdx, clusterIdx). Clusters
    % missing in a given pool contribute zero rate to that condition
    % (handled by leaving meanFR(k,:,c) at the NaN sentinel and replacing
    % NaN with 0 at the smoothing stage).
    [clusters, perPoolIdx] = localUnionClusters(pools);
    Nclust = numel(clusters);
    if Nclust < 2
        error('NGL:calculate_pca_from_pool:tooFewClusters', ...
            'Need at least 2 clusters in the union pool; got %d.', Nclust);
    end

    %% Bin parameters.
    param = struct( ...
        'stepSz',   stepSz_ms,    ...
        'binSize',  binSize_ms,   ...
        'interval', intervalMs,   ...
        'smpRate',  params.smpRate);
    % Nbins is the number of sliding windows that fit in the interval.
    % calcFireRate's interior: windowBorder = interval(1):stepSz:interval(2)
    % -> Nbins = length(windowBorder) - 1.
    edges  = intervalMs(1):stepSz_ms:intervalMs(2);
    Nbins  = numel(edges) - 1;
    timeAxis_s = (edges(1:end-1) + binSize_ms/2) / 1000;   % bin centers in seconds

    %% Build [Nclust x Nbins x Ncond] condition-mean rates.
    % We also stash per-cluster, per-condition trial-by-bin matrices for
    % the bootstrap stage (kept as a cell to handle ragged Ntrials).
    trialMats = cell(Nclust, Ncond);     % {k,c} = [Ntrials_kc x Nbins] or []
    meanFR    = nan(Nclust, Nbins, Ncond);
    for c = 1:Ncond
        pool = pools{c};
        for k = 1:Nclust
            pIdx = perPoolIdx(k, c);
            if pIdx == 0, continue; end
            cl = pool.byCluster(pIdx);
            if isempty(cl.trials), continue; end
            spikes = cl.trials(:);
            empt   = cellfun(@isempty, spikes);
            spikes(empt) = {NaN};
            fr  = calcFireRate(spikes, param, []);   % no baseline -> single output only
            mat = fr{1};                              % [Ntrials x Nbins]
            trialMats{k, c} = mat;
            meanFR(k, :, c) = mean(mat, 1, 'omitnan');
        end
    end

    %% Replace remaining NaN rates with zero (missing clusters in some
    %  conditions, or all-NaN trial sets) so PCA has a finite matrix.
    meanFR(isnan(meanFR)) = 0;

    %% Gaussian smoothing along time.
    binSize_s = stepSz_ms / 1000;        % bin spacing for kernel sizing
    meanFR_sm = smooth_spikes(meanFR, smoothSigma_s, binSize_s);

    %% PCA fit on [Nbins*Ncond x Nclust].
    X = reshape(permute(meanFR_sm, [2 3 1]), Nbins * Ncond, Nclust);
    Kavail = min(K, min(size(X)));
    if Kavail < K
        warning('NGL:calculate_pca_from_pool:fewerComponents', ...
            'Requested %d components but only %d available; using %d.', K, Kavail, Kavail);
    end
    [coeff, scoreMat, ~, ~, explained, mu] = pca(X, 'NumComponents', Kavail);

    traj_mean = reshape(scoreMat, Nbins, Ncond, Kavail);
    traj_mean = permute(traj_mean, [1 3 2]);     % [Nbins x K x Ncond]

    %% Per-session marginal projections.
    sessionKeys   = unique({clusters.sessionKey}, 'stable')';
    Nsess         = numel(sessionKeys);
    sessClustMask = false(Nclust, Nsess);
    for s = 1:Nsess
        sessClustMask(:, s) = strcmp({clusters.sessionKey}, sessionKeys{s});
    end
    nClustPerSess = sum(sessClustMask, 1)';

    traj_session = cell(Ncond, 1);
    for c = 1:Ncond
        per_sess = nan(Nbins, Kavail, Nsess);
        Yc = squeeze(meanFR_sm(:, :, c))';       % [Nbins x Nclust]
        Yc_centered = Yc - mu;                   % broadcast over rows
        for s = 1:Nsess
            sel = sessClustMask(:, s);
            if ~any(sel), continue; end
            per_sess(:, :, s) = Yc_centered(:, sel) * coeff(sel, :);
        end
        traj_session{c} = per_sess;
    end

    %% Bootstrap CI tube (trial bootstrap).
    ciLo = nan(Nbins, Kavail, Ncond);
    ciHi = nan(Nbins, Kavail, Ncond);
    if nBoot > 0
        if ~isempty(params.rngSeed), rng(params.rngSeed); end
        bootScores = nan(nBoot, Nbins, Kavail, Ncond);
        for b = 1:nBoot
            meanFR_b = zeros(Nclust, Nbins, Ncond);
            for c = 1:Ncond
                for k = 1:Nclust
                    mat = trialMats{k, c};
                    if isempty(mat), continue; end
                    Ntr = size(mat, 1);
                    idx = randi(Ntr, Ntr, 1);
                    meanFR_b(k, :, c) = mean(mat(idx, :), 1, 'omitnan');
                end
            end
            meanFR_b(isnan(meanFR_b)) = 0;
            meanFR_b = smooth_spikes(meanFR_b, smoothSigma_s, binSize_s);
            Xb = reshape(permute(meanFR_b, [2 3 1]), Nbins * Ncond, Nclust);
            sB = (Xb - mu) * coeff;              % [Nbins*Ncond x Kavail]
            sB = reshape(sB, Nbins, Ncond, Kavail);
            sB = permute(sB, [1 3 2]);           % [Nbins x Kavail x Ncond]
            bootScores(b, :, :, :) = sB;
        end
        ciLo = squeeze(prctile(bootScores,  5, 1));
        ciHi = squeeze(prctile(bootScores, 95, 1));
        if Ncond == 1
            ciLo = reshape(ciLo, Nbins, Kavail, 1);
            ciHi = reshape(ciHi, Nbins, Kavail, 1);
        end
    end

    %% Pack result.
    result.method        = 'PCA-pool';
    result.coeff         = coeff;
    result.mu            = mu;
    result.explained     = explained(1:Kavail);
    result.timeAxis      = timeAxis_s;
    result.conditions    = condLabels(:);
    result.clusters      = clusters;
    result.meanFR        = meanFR_sm;
    result.traj_mean     = traj_mean;
    result.traj_session  = traj_session;
    result.sessionKeys   = sessionKeys;
    result.nClustPerSess = nClustPerSess;
    result.ciLo          = ciLo;
    result.ciHi          = ciHi;
    result.params        = params;
end

% ------------------------------------------------------------------------
function [clusters, perPoolIdx] = localUnionClusters(pools)
% Build a struct array of unique (subjIdx, sessIdx, clusterIdx) triples
% across all pools. perPoolIdx(k, c) is the row of pools{c}.byCluster
% corresponding to global cluster k, or 0 if that pool has no entry.
    Ncond  = numel(pools);
    keys   = cell(0, 1);
    src    = struct('subjIdx', {}, 'sessIdx', {}, 'sessionKey', {}, ...
                    'clusterIdx', {}, 'area', {});
    perPoolIdxRows = zeros(0, Ncond);

    for c = 1:Ncond
        bc = pools{c}.byCluster;
        for j = 1:numel(bc)
            tag = sprintf('%d_%d_%d', bc(j).subjIdx, bc(j).sessIdx, bc(j).clusterIdx);
            idx = find(strcmp(keys, tag), 1);
            if isempty(idx)
                idx = numel(keys) + 1;
                keys{idx, 1}              = tag;
                src(idx).subjIdx          = bc(j).subjIdx;
                src(idx).sessIdx          = bc(j).sessIdx;
                src(idx).sessionKey       = bc(j).sessionKey;
                src(idx).clusterIdx       = bc(j).clusterIdx;
                src(idx).area             = bc(j).area;
                perPoolIdxRows(idx, :)    = 0;
            end
            perPoolIdxRows(idx, c) = j;
        end
    end

    clusters   = src(:);
    perPoolIdx = perPoolIdxRows;
end
