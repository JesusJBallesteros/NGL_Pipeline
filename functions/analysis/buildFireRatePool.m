function pool = buildFireRatePool(aggregated, alignName, condField, labelValue, labelField)
% buildFireRatePool  Pool per-trial spike vectors across (subj, sess, clust)
%                    that match an (alignment, condition, label) triple.
%
% PURPOSE:
%   Walks aggregated.allspike{x,y} and accumulates every cluster whose
%   spike.(labelField){c} matches `labelValue`. For each matched cluster,
%   pulls the per-trial spike vectors from allneurons{x,y}.(alignName){c},
%   filters by allcondition{x,y}.(condField), and appends them to a
%   flat list (used by PSTH plotters) AND to a per-cluster list (used
%   by PCA / population-dynamics consumers that need per-neuron identity).
%
% USAGE:
%   pool = buildFireRatePool(aggregated, alignName, condField, labelValue, labelField);
%
% INPUTS:
%   aggregated  - struct from loadAggregatedSpikes
%   alignName   - char, alignment fieldname (e.g. 'stim2')
%   condField   - char, condition fieldname (e.g. 'correct')
%   labelValue  - char, cluster-label value to match (e.g. 'good')
%   labelField  - char, which label field carries the value
%                 ('HumanLabel'|'KSLabel'|'bc_unitType'|'phyLabel')
%
% OUTPUT (struct):
%   .cells       {Ntrials x 1} flat list of spike-time vectors (ms).
%                Trials from different clusters are concatenated;
%                identity is lost. PSTH plotters consume this directly.
%   .nSubj       number of distinct subjects contributing >=1 cluster
%   .nSess       number of distinct (subj, sess) pairs contributing
%   .nClust      number of clusters contributing
%   .nTrials     numel(cells)
%   .waveforms   {nClust x 1} cell of mean waveforms (where available)
%   .byCluster   {nClust x 1} struct array, ONE ENTRY PER CLUSTER:
%                  .subjIdx     row index into aggregated
%                  .sessIdx     column index into aggregated
%                  .sessionKey  sprintf('%d_%d', subjIdx, sessIdx)
%                  .clusterIdx  cluster index within its session
%                  .labelField  echoed input
%                  .labelValue  echoed input
%                  .area        spike.ROI{c} if available, else ''
%                  .trials      {Ntrials_c x 1} cell of spike vectors
%                               (condition-filtered, same content as
%                               the slice contributed to .cells)
%                  .nTrials     numel(trials)
%                  .waveform    mean waveform vector, or []
%
% CONTRACT:
%   The flat .cells field is provided strictly for backward
%   compatibility with PSTH callers. PCA callers should use .byCluster.
%
% Last modified 09.06.2026 (Jesus)

    pool = struct( ...
        'cells',     {{}},  ...
        'nSubj',     0,     ...
        'nSess',     0,     ...
        'nClust',    0,     ...
        'nTrials',   0,     ...
        'waveforms', {{}},  ...
        'byCluster', struct([]));

    if ~isfield(aggregated,'allspike') || ~isfield(aggregated,'allneurons') ...
            || ~isfield(aggregated,'allcondition')
        warning('NGL:buildFireRatePool:missingAggField', ...
            'Aggregated file lacks one of allspike/allneurons/allcondition; pool empty.');
        return
    end

    [nSubj, nSess] = size(aggregated.allspike);
    warnedAlign    = false(nSubj, nSess);
    warnedCond     = false(nSubj, nSess);

    cells       = {};
    waveforms   = {};
    byCluster   = struct( ...
        'subjIdx',    {}, ...
        'sessIdx',    {}, ...
        'sessionKey', {}, ...
        'clusterIdx', {}, ...
        'labelField', {}, ...
        'labelValue', {}, ...
        'area',       {}, ...
        'trials',     {}, ...
        'nTrials',    {}, ...
        'waveform',   {});
    subjFlag    = false(nSubj, 1);
    sessFlag    = false(nSubj, nSess);

    for x = 1:nSubj
        for y = 1:nSess
            spk = localSafeIdx(aggregated.allspike,     x, y);
            neu = localSafeIdx(aggregated.allneurons,   x, y);
            cnd = localSafeIdx(aggregated.allcondition, x, y);
            if isempty(spk) || isempty(neu) || isempty(cnd), continue, end
            if ~isstruct(spk) || ~isfield(spk, labelField), continue, end
            if ~isstruct(neu) || ~isfield(neu, alignName)
                if ~warnedAlign(x, y)
                    warning('NGL:buildFireRatePool:missingAlign', ...
                        'allneurons{%d,%d} missing alignment ''%s''; skipping.', ...
                        x, y, alignName);
                    warnedAlign(x, y) = true;
                end
                continue
            end
            % condField can be a real condition field OR the special
            % token 'allInitiated' (= keep non-aborted trials). For the
            % token we require cnd.aborted instead of cnd.allInitiated.
            if strcmpi(condField, 'allInitiated')
                requiredCond = 'aborted';
            else
                requiredCond = condField;
            end
            if ~isstruct(cnd) || ~isfield(cnd, requiredCond)
                if ~warnedCond(x, y)
                    warning('NGL:buildFireRatePool:missingCond', ...
                        'allcondition{%d,%d} missing field ''%s''; skipping.', ...
                        x, y, requiredCond);
                    warnedCond(x, y) = true;
                end
                continue
            end

            mask = applyTrialFilter(cnd, condField);
            if ~isfield(spk,'label') || isempty(spk.label), continue, end
            Nclust = numel(spk.label);

            for c = 1:Nclust
                if c > numel(spk.(labelField)), continue, end
                lab = spk.(labelField){c};
                if isempty(lab) || ~strcmp(char(string(lab)), labelValue), continue, end
                if c > numel(neu.(alignName)), continue, end
                trialCells = neu.(alignName){c};
                if isempty(trialCells), continue, end
                if numel(trialCells) ~= numel(mask)
                    warning('NGL:buildFireRatePool:shapeMismatch', ...
                        'Trial count mismatch at (%d,%d) cluster %d (neurons %d vs condition %d); skipping.', ...
                        x, y, c, numel(trialCells), numel(mask));
                    continue
                end
                filtered = trialCells(mask);
                filtered = filtered(:);
                cells    = [cells; filtered];

                subjFlag(x)    = true;
                sessFlag(x, y) = true;

                % Per-cluster waveform (mean across spikes if a matrix).
                wf = [];
                if isfield(spk,'waveform') && c <= numel(spk.waveform) && ~isempty(spk.waveform{c})
                    raw = spk.waveform{c};
                    if iscell(raw), raw = raw{1}; end
                    if isnumeric(raw) && ~isempty(raw)
                        if isvector(raw)
                            wf = raw(:);
                        else
                            wf = mean(raw, 2, 'omitnan');
                        end
                    end
                end
                if ~isempty(wf), waveforms{end+1, 1} = wf; end %#ok<AGROW>

                % Per-cluster area tag (set by loadSpikes via opt.area).
                area = '';
                if isfield(spk,'ROI') && iscell(spk.ROI) && c <= numel(spk.ROI)
                    raw = spk.ROI{c};
                    if ~isempty(raw), area = char(string(raw)); end
                end

                k = numel(byCluster) + 1;
                byCluster(k).subjIdx    = x;
                byCluster(k).sessIdx    = y;
                byCluster(k).sessionKey = sprintf('%d_%d', x, y);
                byCluster(k).clusterIdx = c;
                byCluster(k).labelField = labelField;
                byCluster(k).labelValue = labelValue;
                byCluster(k).area       = area;
                byCluster(k).trials     = filtered;
                byCluster(k).nTrials    = numel(filtered);
                byCluster(k).waveform   = wf;
            end
        end
    end

    pool.cells     = cells;
    pool.nSubj     = sum(subjFlag);
    pool.nSess     = sum(sessFlag(:));
    pool.nClust    = numel(byCluster);
    pool.nTrials   = numel(cells);
    pool.waveforms = waveforms;
    pool.byCluster = byCluster;
end

function v = localSafeIdx(arr, x, y)
    [nx, ny] = size(arr);
    if x > nx || y > ny, v = []; else, v = arr{x, y}; end
end
