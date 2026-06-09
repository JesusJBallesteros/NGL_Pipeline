function pool = buildSingleSessionPool(neurons, condition, alignName, condField, labelValue, labelField)
% buildSingleSessionPool  Single-session sibling of buildFireRatePool.
%
% PURPOSE:
%   Produce a pool struct in the SAME schema as buildFireRatePool out of
%   ONE session's (neurons, condition) pair. Lets calculate_pca_from_pool
%   and plot_pca_state_space process per-session data through the exact
%   same path used by the cross-subject NGL04_PCA, including the
%   single-trial overlay (the per-trial Churchland/Mante view).
%
% USAGE:
%   pool = buildSingleSessionPool(neurons, condition, alignName, ...
%                                 condField, labelValue, labelField);
%
% INPUTS:
%   neurons    - struct from sort2trials. Required:
%                  .<alignName>{c,1}   {Ntrials x 1} cell of spike vectors (ms)
%                  .<labelField>{c}    per-cluster label value (char)
%                  .ROI{c}             (optional) area tag
%   condition  - per-trial condition struct. Must contain .aborted for the
%                'allInitiated' token; else the named field.
%   alignName  - char, alignment fieldname on `neurons`.
%   condField  - char, condition fieldname or 'allInitiated'.
%   labelValue - char, cluster-label value to match.
%   labelField - char, which label field carries the value
%                ('HumanLabel'|'KSLabel'|'bc_unitType'|'phyLabel').
%
% OUTPUT (struct, matches buildFireRatePool's schema):
%   .cells       {Ntrials x 1} flat list of spike vectors
%   .nSubj       always 1
%   .nSess       always 1
%   .nClust
%   .nTrials
%   .waveforms   {nClust x 1}, empty if neurons has no waveform field
%   .byCluster   struct array, one entry per matched cluster:
%                  .subjIdx = 1
%                  .sessIdx = 1
%                  .sessionKey = '1_1'
%                  .clusterIdx = c (index into neurons.<labelField>)
%                  .labelField, .labelValue
%                  .area
%                  .trials, .nTrials
%                  .waveform
%
% NOTES:
%   - Returns an empty-but-valid pool (nClust = 0) when no cluster in this
%     session matches labelValue. calculate_pca_from_pool detects the
%     empty case upstream and errors with a clear message.
%   - Per-cluster waveforms are not on `neurons`; they live on the spike
%     struct loaded by loadSpikes. Pass them in via an extension if
%     needed; for now the waveform field stays empty here.
%
% SEE ALSO:
%   buildFireRatePool, calculate_pca_from_pool, plot_pca_state_space,
%   applyTrialFilter.
%
% Last modified 09.06.2026 (Jesus)

    pool = struct( ...
        'cells',     {{}},  ...
        'nSubj',     1,     ...
        'nSess',     1,     ...
        'nClust',    0,     ...
        'nTrials',   0,     ...
        'waveforms', {{}},  ...
        'byCluster', struct([]));

    if ~isstruct(neurons)
        warning('NGL:buildSingleSessionPool:badNeurons', ...
            'neurons must be a struct; returning empty pool.');
        return
    end
    if ~isfield(neurons, alignName)
        warning('NGL:buildSingleSessionPool:missingAlign', ...
            'neurons has no alignment ''%s''; returning empty pool.', alignName);
        return
    end
    if ~isfield(neurons, labelField)
        warning('NGL:buildSingleSessionPool:missingLabel', ...
            'neurons has no label field ''%s''; returning empty pool.', labelField);
        return
    end

    % Condition mask via shared helper (handles 'allInitiated').
    if strcmpi(condField, 'allInitiated')
        requiredCond = 'aborted';
    else
        requiredCond = condField;
    end
    if ~isstruct(condition) || ~isfield(condition, requiredCond)
        warning('NGL:buildSingleSessionPool:missingCond', ...
            'condition has no field ''%s''; returning empty pool.', requiredCond);
        return
    end
    mask = applyTrialFilter(condition, condField);

    labels   = neurons.(labelField);
    alignArr = neurons.(alignName);
    Nclust   = numel(labels);

    cells     = {};
    waveforms = {};
    byCluster = struct( ...
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

    for c = 1:Nclust
        if c > numel(labels), continue, end
        lab = labels{c};
        if isempty(lab) || ~strcmp(char(string(lab)), labelValue), continue, end
        if c > numel(alignArr), continue, end
        trialCells = alignArr{c};
        if isempty(trialCells), continue, end
        if numel(trialCells) ~= numel(mask)
            warning('NGL:buildSingleSessionPool:shapeMismatch', ...
                'Trial count mismatch at cluster %d (neurons %d vs condition %d); skipping.', ...
                c, numel(trialCells), numel(mask));
            continue
        end

        filtered = trialCells(mask);
        filtered = filtered(:);
        cells    = [cells; filtered]; %#ok<AGROW>

        area = '';
        if isfield(neurons,'ROI') && iscell(neurons.ROI) && c <= numel(neurons.ROI)
            raw = neurons.ROI{c};
            if ~isempty(raw), area = char(string(raw)); end
        end

        k = numel(byCluster) + 1;
        byCluster(k).subjIdx    = 1;
        byCluster(k).sessIdx    = 1;
        byCluster(k).sessionKey = '1_1';
        byCluster(k).clusterIdx = c;
        byCluster(k).labelField = labelField;
        byCluster(k).labelValue = labelValue;
        byCluster(k).area       = area;
        byCluster(k).trials     = filtered;
        byCluster(k).nTrials    = numel(filtered);
        byCluster(k).waveform   = [];
    end

    pool.cells     = cells;
    pool.nClust    = numel(byCluster);
    pool.nTrials   = numel(cells);
    pool.waveforms = waveforms;
    pool.byCluster = byCluster;
end
