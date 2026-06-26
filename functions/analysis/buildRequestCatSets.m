function catSets = buildRequestCatSets(aggregated, opt)
% buildRequestCatSets  Build the category sets used to categorise tokens
%                      in a request cell (NGL04_fireRate / NGL04_PCA).
%
% PURPOSE:
%   Walks the aggregated data once to enumerate every distinct
%   condition fieldname and every distinct cluster-label value seen
%   across all (subj, sess) cells, then combines them with the
%   user-configured alignment list and label-priority order. Returned
%   struct is consumed by parseFireRateRequest.
%
% USAGE:
%   catSets = buildRequestCatSets(aggregated, opt);
%
% CONTRACT:
%   `aggregated` MUST be the FLAT single-area shape produced by
%   loadAggregatedSpikes(input, area). Each cell of aggregated.allspike
%   is a spike struct with HumanLabel / KSLabel / bc_unitType / phyLabel
%   directly as fields (NOT nested by area). This was always the assumed
%   shape; up to 26.06.2026 multi-area runs accidentally passed area-
%   nested cells in here and labelPool came back empty. The NGL03 26.06
%   refactor writes one aggregated_<area>.mat per area so callers always
%   see this flat shape now and there is no multi-area branch in here.
%
% OUTPUT FIELDS:
%   .alignments     opt.alignto verbatim
%   .conditions     cell of unique fieldnames from aggregated.allcondition
%   .labelPool      struct with .HumanLabel/.KSLabel/.bc_unitType/.phyLabel,
%                   each a cell of unique string values seen across
%                   aggregated.allspike
%   .labelPriority  resolution order for label tokens (defaults to
%                   {'HumanLabel','KSLabel','bc_unitType'})
%
% Last modified 26.06.2026 (Jesus) - docstring nails down flat-shape
%                                     contract; per-area aggregated
%                                     layout from NGL03 means there is
%                                     never anything else to feed in.

    catSets = struct();
    catSets.alignments    = opt.alignto;
    catSets.conditions    = localCollectConditionFields(aggregated);
    % 'allInitiated' is a magic condition token recognised by
    % applyTrialFilter / buildFireRatePool (= ~aborted). Inject it so
    % parseFireRateRequest accepts it as a valid condition entry.
    if ~ismember('allInitiated', catSets.conditions)
        catSets.conditions = [{'allInitiated'}; catSets.conditions(:)];
    end
    catSets.labelPool     = localCollectClusterLabels(aggregated);
    if isfield(opt,'fireRatePlot') && isfield(opt.fireRatePlot,'labelPriority')
        catSets.labelPriority = opt.fireRatePlot.labelPriority;
    else
        catSets.labelPriority = {'HumanLabel','KSLabel','bc_unitType'};
    end
end

function fields = localCollectConditionFields(aggregated)
% Union of fieldnames seen across aggregated.allcondition cells.
    fields = {};
    if ~isfield(aggregated, 'allcondition'), return; end
    cells = aggregated.allcondition;
    for k = 1:numel(cells)
        c = cells{k};
        if isempty(c) || ~isstruct(c), continue, end
        fields = union(fields, fieldnames(c));
    end
    fields = fields(:);
end

function labelPool = localCollectClusterLabels(aggregated)
% Build a struct keyed by label-field name, each value the cell of
% unique non-empty string values observed across aggregated.allspike.
    labelPool = struct( ...
        'HumanLabel',  {{}}, ...
        'KSLabel',     {{}}, ...
        'bc_unitType', {{}}, ...
        'phyLabel',    {{}});
    if ~isfield(aggregated, 'allspike'), return; end
    cells = aggregated.allspike;
    poolFields = fieldnames(labelPool)';
    for k = 1:numel(cells)
        spk = cells{k};
        if isempty(spk) || ~isstruct(spk), continue, end
        for f = poolFields
            field = f{1};
            if ~isfield(spk, field) || ~iscell(spk.(field)), continue, end
            vals = spk.(field);
            vals = vals(~cellfun(@isempty, vals));
            valsChar = cellfun(@(v) char(string(v)), vals, 'uni', false);
            labelPool.(field) = union(labelPool.(field), valsChar);
        end
    end
end
