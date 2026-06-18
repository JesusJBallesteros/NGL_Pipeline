function flat = flattenAggregatedForArea(aggregated, areaName)
% flattenAggregatedForArea  Return a per-area "view" of an aggregated struct.
%
% PURPOSE:
%   Multi-area NGL01/NGL02 runs store per-(subj,sess) data as a struct
%   keyed by area, e.g.
%       aggregated.allspike{x,y}.NCL.HumanLabel{c}
%       aggregated.allneurons{x,y}.NCL.stim2{c,1}{i,1}
%   The existing helpers (buildRequestCatSets, buildFireRatePool, ...)
%   expect the flat shape used by single-area runs, i.e.
%       aggregated.allspike{x,y}.HumanLabel{c}
%       aggregated.allneurons{x,y}.stim2{c,1}{i,1}
%
%   This helper produces a view in which every cell that carries the
%   per-area nesting is drilled into for `areaName`, leaving everything
%   else untouched. Downstream code sees a normal flat aggregated and
%   needs no awareness of multi-area mode.
%
% USAGE:
%   flat = flattenAggregatedForArea(aggregated, 'NCL');
%
% INPUTS:
%   aggregated  - struct from loadAggregatedSpikes; cells may be either
%                 flat (single-area) or per-area-nested (multi-area).
%   areaName    - char, the area to drill into (one of input.Areas).
%
% OUTPUT:
%   flat        - struct of the same fieldnames as aggregated, with each
%                 cell entry replaced by entry.(areaName) when that
%                 sub-field exists. Cells that don't carry area-nesting
%                 (e.g. allcondition, alltrialdef) are passed through.
%
% NOTES:
%   - Cells that have NO matching area sub-field become empty so that
%     downstream helpers see them as "no data for this (subj, sess)".
%     That preserves the standard "skip empty cell" semantics.
%   - Non-cell fields (anything that isn't a cell array of per-(subj,sess)
%     structs) are copied verbatim.
%
% SEE ALSO:
%   detectMultiAreaFields, loadAggregatedSpikes, buildFireRatePool.
%
% Last modified 18.06.2026 (Jesus)

    assert(isstruct(aggregated), 'NGL:flattenAggregated', ...
        'aggregated must be a struct.');
    assert(ischar(areaName) && ~isempty(areaName), ...
        'NGL:flattenAggregated', 'areaName must be a non-empty char.');

    flat = struct();
    fn   = fieldnames(aggregated);
    for k = 1:numel(fn)
        f     = fn{k};
        value = aggregated.(f);
        if ~iscell(value)
            flat.(f) = value;
            continue
        end
        out = cell(size(value));
        for j = 1:numel(value)
            c = value{j};
            if isstruct(c) && isfield(c, areaName)
                out{j} = c.(areaName);
            elseif isstruct(c) && ~localLooksNested(c)
                % Cell is a flat struct (already single-area shape) —
                % pass through. Lets one aggregated mix area-nested and
                % flat fields without crashing.
                out{j} = c;
            else
                % Either empty, or a struct that's nested by area but
                % doesn't carry this particular area -> empty for this
                % (subj, sess).
                out{j} = [];
            end
        end
        flat.(f) = out;
    end
end

function tf = localLooksNested(c)
% Heuristic: a "nested" cell value is a struct whose fields are all
% themselves structs (i.e. an area->payload map). Single-area payloads
% have numeric / cell / char leaves, so they fail this test.
    fns = fieldnames(c);
    if isempty(fns), tf = false; return; end
    tf = all(cellfun(@(f) isstruct(c.(f)), fns));
end
