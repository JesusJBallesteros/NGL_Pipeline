function FT_data = ensureChanArea(FT_data, areaMap)
% ensureChanArea  Guarantee FT_data.chanArea{i} is populated per channel.
%
% PURPOSE:
%   Multi-area LFP analyses want to restrict processing to a subset of
%   channels by area label without needing a separate FT file per area
%   (Q3 answer, 26.06.2026 LFP audit: option (b) - one FT file, per-
%   channel area tags). This helper is the single source of truth for
%   attaching the tags:
%       FT_data.chanArea{i} = <area label>   for i = 1:numel(FT_data.label)
%
%   Idempotent: if the field already exists with the right length, the
%   function is a no-op. Legacy FT files created before 2026-06 have no
%   chanArea; this helper backfills at load time (called from NGL02_LFP
%   / NGL07_LFPanalysis).
%
% USAGE:
%   FT_data = ensureChanArea(FT_data, areaMap);    % explicit map (or [])
%
% INPUTS:
%   FT_data - FieldTrip data struct with .label {nChan x 1 cell}.
%   areaMap - either:
%               * struct from buildAreaMap (fields .uniqueAreas cell,
%                 .chanIdx_per_area cell of numeric vectors)
%               * empty ([] or struct()) -> every channel gets 'main'.
%
% OUTPUT:
%   FT_data - same struct, guaranteed to carry .chanArea {nChan x 1 cell}
%             of area-label chars.
%
% NOTES:
%   * If a channel index appears in multiple areas' chanIdx_per_area
%     (should not happen, but be defensive), the LAST area wins and a
%     warning is emitted.
%   * If a channel index falls outside every area's coverage, its label
%     becomes 'unassigned'. Downstream code that filters by area should
%     drop 'unassigned' channels or handle them explicitly.
%
% SEE ALSO:
%   buildAreaMap (source of areaMap struct), MAT2FieldTrip (bakes chanArea
%   into new FT files at NGL01 save time), NGL02_LFP (backfills at load).
%
% Last modified 26.06.2026 (Jesus) - new helper (LFP Pass 2).

    assert(isstruct(FT_data) && isfield(FT_data, 'label'), ...
        'NGL:ensureChanArea:badFT', 'FT_data must be a FieldTrip struct with .label.');
    nCh = numel(FT_data.label);

    % Idempotency short-circuit.
    if isfield(FT_data, 'chanArea') && iscell(FT_data.chanArea) ...
            && numel(FT_data.chanArea) == nCh ...
            && all(cellfun(@(c) ischar(c) && ~isempty(c), FT_data.chanArea))
        return
    end

    % Default fill = 'main'.
    chanArea = repmat({'main'}, nCh, 1);

    % Multi-area case.
    hasMap = ~isempty(areaMap) && isstruct(areaMap) ...
             && isfield(areaMap, 'uniqueAreas') && ~isempty(areaMap.uniqueAreas) ...
             && isfield(areaMap, 'chanIdx_per_area');
    if hasMap
        % Track which channels remain unassigned so we can warn.
        assigned = false(nCh, 1);
        for a = 1:numel(areaMap.uniqueAreas)
            areaName = areaMap.uniqueAreas{a};
            idx      = areaMap.chanIdx_per_area{a};
            idx      = idx(idx >= 1 & idx <= nCh);
            if any(assigned(idx))
                warning('NGL:ensureChanArea:overlap', ...
                    ['Channel(s) [%s] appear in more than one area''s chanIdx_per_area. ', ...
                     'Last-writer-wins: area ''%s'' is being applied.'], ...
                    num2str(idx(assigned(idx))'), areaName);
            end
            chanArea(idx) = {areaName};
            assigned(idx) = true;
        end
        if any(~assigned)
            chanArea(~assigned) = {'unassigned'};
        end
    end

    FT_data.chanArea = chanArea;
end
