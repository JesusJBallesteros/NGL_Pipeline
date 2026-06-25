function sub = restrictAggregatedToSubject(aggregated, subjIdx)
% restrictAggregatedToSubject  Return a 1-row "view" of aggregated for one subject.
%
% PURPOSE:
%   loadAggregatedSpikes returns cell arrays sized (nSubj x maxSess) for
%   every field (allspike, allneurons, allcondition, ...). When NGL04
%   iterates per subject, downstream helpers (buildRequestCatSets,
%   buildFireRatePool) expect their input to be "the data they should
%   look at", not "the full study cube". This helper takes one row of
%   each cell field and returns a struct that those helpers consume
%   without any awareness of multi-subject mode.
%
% USAGE:
%   sub = restrictAggregatedToSubject(aggregated, 2);   % keep row 2
%
% INPUTS:
%   aggregated - struct from loadAggregatedSpikes.
%   subjIdx    - 1-based subject row index.
%
% OUTPUT:
%   sub        - struct with the same fieldnames as aggregated. Cell
%                fields keep only their `subjIdx`-th row (shape
%                1 x maxSess). Non-cell fields are copied verbatim.
%                If subjIdx is out of range, every cell field comes back
%                empty (cell(0,0)).
%
% Last modified 23.06.2026 (Jesus)

    assert(isstruct(aggregated), 'NGL:restrictAgg', 'aggregated must be a struct.');
    assert(isnumeric(subjIdx) && isscalar(subjIdx) && subjIdx >= 1 && subjIdx == floor(subjIdx), ...
        'NGL:restrictAgg', 'subjIdx must be a positive integer.');

    sub = struct();
    fns = fieldnames(aggregated);
    for k = 1:numel(fns)
        f = fns{k};
        v = aggregated.(f);
        if ~iscell(v)
            sub.(f) = v;
            continue
        end
        [nx, ~] = size(v);
        if subjIdx <= nx
            sub.(f) = v(subjIdx, :);
        else
            sub.(f) = cell(0, 0);
        end
    end
end
