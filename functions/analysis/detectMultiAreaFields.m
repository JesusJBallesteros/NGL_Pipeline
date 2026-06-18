function areas = detectMultiAreaFields(aggregated, candidateAreas)
% detectMultiAreaFields  Return the list of areas present in aggregated.allspike.
%
% PURPOSE:
%   Used by NGL04_fireRate / NGL04_PCA to decide whether to iterate over
%   areas. If aggregated.allspike's cells nest per-area sub-structs
%   matching any entry in candidateAreas (typically input.Areas), this
%   returns the unique area names present. Otherwise returns {} (the
%   data is flat / single-area and the script processes it in one pass).
%
% USAGE:
%   areas = detectMultiAreaFields(aggregated, input.Areas);
%
% INPUTS:
%   aggregated      - struct from loadAggregatedSpikes.
%   candidateAreas  - cell of area names to check for (e.g. input.Areas
%                     such as {'NCL','NCL','STR'}).
%
% OUTPUT:
%   areas - cell of unique area names actually present in aggregated.
%           Empty if no nesting is detected.
%
% Last modified 18.06.2026 (Jesus)

    areas = {};
    if ~isfield(aggregated, 'allspike'),       return; end
    if nargin < 2 || isempty(candidateAreas),  return; end
    if ~iscell(candidateAreas),                return; end

    uniqueCands = unique(candidateAreas, 'stable');

    for k = 1:numel(aggregated.allspike)
        spk = aggregated.allspike{k};
        if isempty(spk) || ~isstruct(spk), continue, end
        fns     = fieldnames(spk);
        present = intersect(fns, uniqueCands, 'stable');
        if ~isempty(present)
            % Preserve canonical order from candidateAreas.
            areas = uniqueCands(ismember(uniqueCands, present));
            return
        end
        % First non-empty allspike cell doesn't carry any candidate area
        % as a fieldname -> data is flat. Stop scanning.
        return
    end
end
