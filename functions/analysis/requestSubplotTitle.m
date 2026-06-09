function s = requestSubplotTitle(parsed, aIdx)
% requestSubplotTitle  Subplot title for an aligned-data figure.
%
%   Shows the current alignment plus any non-alignment factor that is
%   NOT varying (i.e. fixed across the whole figure). Varying non-
%   alignment factors (condition, label) belong in the legend, not the
%   title.
%
% Last modified 09.06.2026 (Jesus)

    parts = {};
    parts{end+1} = parsed.alignment{aIdx};
    if numel(parsed.condition) == 1, parts{end+1} = parsed.condition{1}; end
    if numel(parsed.label)     == 1, parts{end+1} = parsed.label{1};     end
    s = strjoin(parts, ' | ');
end
