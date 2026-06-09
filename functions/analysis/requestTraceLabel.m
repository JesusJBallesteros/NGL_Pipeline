function s = requestTraceLabel(varying, parsed, cIdx, lIdx)
% requestTraceLabel  Legend label for one overlaid trace.
%
%   Returns the parts of (condition, label) that VARY within a subplot.
%   Alignment is never part of a trace label since alignment varies
%   between subplots, not within them. If neither condition nor label
%   varies, falls back to whichever single value is present so the
%   legend still shows context.
%
% Last modified 09.06.2026 (Jesus)

    parts = {};
    if any(strcmp(varying, 'condition'))
        parts{end+1} = parsed.condition{cIdx};
    end
    if any(strcmp(varying, 'label'))
        parts{end+1} = parsed.label{lIdx};
    end
    if isempty(parts)
        if numel(parsed.condition) == 1, parts{end+1} = parsed.condition{1}; end
        if numel(parsed.label)     == 1, parts{end+1} = parsed.label{1};     end
    end
    s = strjoin(parts, ' | ');
end
