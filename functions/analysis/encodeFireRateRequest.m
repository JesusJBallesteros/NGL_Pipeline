function s = encodeFireRateRequest(request)
% encodeFireRateRequest  Filesystem-safe encoding of a 3-cell request.
%
% USAGE:
%   s = encodeFireRateRequest({'correct vs incorrect', 'good', 'stim2'});
%   -> 'correct_vs_incorrect__good__stim2'
%
% Spaces become underscores; ' vs ' becomes '_vs_'; the three entries
% are joined by '__'. Stable across calls so it is safe to use as a
% cache key for buildFireRatePool outputs.
%
% Last modified 09.06.2026 (Jesus)

    raw = cellfun(@(c) strrep(c, ' vs ', '_vs_'), request, 'uni', 0);
    s   = strrep(strjoin(raw, '__'), ' ', '_');
end
