function key = fireRatePoolCacheKey(alignName, condField, labelValue, labelField, area)
% fireRatePoolCacheKey  Stable, filesystem-safe key for one cached pool.
%
%   Per-pool cache files live at <cacheDir>/<key>.mat. The key encodes
%   the (alignment, condition, label_value, label_field) quadruple that
%   buildFireRatePool uses as its filter, plus an optional area tag for
%   multi-area runs. The same pool is reused by NGL04_fireRate and
%   NGL04_PCA so the cache key must be consistent across both scripts.
%
% USAGE:
%   key = fireRatePoolCacheKey('stim2', 'correct', 'good', 'HumanLabel');
%   key = fireRatePoolCacheKey('stim2', 'correct', 'good', 'HumanLabel', 'NCL');
%
% INPUTS:
%   alignName, condField, labelValue, labelField - the standard quadruple.
%   area  - (optional) area tag (one of input.Areas). '' or omitted ->
%           single-area key, matches the legacy filename so old caches
%           keep being reused. A non-empty area gets suffixed as
%           '__area_<NAME>' so per-area pools don't collide.
%
% Last modified 18.06.2026 (Jesus) - added optional `area` for multi-area.

    if nargin < 5, area = ''; end

    key = sprintf('%s__%s__%s_%s', ...
        sanitize(alignName), sanitize(condField), ...
        sanitize(labelValue), sanitize(labelField));

    if ~isempty(area)
        key = sprintf('%s__area_%s', key, sanitize(area));
    end
end

function s = sanitize(s)
    s = regexprep(char(string(s)), '[^\w\-]', '_');
    if isempty(s), s = '_'; end
end
