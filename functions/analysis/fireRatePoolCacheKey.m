function key = fireRatePoolCacheKey(alignName, condField, labelValue, labelField)
% fireRatePoolCacheKey  Stable, filesystem-safe key for one cached pool.
%
%   Per-pool cache files live at <cacheDir>/<key>.mat. The key encodes
%   the (alignment, condition, label_value, label_field) quadruple that
%   buildFireRatePool uses as its filter. The same pool is reused by
%   NGL04_fireRate and NGL04_PCA.
%
% USAGE:
%   key = fireRatePoolCacheKey('stim2', 'correct', 'good', 'HumanLabel');
%
% Last modified 09.06.2026 (Jesus)

    key = sprintf('%s__%s__%s_%s', ...
        sanitize(alignName), sanitize(condField), ...
        sanitize(labelValue), sanitize(labelField));
end

function s = sanitize(s)
    s = regexprep(char(string(s)), '[^\w\-]', '_');
    if isempty(s), s = '_'; end
end
