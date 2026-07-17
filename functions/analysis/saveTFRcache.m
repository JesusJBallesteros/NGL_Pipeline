function saveTFRcache(cacheDir, key, TFR, cfg, provenance)
% saveTFRcache  Persist a TFR + cfg + provenance under a cache key.
%
%   Sibling of saveFireRatePoolCache. File:
%     <cacheDir>/<key>.mat   with variables (TFR, cfg, provenance).
%   Empty cacheDir is a no-op so callers can disable caching.
%
% USAGE:
%   saveTFRcache(cacheDir, key, TFR, cfg, provenance);
%
% Last modified 26.06.2026 (Jesus) - new helper (LFP Pass 2).

    if isempty(cacheDir), return; end
    if ~isfolder(cacheDir), mkdir(cacheDir); end
    if nargin < 5, provenance = struct(); end
    fpath = fullfile(cacheDir, [key '.mat']);
    try
        save(fpath, 'TFR', 'cfg', 'provenance', '-v7.3');
    catch ME
        warning('NGL:saveTFRcache:writeFail', ...
            'Could not write TFR cache %s: %s', fpath, ME.message);
    end
end
