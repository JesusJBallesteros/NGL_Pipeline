function saveFireRatePoolCache(cacheDir, key, pool)
% saveFireRatePoolCache  Persist a buildFireRatePool output to disk.
%
%   Cache file: <cacheDir>/<key>.mat (variable name: pool). Directory
%   is created on demand. Empty `cacheDir` is a no-op so the orchestrator
%   can pass '' to disable caching.
%
% Last modified 09.06.2026 (Jesus)

    if isempty(cacheDir), return; end
    if ~isfolder(cacheDir), mkdir(cacheDir); end
    fpath = fullfile(cacheDir, [key '.mat']);
    try
        save(fpath, 'pool', '-v7.3');
    catch ME
        warning('NGL:saveFireRatePoolCache:writeFail', ...
            'Could not write cache %s: %s', fpath, ME.message);
    end
end
