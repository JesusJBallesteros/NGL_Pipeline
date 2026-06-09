function [pool, hit] = loadFireRatePoolCache(cacheDir, key, sourceFile)
% loadFireRatePoolCache  Read a cached pool struct keyed by `key`.
%
% PURPOSE:
%   Returns the pool produced by a previous buildFireRatePool call, if
%   available and not stale. Stale means the source aggregated file is
%   NEWER than the cache entry (NGL03 must have been re-run since).
%
% USAGE:
%   [pool, hit] = loadFireRatePoolCache(cacheDir, key);
%   [pool, hit] = loadFireRatePoolCache(cacheDir, key, sourceFile);
%
% INPUTS:
%   cacheDir   - char, directory holding <key>.mat files. Empty / missing
%                directory yields hit = false.
%   key        - char, from fireRatePoolCacheKey(...).
%   sourceFile - (optional) absolute path to aggregated.mat. If supplied
%                and exists, its modification time is compared with the
%                cache entry's; older cache => treated as miss.
%
% OUTPUTS:
%   pool       - the pool struct on hit, [] on miss.
%   hit        - logical.
%
% Last modified 09.06.2026 (Jesus)

    pool = [];
    hit  = false;
    if isempty(cacheDir) || ~isfolder(cacheDir), return; end
    fpath = fullfile(cacheDir, [key '.mat']);
    if ~isfile(fpath), return; end

    if nargin >= 3 && ~isempty(sourceFile) && isfile(sourceFile)
        sInfo = dir(sourceFile);
        cInfo = dir(fpath);
        if ~isempty(sInfo) && ~isempty(cInfo) && cInfo.datenum < sInfo.datenum
            return  % stale; treat as miss
        end
    end

    try
        S = load(fpath, 'pool');
        if isfield(S, 'pool')
            pool = S.pool;
            hit  = true;
        end
    catch ME
        warning('NGL:loadFireRatePoolCache:loadFail', ...
            'Could not read cache %s: %s', fpath, ME.message);
    end
end
