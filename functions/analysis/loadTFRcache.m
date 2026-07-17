function [TFR, cfg, hit] = loadTFRcache(cacheDir, key, sourceFile)
% loadTFRcache  Read a cached TFR struct keyed by `key`.
%
% PURPOSE:
%   Sibling of loadFireRatePoolCache on the spike side. Returns the
%   TFR + cfg persisted by saveTFRcache, if the entry exists and is not
%   stale relative to sourceFile.
%
% USAGE:
%   [TFR, cfg, hit] = loadTFRcache(cacheDir, key);
%   [TFR, cfg, hit] = loadTFRcache(cacheDir, key, sourceFile);
%
% INPUTS:
%   cacheDir   - directory holding <key>.mat files (empty / missing dir
%                yields hit = false).
%   key        - from tfrCacheKey(...).
%   sourceFile - (optional) absolute path to the FT file this TFR was
%                derived from. If supplied and newer than the cache,
%                the cache is treated as stale (miss).
%
% OUTPUTS:
%   TFR - struct on hit, [] on miss.
%   cfg - ft_freqanalysis cfg used, or [].
%   hit - logical.
%
% Last modified 26.06.2026 (Jesus) - new helper (LFP Pass 2).

    TFR = [];  cfg = [];  hit = false;
    if isempty(cacheDir) || ~isfolder(cacheDir), return; end
    fpath = fullfile(cacheDir, [key '.mat']);
    if ~isfile(fpath), return; end

    if nargin >= 3 && ~isempty(sourceFile) && isfile(sourceFile)
        sInfo = dir(sourceFile);
        cInfo = dir(fpath);
        if ~isempty(sInfo) && ~isempty(cInfo) && cInfo.datenum < sInfo.datenum
            return  % stale
        end
    end

    try
        S = load(fpath, 'TFR', 'cfg');
        if isfield(S, 'TFR')
            TFR = S.TFR;
            if isfield(S, 'cfg'), cfg = S.cfg; end
            hit = true;
        end
    catch ME
        warning('NGL:loadTFRcache:loadFail', ...
            'Could not read TFR cache %s: %s', fpath, ME.message);
    end
end
