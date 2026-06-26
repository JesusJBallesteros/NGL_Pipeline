function [parsed, hit] = loadParsedRequestCache(cacheDir, request, area, sourceFile)
% loadParsedRequestCache  Read a cached `parsed` struct keyed by request and area.
%
% PURPOSE:
%   parseFireRateRequest needs catSets, which is built from aggregated
%   data. Caching the parsed result per (subject, request, area) lets
%   NGL04 skip the per-area aggregated load entirely when every pool
%   the request needs is also already cached. Same staleness contract
%   as loadFireRatePoolCache: if sourceFile is provided and newer than
%   the cache, treat as miss.
%
% USAGE:
%   [parsed, hit] = loadParsedRequestCache(cacheDir, request, area);
%   [parsed, hit] = loadParsedRequestCache(cacheDir, request, area, sourceFile);
%
% INPUTS:
%   cacheDir    - char, the per-subject cache folder.
%   request     - 1x3 cell (same one passed to parseFireRateRequest).
%   area        - char, area tag (one of input.Areas, or '' / 'main' for
%                 single-area). Different areas may produce different
%                 parsed results (each area has its own label vocabulary),
%                 so the cache key includes it.
%   sourceFile  - (optional) absolute path to the per-area aggregated
%                 file for mtime staleness check.
%
% OUTPUTS:
%   parsed - the parsed struct on hit, [] on miss.
%   hit    - logical.
%
% Last modified 26.06.2026 (Jesus) - per-area cache key (added `area` arg).

    parsed = [];
    hit    = false;
    if nargin < 3, area = ''; end
    if isempty(cacheDir) || ~isfolder(cacheDir), return; end
    fpath = fullfile(cacheDir, [localParsedCacheName(request, area) '.mat']);
    if ~isfile(fpath), return; end

    if nargin >= 4 && ~isempty(sourceFile) && isfile(sourceFile)
        sInfo = dir(sourceFile);
        cInfo = dir(fpath);
        if ~isempty(sInfo) && ~isempty(cInfo) && cInfo.datenum < sInfo.datenum
            return  % stale
        end
    end

    try
        S = load(fpath, 'parsed');
        if isfield(S, 'parsed')
            parsed = S.parsed;
            hit    = true;
        end
    catch ME
        warning('NGL:loadParsedRequestCache:loadFail', ...
            'Could not read parsed cache %s: %s', fpath, ME.message);
    end
end

function name = localParsedCacheName(request, area)
% Filename body: <encoded-request>[__area_<NAME>]__parsed
    base = encodeFireRateRequest(request);
    if isempty(area)
        name = [base '__parsed'];
    else
        ak   = char(area);
        ak   = regexprep(ak, '[^\w\-]', '_');
        name = [base '__area_' ak '__parsed'];
    end
end
