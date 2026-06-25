function [parsed, hit] = loadParsedRequestCache(cacheDir, request, sourceFile)
% loadParsedRequestCache  Read a cached `parsed` struct keyed by request.
%
% PURPOSE:
%   parseFireRateRequest needs catSets, which is built from aggregated
%   data. Caching the parsed result per (subject, request) lets NGL04
%   skip the expensive aggregated.mat load entirely when every pool the
%   request needs is also already cached. Same staleness contract as
%   loadFireRatePoolCache: if sourceFile is provided and newer than the
%   cache, treat as miss.
%
% USAGE:
%   [parsed, hit] = loadParsedRequestCache(cacheDir, request);
%   [parsed, hit] = loadParsedRequestCache(cacheDir, request, sourceFile);
%
% INPUTS:
%   cacheDir    - char, the per-subject cache folder.
%   request     - 1x3 cell (same one passed to parseFireRateRequest).
%   sourceFile  - (optional) absolute path to aggregated.mat for
%                 mtime staleness check.
%
% OUTPUTS:
%   parsed - the parsed struct on hit, [] on miss.
%   hit    - logical.
%
% Last modified 23.06.2026 (Jesus)

    parsed = [];
    hit    = false;
    if isempty(cacheDir) || ~isfolder(cacheDir), return; end
    fpath = fullfile(cacheDir, [encodeFireRateRequest(request) '__parsed.mat']);
    if ~isfile(fpath), return; end

    if nargin >= 3 && ~isempty(sourceFile) && isfile(sourceFile)
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
