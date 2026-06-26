function saveParsedRequestCache(cacheDir, request, area, parsed)
% saveParsedRequestCache  Persist a parseFireRateRequest result.
%
%   Sibling of saveFireRatePoolCache. File:
%   <cacheDir>/<encoded-request>[__area_<NAME>]__parsed.mat
%   (variable name: parsed). Empty cacheDir is a no-op so callers can
%   disable caching.
%
% USAGE:
%   saveParsedRequestCache(cacheDir, request, area, parsed)
%
% INPUTS:
%   cacheDir - char, per-subject cache folder.
%   request  - 1x3 cell, same one passed to parseFireRateRequest.
%   area     - char, area tag (or '' for single-area). Included in the
%              filename so per-area parsed results don't collide.
%   parsed   - struct from parseFireRateRequest.
%
% Last modified 26.06.2026 (Jesus) - per-area cache key (added `area` arg).

    if isempty(cacheDir), return; end
    if ~isfolder(cacheDir), mkdir(cacheDir); end
    fpath = fullfile(cacheDir, [localParsedCacheName(request, area) '.mat']);
    try
        save(fpath, 'parsed', '-v7.3');
    catch ME
        warning('NGL:saveParsedRequestCache:writeFail', ...
            'Could not write parsed cache %s: %s', fpath, ME.message);
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
