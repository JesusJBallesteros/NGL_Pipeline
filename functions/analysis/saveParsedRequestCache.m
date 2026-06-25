function saveParsedRequestCache(cacheDir, request, parsed)
% saveParsedRequestCache  Persist a parseFireRateRequest result.
%
%   Sibling of saveFireRatePoolCache. File:
%   <cacheDir>/<encoded-request>__parsed.mat (variable name: parsed).
%   Empty cacheDir is a no-op so callers can disable caching.
%
% Last modified 23.06.2026 (Jesus)

    if isempty(cacheDir), return; end
    if ~isfolder(cacheDir), mkdir(cacheDir); end
    fpath = fullfile(cacheDir, [encodeFireRateRequest(request) '__parsed.mat']);
    try
        save(fpath, 'parsed', '-v7.3');
    catch ME
        warning('NGL:saveParsedRequestCache:writeFail', ...
            'Could not write parsed cache %s: %s', fpath, ME.message);
    end
end
