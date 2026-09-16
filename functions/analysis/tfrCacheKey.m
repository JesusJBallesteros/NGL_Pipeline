function key = tfrCacheKey(alignName, method, freqInterest, area, extra)
% tfrCacheKey  Stable, filesystem-safe key for one cached LFP TFR.
%
%   Per-TFR cache files live at <cacheDir>/<key>.mat. The key encodes
%   (alignment, TFR method, frequency signature, area, extra) so a rerun
%   with different settings never clobbers or reuses another run's data.
%   Same idiom / helpers as fireRatePoolCacheKey on the spike side.
%
% USAGE:
%   key = tfrCacheKey('stim2', 'wavelet', {[4:1:30] [30:2:150]}, 'NCL');
%   key = tfrCacheKey('stim2', 'superlet', {...}, 'main', struct('artZ', 10));
%
% INPUTS:
%   alignName    - char, alignment tag (as in opt.alignto).
%   method       - char, opt.TFRmethod ('wavelet'|'mtmconvol'|'superlet').
%   freqInterest - cell of numeric vectors (opt.freqInterest). Hashed to
%                  a compact signature: '<band>Nx<min>-<max>-<npts>'.
%   area         - char, area tag (from FT_data.chanArea; single-area
%                  runs pass 'main').
%   extra        - (optional) struct of extra key-forming fields
%                  (artZvalue, timeResol, superletOrder, ...).
%                  Its fieldnames + values are stringified and appended.
%
% Last modified 26.06.2026 (Jesus) - new helper (LFP Pass 2).

    if nargin < 5, extra = struct(); end
    freqSig = localFreqSignature(freqInterest);
    key = sprintf('TFR__%s__%s__%s__area_%s', ...
        sanitize(alignName), sanitize(method), ...
        freqSig, sanitize(area));
    if isstruct(extra) && ~isempty(fieldnames(extra))
        key = [key '__' localExtraSignature(extra)];
    end

    % Windows refuses a path over 260 characters, and MATLAB reports that
    % refusal as "the file appears to be corrupt" - so a long key does not look
    % like a naming problem at all, it looks like a broken cache. Keep the
    % readable head and fold the rest into a hash: still unique, still stable
    % across runs, and short enough to live under a deep study folder.
    maxLen = 120;
    if numel(key) > maxLen
        h = localHash(key);
        key = [key(1:maxLen - numel(h) - 2) '__' h];
    end
end

function h = localHash(s)
% FNV-1a, 32 bit: a few lines, no toolbox, and the same value on every
% machine and MATLAB version - which a cache key has to be.
    hash = uint64(2166136261);
    for c = uint64(double(s))
        hash = bitxor(hash, c);
        hash = bitand(hash * uint64(16777619), uint64(4294967295));
    end
    h = lower(dec2hex(hash, 8));
end

function sig = localFreqSignature(freqInterest)
% Compact one-liner per freq band: '<nBands>_<min>-<max>-<npts>_...'.
    if ~iscell(freqInterest) || isempty(freqInterest)
        sig = 'freqEMPTY';
        return
    end
    parts = cell(1, numel(freqInterest));
    for k = 1:numel(freqInterest)
        b = freqInterest{k};
        if isempty(b) || ~isnumeric(b)
            parts{k} = 'nan';
        else
            parts{k} = sprintf('%g-%g-%dp', min(b), max(b), numel(b));
        end
    end
    sig = ['freq' num2str(numel(parts)) 'x' strjoin(parts, '_')];
    sig = sanitize(sig);
end

function sig = localExtraSignature(extra)
    fns = fieldnames(extra);
    parts = cell(1, numel(fns));
    for k = 1:numel(fns)
        v = extra.(fns{k});
        parts{k} = [fns{k} '_' localValStr(v)];
    end
    sig = sanitize(strjoin(parts, '__'));
end

function s = localValStr(v)
    if isnumeric(v)
        s = num2str(v(:).', '%g_');
        s = regexprep(s, '_+$', '');
    elseif islogical(v)
        s = char('F' + logical(v)*('T'-'F'));
    elseif ischar(v)
        s = v;
    elseif iscell(v)
        s = sprintf('cell%d', numel(v));
    else
        s = 'x';
    end
end

function s = sanitize(s)
    s = regexprep(char(string(s)), '[^\w\-]', '_');
    if isempty(s), s = '_'; end
end
