function stem = lfpResultName(kind, opt, varargin)
% lfpResultName  The file-name rule for every LFP output, in one place.
%
% PURPOSE:
%   Names that sort sensibly and can be parsed back: the session first, then
%   the analysis, then what distinguishes this instance of it. Used by both
%   saveLFPresult and saveLFPfigure so the .mat and the .png of one analysis
%   differ only in extension.
%
% USAGE:
%   stem = lfpResultName('TFR', opt, 'area', 'NCL', 'align', 'stim2')
%     -> 'ABC_20260101_LFP_TFR_NCL_stim2'
%
% INPUTS:
%   kind - char, the analysis name.
%   opt  - options struct; uses opt.SavFileName as the session stem.
%   Name/value pairs: 'area', 'align', 'tags' (cellstr), 'session' (override).
%
% OUTPUT:
%   stem - char, no extension. Parts are joined with '_'; anything unsafe for
%          a file name is replaced with '-' so a condition label like
%          'correct vs error' cannot break the path.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 0).

    p = inputParser;
    p.addParameter('area',    '');
    p.addParameter('align',   '');
    p.addParameter('tags',    {});
    p.addParameter('session', '');
    p.parse(varargin{:});
    a = p.Results;

    session = a.session;
    if isempty(session)
        session = localOpt(opt, {'SavFileName'}, 'session');
    end
    tags = a.tags;
    if ischar(tags) || isstring(tags), tags = {char(tags)}; end

    parts = [{session, 'LFP', char(kind), char(a.area), char(a.align)}, ...
             cellfun(@char, tags(:)', 'UniformOutput', false)];
    parts = parts(~cellfun(@isempty, parts));
    stem = strjoin(cellfun(@localSafe, parts, 'UniformOutput', false), '_');
end

function s = localSafe(s)
    s = regexprep(char(s), '[^A-Za-z0-9\-\.]+', '-');
    s = regexprep(s, '-+', '-');
    s = regexprep(s, '^-|-$', '');
end

function v = localOpt(opt, path, default)
    v = default;
    s = opt;
    for k = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{k}), return; end
        s = s.(path{k});
    end
    if ~isempty(s) || ischar(s), v = s; end
end
