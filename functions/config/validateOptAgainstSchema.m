function validateOptAgainstSchema(opt, schema)
% validateOptAgainstSchema  Per-entry validation of an opt struct.
%
% PURPOSE:
%   Schema-driven replacement for the long inline-assert block that used
%   to live in set_default.m Section 2. Walks every schema entry, reads
%   the corresponding leaf field of opt, and calls its validator.
%   Throws NGL:invalidOption on the first failure with a message that
%   identifies the dotted path, echoes the schema's help text, and shows
%   the offending value.
%
% USAGE:
%   validateOptAgainstSchema(opt, optSchema());
%
% INPUTS:
%   opt     - resolved opt struct (post-applyOverrides).
%   schema  - struct array from optSchema().
%
% NOTES:
%   - Cross-field rules (e.g. stepSz_ms <= binSize_ms,
%     aggregateSubjects implies aggregateSessions) belong in
%     optPostChecks.m, NOT here. Keep this function purely per-value.
%   - Identifier is always NGL:invalidOption for compatibility with
%     callers that pattern-match on identifier.
%
% Last modified 09.06.2026 (Jesus)

    for k = 1:numel(schema)
        path = schema(k).name;
        v    = localGetNested(opt, path);
        try
            ok = schema(k).validator(v);
        catch ME
            error('NGL:invalidOption', ...
                'opt.%s: validator threw an error (%s). Schema says: %s', ...
                path, ME.message, schema(k).help);
        end
        if ~ok
            error('NGL:invalidOption', ...
                'opt.%s failed validation: %s (got %s).', ...
                path, schema(k).help, localShort(v));
        end
    end
end

function v = localGetNested(s, dottedPath)
    parts = strsplit(dottedPath, '.');
    v = s;
    for k = 1:numel(parts)
        if ~isstruct(v) || ~isfield(v, parts{k})
            error('NGL:invalidOption', ...
                'opt.%s is missing (or its parent is not a struct).', dottedPath);
        end
        v = v.(parts{k});
    end
end

function s = localShort(v)
    try
        if isempty(v)
            s = '[]';
        elseif ischar(v)
            s = sprintf('''%s''', v);
        elseif iscell(v)
            s = sprintf('{1x%d cell}', numel(v));
        elseif islogical(v) || isnumeric(v)
            if isscalar(v), s = mat2str(v);
            else,           s = sprintf('%s %s', mat2str(size(v)), class(v));
            end
        elseif isstruct(v)
            s = sprintf('struct(%s)', strjoin(fieldnames(v), ', '));
        else
            s = ['<' class(v) '>'];
        end
    catch
        s = '<?>';
    end
end
