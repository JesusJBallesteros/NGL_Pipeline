function opts = generateDefaultsFromSchema(schema)
% generateDefaultsFromSchema  Build the canonical defaults struct from
%                              an option schema (struct array of optEntry).
%
% PURPOSE:
%   Schema-driven replacement for the long hand-maintained block in
%   default_opt.m. Walks each schema entry and writes its default value
%   to the corresponding dotted path on the output struct, building
%   nested sub-structs as needed.
%
% USAGE:
%   opts = generateDefaultsFromSchema(optSchema());
%
% INPUT:
%   schema  - struct array from optSchema(), each entry has .name and
%             .default (see optEntry.m).
%
% OUTPUT:
%   opts    - struct shaped exactly like the legacy default_opt() return.
%             verifySchemaParity asserts this field-for-field.
%
% NOTES:
%   - Order of insertion matches schema order. Fieldname ordering inside
%     nested structs follows schema order too.
%   - For nested-struct defaults we always insert via dotted-name entries
%     (e.g. 'gwfparams.wfWin', 'gwfparams.nWf', ...) rather than one entry
%     with a struct default. This keeps validators per-field and lets
%     opt_help describe each leaf individually.
%
% Last modified 09.06.2026 (Jesus)

    opts = struct();
    for k = 1:numel(schema)
        opts = setNested(opts, schema(k).name, schema(k).default);
    end
end

function s = setNested(s, dottedPath, val)
    parts = strsplit(dottedPath, '.');
    if numel(parts) == 1
        s.(parts{1}) = val;
        return
    end
    head = parts{1};
    rest = strjoin(parts(2:end), '.');
    if ~isfield(s, head) || ~isstruct(s.(head))
        s.(head) = struct();
    end
    s.(head) = setNested(s.(head), rest, val);
end
