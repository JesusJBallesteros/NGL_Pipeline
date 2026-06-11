function e = optEntry(name, default, validator, group, help)
% optEntry  Constructor for one entry in the NGL option schema.
%
% PURPOSE:
%   The schema returned by optSchema() is a struct array; each entry
%   bundles everything we need to know about ONE option (canonical name,
%   default, validator, group, help text). Schema consumers
%   (generateDefaultsFromSchema, validateOptAgainstSchema, opt_help)
%   walk the array.
%
% USAGE:
%   e = optEntry('popDyn.pcaConditions', {'allInitiated'}, ...
%                @(v) iscell(v) && ~isempty(v) && all(cellfun(@ischar, v)), ...
%                'PopDyn', ...
%                'cell of condition tokens for per-session PCA iteration');
%
% FIELDS:
%   .name       char, dotted path (e.g. 'popDyn.pcaConditions').
%   .default    any value (cell, struct, numeric, logical, char). Used by
%               generateDefaultsFromSchema to populate opts.
%   .validator  @(v) -> logical. Returns true iff v is acceptable for this
%               option. Cross-field rules belong in optPostChecks, NOT here.
%   .group      char, used by opt_help for grouping. Free-form.
%   .help       char, one-line description shown by opt_help and in any
%               generated SetAndRunMe template.
%
% NOTES:
%   - Field assignment is done via direct dot-syntax (not struct(...))
%     so cell / struct defaults stay verbatim without {{...}} gymnastics.
%
% Last modified 09.06.2026 (Jesus)

    e = struct();
    e.name      = name;
    e.default   = default;
    e.validator = validator;
    e.group     = group;
    e.help      = help;
end
