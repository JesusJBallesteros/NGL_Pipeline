function merged = applyOverrides(base, overrides)
% applyOverrides  Deep-merge `overrides` onto `base`, returning a new struct.
%
% PURPOSE:
%   The current set_default Section 1 does a FLAT merge that
%   wholesale-replaces nested struct fields. As a result, setting
%   `opt.popDyn.do = true` from NGL_SetAndRunMe replaces the entire
%   default popDyn struct with a 1-field {do:true}, and later
%   validators crash because opt.popDyn.pca etc. are gone.
%
%   applyOverrides recurses into nested structs so user-provided fields
%   override defaults at any depth, while unprovided fields keep their
%   defaults. Strictly more permissive than the flat merge; existing
%   SetAndRunMe files that set full structs continue to work.
%
% USAGE:
%   opt = applyOverrides(default_opt(), userOpt);
%
% INPUTS:
%   base       - struct (typically defaults).
%   overrides  - struct (typically user opts). Empty struct is fine.
%
% OUTPUT:
%   merged     - struct with .field == overrides.field where set, else
%                base.field. Nested structs are recursed into.
%
% MERGE RULES:
%   - Numeric / logical / char / cell values: overrides wins.
%   - Struct values: recursive merge.
%   - Field present in base but absent in overrides: base wins.
%   - Field present in overrides but absent in base: warned and dropped
%     (matches the "unknown option" behaviour of legacy set_default).
%
% Last modified 09.06.2026 (Jesus)

    if isempty(overrides) || ~isstruct(overrides)
        merged = base;
        return
    end
    merged = base;
    fns = fieldnames(overrides);
    for k = 1:numel(fns)
        f = fns{k};
        if ~isfield(base, f)
            warning('NGL:unknownOption', ...
                'Option ''%s'' is not recognised and will be ignored.', f);
            continue
        end
        if isstruct(base.(f)) && isstruct(overrides.(f))
            merged.(f) = applyOverrides(base.(f), overrides.(f));
        else
            merged.(f) = overrides.(f);
        end
    end
end
