function [ok, report] = verifySchemaParity(legacyDefaultsFile)
% verifySchemaParity  Confirm that the schema-driven defaults match the
%                     legacy hand-maintained default_opt return.
%
% PURPOSE:
%   Run this once before swapping default_opt.m to the schema wrapper.
%   Walks every leaf field of the legacy struct and the schema-built
%   struct in parallel; reports any divergence. Returns true if and only
%   if the two structs are field-for-field identical.
%
% USAGE:
%   ok = verifySchemaParity()                  % uses default_opt() as legacy
%   ok = verifySchemaParity('default_opt_old.m') % compare to a frozen copy
%   [ok, report] = verifySchemaParity()        % capture textual report
%
% INPUT:
%   legacyDefaultsFile  - (optional) path to a legacy default_opt-style
%                         function file. Defaults to 'default_opt'. Useful
%                         when comparing against a frozen pre-swap copy.
%
% OUTPUT:
%   ok       - logical: true iff fully equal.
%   report   - cellstr of issue messages (empty when ok=true).
%
% NOTES:
%   - This is a one-off verification helper, not part of the runtime
%     pipeline. Safe to delete once the swap is in.
%
% Last modified 09.06.2026 (Jesus)

    if nargin < 1 || isempty(legacyDefaultsFile)
        legacyFn = @default_opt;
    else
        [p, n, ~] = fileparts(legacyDefaultsFile);
        if ~isempty(p), addpath(p); end
        legacyFn = str2func(n);
    end

    legacy = legacyFn();
    schema = optSchema();
    fresh  = generateDefaultsFromSchema(schema);

    report = {};
    report = localCompare(legacy, fresh, '', report);
    ok = isempty(report);

    if ok
        fprintf('verifySchemaParity: PASS - schema defaults match legacy default_opt.\n');
    else
        fprintf('verifySchemaParity: %d divergence(s):\n', numel(report));
        for k = 1:numel(report)
            fprintf('  %s\n', report{k});
        end
    end
end

function r = localCompare(a, b, prefix, r)
    % Recurse into structs; for leaves use isequaln.
    if isstruct(a) && isstruct(b)
        af = fieldnames(a);
        bf = fieldnames(b);
        onlyA = setdiff(af, bf);
        onlyB = setdiff(bf, af);
        for k = 1:numel(onlyA)
            r{end+1} = sprintf('missing in schema: %s%s', prefix, onlyA{k}); %#ok<AGROW>
        end
        for k = 1:numel(onlyB)
            r{end+1} = sprintf('extra in schema:   %s%s', prefix, onlyB{k}); %#ok<AGROW>
        end
        common = intersect(af, bf);
        for k = 1:numel(common)
            r = localCompare(a.(common{k}), b.(common{k}), ...
                             [prefix common{k} '.'], r);
        end
        return
    end
    if isstruct(a) || isstruct(b)
        r{end+1} = sprintf('type mismatch at %s: %s vs %s', ...
                           strip_trailing_dot(prefix), class(a), class(b));
        return
    end
    if ~isequaln(a, b)
        r{end+1} = sprintf('value differs at %s: legacy=%s, schema=%s', ...
                           strip_trailing_dot(prefix), ...
                           localShort(a), localShort(b));
    end
end

function s = strip_trailing_dot(s)
    if ~isempty(s) && s(end) == '.', s = s(1:end-1); end
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
        else
            s = ['<' class(v) '>'];
        end
    catch
        s = '<?>';
    end
end
