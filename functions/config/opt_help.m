function opt_help(name)
% opt_help  Print help text for one option (or list all by group).
%
% USAGE:
%   opt_help                        % print every option grouped by .group
%   opt_help('popDyn.pcaConditions')% print one option's help + default
%   opt_help('PopDyn')              % print every entry in a group
%
% Last modified 09.06.2026 (Jesus)

    schema = optSchema();
    if nargin < 1 || isempty(name)
        localListAll(schema);
        return
    end

    % Exact dotted-name match first.
    hit = find(strcmp({schema.name}, name), 1);
    if ~isempty(hit)
        localPrintEntry(schema(hit));
        return
    end

    % Group match.
    grp = find(strcmpi({schema.group}, name));
    if ~isempty(grp)
        fprintf('Group: %s  (%d options)\n', schema(grp(1)).group, numel(grp));
        for k = grp
            localPrintEntry(schema(k));
        end
        return
    end

    % Substring match as a fallback (helps with typos / partial names).
    sub = find(contains({schema.name}, name, 'IgnoreCase', true));
    if ~isempty(sub)
        fprintf('No exact match for ''%s''. %d substring match(es):\n', name, numel(sub));
        for k = sub
            localPrintEntry(schema(k));
        end
        return
    end

    fprintf('opt_help: ''%s'' not found in schema.\n', name);
end

function localListAll(schema)
    [groups, ~, gidx] = unique({schema.group}, 'stable');
    for g = 1:numel(groups)
        idx = find(gidx == g);
        fprintf('\n[%s]\n', groups{g});
        for k = idx(:)'
            fprintf('  opt.%-32s  %s\n', schema(k).name, schema(k).help);
        end
    end
    fprintf('\nCall opt_help(''<name>'') for a single entry.\n');
end

function localPrintEntry(e)
    fprintf('opt.%s\n', e.name);
    fprintf('  group   : %s\n', e.group);
    fprintf('  default : %s\n', localFormatVal(e.default));
    fprintf('  help    : %s\n', e.help);
end

function s = localFormatVal(v)
    if islogical(v) || isnumeric(v)
        if isempty(v),       s = '[]';
        elseif isscalar(v),  s = mat2str(v);
        else,                s = sprintf('%s (%s)', mat2str(size(v)), class(v));
        end
    elseif ischar(v)
        s = sprintf('''%s''', v);
    elseif iscell(v)
        s = sprintf('{%s}', strjoin(cellfun(@(c) localFormatVal(c), v(:)', 'uni', false), ', '));
    elseif isstruct(v)
        s = sprintf('struct(%s)', strjoin(fieldnames(v), ', '));
    else
        s = ['<' class(v) '>'];
    end
end
