function parsed = parseFireRateRequest(request, catSets)
% parseFireRateRequest  Parse a 3-cell `request` into a parsed struct.
%
% PURPOSE:
%   Convert the user-facing 3-cell `request` used by NGL04_fireRate
%   (and NGL04_PCA) into a struct with explicit per-factor level lists.
%   Categories (alignment / condition / label) are inferred from token
%   content using `catSets` (built by buildRequestCatSets), so the
%   order of the three entries is flexible. Each entry may contain a
%   single ' vs ' to define a comparison; both sides of the comparison
%   must resolve to the same category.
%
% USAGE:
%   parsed = parseFireRateRequest(request, catSets);
%
% INPUT:
%   request  - 1x3 cell. Each entry one of:
%                * single value  (e.g. 'correct', 'good', 'stim2')
%                * comparison    (e.g. 'correct vs incorrect')
%   catSets  - struct from buildRequestCatSets with:
%                .alignments     cell of alignment names
%                .conditions     cell of condition fieldnames
%                .labelPool      struct keyed by label-field name
%                .labelPriority  resolution order for label tokens
%
% OUTPUT:
%   parsed   - struct with:
%                .alignment  {1x1 or 1x2} cell of alignment names
%                .condition  {1x1 or 1x2} cell of condition field names
%                .label      {1x1 or 1x2} cell of label values
%                .labelField {1x1 or 1x2} cell of label-field names
%                            (resolved per label entry via labelPriority)
%                .varying    cell of factor names with >1 level
%                            (e.g. {'alignment','condition'})
%
% Last modified 09.06.2026 (Jesus)

    assert(iscell(request) && numel(request) == 3, ...
        'NGL:parseFireRateRequest:badRequest', ...
        'request must be a 3-element cell.');

    entries     = cell(3, 1);
    entryCats   = cell(3, 1);
    entryFields = cell(3, 1);
    for k = 1:3
        s     = strtrim(request{k});
        parts = regexp(s, '\s+vs\s+', 'split', 'ignorecase');
        parts = cellfun(@strtrim, parts, 'uni', false);
        entries{k}     = parts;
        entryCats{k}   = cell(1, numel(parts));
        entryFields{k} = cell(1, numel(parts));
        for p = 1:numel(parts)
            [cat, field] = localCategoriseToken(parts{p}, catSets);
            entryCats{k}{p}   = cat;
            entryFields{k}{p} = field;
        end
        if ~all(strcmp(entryCats{k}, entryCats{k}{1}))
            error('NGL:parseFireRateRequest:badRequest', ...
                'Entry ''%s'' mixes categories (%s); both sides of ''vs'' must be the same kind.', ...
                request{k}, strjoin(entryCats{k}, ', '));
        end
    end

    byCat       = struct('alignment', {{}}, 'condition', {{}}, 'label', {{}});
    labelFields = {};
    for k = 1:3
        cat = entryCats{k}{1};
        if ~isempty(byCat.(cat))
            error('NGL:parseFireRateRequest:badRequest', ...
                ['Two entries both resolve to category ''%s'': ''%s'' and ''%s''. ', ...
                 'Each of the three slots must be a different category.'], ...
                cat, strjoin(byCat.(cat), ' vs '), strjoin(entries{k}, ' vs '));
        end
        byCat.(cat) = entries{k};
        if strcmp(cat, 'label'), labelFields = entryFields{k}; end
    end

    for c = {'alignment','condition','label'}
        if isempty(byCat.(c{1}))
            error('NGL:parseFireRateRequest:badRequest', ...
                'No entry in request resolves to category ''%s''.', c{1});
        end
    end

    parsed = struct();
    parsed.alignment  = byCat.alignment;
    parsed.condition  = byCat.condition;
    parsed.label      = byCat.label;
    parsed.labelField = labelFields;
    parsed.varying    = {};
    catNames = {'alignment','condition','label'};
    for c = catNames
        if numel(byCat.(c{1})) > 1
            parsed.varying{end+1} = c{1};
        end
    end
end

function [cat, field] = localCategoriseToken(token, catSets)
% Priority for ambiguous tokens: alignment > condition > label, then
% labelPool inspected in catSets.labelPriority order.
    field = '';
    if any(strcmp(token, catSets.alignments))
        cat = 'alignment'; return
    end
    if any(strcmp(token, catSets.conditions))
        cat = 'condition'; return
    end
    for f = catSets.labelPriority
        fld = f{1};
        if isfield(catSets.labelPool, fld) && ...
                any(strcmp(token, catSets.labelPool.(fld)))
            cat   = 'label';
            field = fld;
            return
        end
    end
    error('NGL:parseFireRateRequest:unknownToken', ...
        ['Cannot categorise ''%s''. It is not in opt.alignto, not a ', ...
         'condition fieldname, and not a known cluster label value ', ...
         '(checked in priority %s).'], token, ...
        strjoin(catSets.labelPriority, ' > '));
end
