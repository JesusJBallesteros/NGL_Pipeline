function spec = parseTrialContrast(request, condition, nTrials)
% parseTrialContrast  Turn 'correct vs error' into the two trial masks it means.
%
% PURPOSE:
%   Contrasts are written the way they are spoken - 'correct vs error',
%   'stim2 vs ~stim2' - and resolved here against the session's condition
%   struct. One parser, so a contrast means the same thing in a figure title,
%   a file name and a statistical test.
%
% USAGE:
%   spec = parseTrialContrast('correct vs error', condition, nTrials)
%   spec = parseTrialContrast('correct vs ~correct', condition, nTrials)
%   spec = parseTrialContrast(struct('A', maskA, 'B', maskB, ...
%                                    'labelA','early','labelB','late'), [], nTrials)
%
% INPUTS:
%   request   - char 'A vs B', where each side is a field of `condition`,
%               optionally negated with '~', or the word 'all'. Also accepts a
%               struct with .A / .B logical masks (and optional .labelA /
%               .labelB), for contrasts a condition field cannot express.
%   condition - condition struct from NGL02_postPhy; each field a per-trial
%               vector. A trial belongs to a side when the field is non-zero,
%               matching the quick-look TFR's convention.
%   nTrials   - trials in the TFR the masks will index.
%
% OUTPUT (struct):
%   .A, .B        logical [nTrials x 1] masks
%   .labelA/.labelB, .label ('A-vs-B'), .request (as given)
%   .nA, .nB      trial counts
%   .overlap      trials in both sides (see NOTES)
%
% NOTES:
%   * Condition vectors shorter or longer than the TFR's trial count are a
%     real hazard: trials get dropped by artifact rejection after the
%     condition struct is written. A mismatch raises rather than silently
%     truncating, because silently comparing the wrong trials is worse than
%     stopping.
%   * Overlapping sides are reported, not fixed. The 'trials' design treats
%     the two groups as independent samples, which a shared trial violates;
%     whether that matters is the analyst's call, but it should not be
%     invisible.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 1).

    if isstruct(request)
        spec = localFromStruct(request, nTrials);
        return
    end

    assert(ischar(request) || isstring(request), 'parseTrialContrast:request', ...
        'request must be ''A vs B'' or a struct with .A / .B masks.');
    txt = char(request);
    parts = regexp(txt, '\s+vs\.?\s+', 'split', 'ignorecase');
    assert(numel(parts) == 2, 'parseTrialContrast:format', ...
        ['request must read ''A vs B'' (got ''%s''). Each side is a condition ', ...
         'field, ''~field'' for its complement, or ''all''.'], txt);

    [spec.A, spec.labelA] = localSide(strtrim(parts{1}), condition, nTrials);
    [spec.B, spec.labelB] = localSide(strtrim(parts{2}), condition, nTrials);
    spec.request = txt;
    spec.label   = sprintf('%s-vs-%s', spec.labelA, spec.labelB);
    spec = localCounts(spec);
end

% ---------------- helpers ----------------
function [mask, label] = localSide(side, condition, nTrials)
    negate = startsWith(side, '~') || startsWith(side, '!');
    field  = strtrim(erase(side, {'~', '!'}));

    if strcmpi(field, 'all')
        mask = true(nTrials, 1);
    else
        assert(isstruct(condition) && isfield(condition, field), ...
            'parseTrialContrast:noField', ...
            'condition.%s does not exist. Available: %s', field, ...
            localAvailable(condition));
        v = condition.(field);
        v = v(:);
        assert(numel(v) == nTrials, 'parseTrialContrast:length', ...
            ['condition.%s has %d entries but the TFR has %d trials. They must ', ...
             'be the same trials in the same order - artifact rejection after ', ...
             'the condition struct was written is the usual cause.'], ...
            field, numel(v), nTrials);
        mask = v ~= 0 & ~isnan(v);
    end
    if negate
        mask = ~mask;
        label = ['not-' field];
    else
        label = field;
    end
end

function spec = localFromStruct(request, nTrials)
    assert(isfield(request, 'A') && isfield(request, 'B'), ...
        'parseTrialContrast:structFields', ...
        'a struct request needs .A and .B logical masks.');
    spec.A = localMask(request.A, nTrials, 'A');
    spec.B = localMask(request.B, nTrials, 'B');
    spec.labelA = localField(request, 'labelA', 'A');
    spec.labelB = localField(request, 'labelB', 'B');
    spec.request = localField(request, 'request', ...
                              sprintf('%s vs %s', spec.labelA, spec.labelB));
    spec.label = sprintf('%s-vs-%s', spec.labelA, spec.labelB);
    spec = localCounts(spec);
end

function m = localMask(m, nTrials, name)
    m = logical(m(:));
    assert(numel(m) == nTrials, 'parseTrialContrast:maskLength', ...
        'mask %s has %d entries but the TFR has %d trials.', name, numel(m), nTrials);
end

function spec = localCounts(spec)
    spec.nA = sum(spec.A);
    spec.nB = sum(spec.B);
    spec.overlap = sum(spec.A & spec.B);
    assert(spec.nA > 0 && spec.nB > 0, 'parseTrialContrast:emptySide', ...
        ['contrast ''%s'' leaves one side empty (%d vs %d trials); nothing to ', ...
         'compare.'], spec.request, spec.nA, spec.nB);
    if spec.overlap > 0
        warning('parseTrialContrast:overlap', ...
            ['contrast ''%s'': %d trial(s) belong to both sides. An independent-', ...
             'samples test assumes they do not overlap.'], spec.request, spec.overlap);
    end
end

function v = localField(s, f, default)
    if isfield(s, f) && ~isempty(s.(f)), v = char(s.(f)); else, v = default; end
end

function txt = localAvailable(condition)
    if ~isstruct(condition), txt = '(no condition struct loaded)'; return; end
    f = fieldnames(condition);
    if numel(f) > 12, f = [f(1:12); {'...'}]; end
    txt = strjoin(f', ', ');
end
