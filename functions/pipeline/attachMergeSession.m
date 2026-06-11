function condition = attachMergeSession(condition, mergedFolder)
% attachMergeSession  Attach per-trial 'A'/'B' tags from a merged session.
%
% PURPOSE:
%   Downstream consumer of mergeMeta.mat. When a session was produced by
%   NGL_mergeSessionsINTAN, this helper reads mergeMeta and adds
%       condition.session   1xNtrials cellstr of 'A' or 'B'
%       condition.mergeMeta the loaded mergeMeta struct (for reference)
%   to an existing condition struct. Idempotent and no-op if mergeMeta.mat
%   is absent (i.e. the session is a normal single recording).
%
% USAGE (from conditions_script or NGL02_postPhy):
%   condition = attachMergeSession(condition, opt.FolderProcDataMat);
%
% INPUTS:
%   condition    - struct with per-trial fields. Must contain at least
%                  one Ntrials-length vector field; the helper uses the
%                  first such field to determine Ntrials.
%   mergedFolder - char, folder where mergeMeta.mat might live (typically
%                  the session's preprocessing folder).
%
% OUTPUT:
%   condition    - same struct, with .session and .mergeMeta added if
%                  mergeMeta.mat was found and the trial counts line up.
%
% BEHAVIOUR:
%   - mergeMeta.mat missing                    -> condition returned unchanged.
%   - mergeMeta.nTrialsA + nTrialsB != Ntrials -> warning + unchanged.
%     This catches the case where the trial table got renumbered after
%     the merge in a way that breaks the boundary assumption.
%   - Already has condition.session            -> no-op (idempotent).
%
% Last modified 09.06.2026 (Jesus)

    if isfield(condition, 'session') && ~isempty(condition.session)
        return  % already tagged
    end

    metaPath = fullfile(mergedFolder, 'mergeMeta.mat');
    if ~isfile(metaPath)
        return  % not a merged session
    end

    S = load(metaPath, 'mergeMeta');
    if ~isfield(S, 'mergeMeta')
        return
    end
    mm = S.mergeMeta;

    Ntrials = localFirstVectorLength(condition);
    if isempty(Ntrials)
        warning('NGL:attachMergeSession:noTrialCount', ...
            'Cannot determine Ntrials from condition struct; skipping merge tagging.');
        return
    end

    expected = mm.nTrialsA + mm.nTrialsB;
    if expected ~= Ntrials
        warning('NGL:attachMergeSession:countMismatch', ...
            ['mergeMeta reports nTrialsA+nTrialsB=%d but condition has %d trials. ', ...
             'Skipping session tagging — boundary may have shifted after merge.'], ...
            expected, Ntrials);
        return
    end

    condition.session   = [repmat({'A'}, mm.nTrialsA, 1); ...
                           repmat({'B'}, mm.nTrialsB, 1)]';
    condition.mergeMeta = mm;
end

function n = localFirstVectorLength(condition)
    n  = [];
    fn = fieldnames(condition);
    for k = 1:numel(fn)
        v = condition.(fn{k});
        if (isnumeric(v) || islogical(v) || iscell(v)) && isvector(v)
            n = numel(v);
            return
        end
    end
end
