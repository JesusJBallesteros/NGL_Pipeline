function input = guardSyntheticSubjects(input, searchRoot, askedForAll, available)
%GUARDSYNTHETICSUBJECTS  Keep synthetic (fixture) data out of real results.
%
% PURPOSE:
%   Synthetic data is built to be analysed, so it looks like real data to every
%   stage. That is what makes it useful and what makes it dangerous: one
%   fixture subject swept into a group mean changes a number that nobody will
%   think to question. This is the one place where the two are separated, and
%   it sits in set_default because that is where every stage resolves its
%   subject list - NGL03 and NGL04 included.
%
% USAGE:
%   input = guardSyntheticSubjects(input, searchRoot, askedForAll, available)
%   Called by set_default; not normally called directly.
%
% INPUTS:
%   input       - the input struct, with .subjects already resolved to a
%                 dir-struct array.
%   searchRoot  - folder the subject list was resolved against.
%   askedForAll - true when the user asked for 'all' rather than naming them.
%   available   - every subject folder found under searchRoot.
%
% OUTPUT:
%   input - with synthetic subjects removed when 'all' was asked for, and
%           .fixtureRun set to true when the run is a synthetic one.
%
% THE RULES:
%   1. 'all' silently skips synthetic subjects.
%   2. Naming them explicitly runs them, loudly.
%   3. Synthetic and real together is refused, as is a synthetic run inside a
%      study that also holds real subjects (study-level files would collide).
%
% SEE ALSO:
%   isSyntheticSubject, tests/fixture/NGL_RunFixture.
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).
    if isempty(input.subjects), return; end

    names = {input.subjects.name};
    isSynth = false(size(names));
    for k = 1:numel(names)
        isSynth(k) = isSyntheticSubject(fullfile(searchRoot, names{k}));
    end

    if askedForAll
        if any(isSynth)
            fprintf(['set_default: skipping synthetic subject(s) %s - ', ...
                     'name them explicitly to analyse the fixture.\n'], ...
                    strjoin(names(isSynth), ', '));
        end
        input.subjects = input.subjects(~isSynth);
        input.fixtureRun = false;
        return
    end

    if any(isSynth) && any(~isSynth)
        error('NGL:mixedSyntheticSubjects', ...
            ['Synthetic subject(s) %s were requested together with real ', ...
             'subject(s) %s. Fixture data must never be aggregated with real ', ...
             'data: run it on its own, in its own study folder.'], ...
            strjoin(names(isSynth), ', '), strjoin(names(~isSynth), ', '));
    end

    input.fixtureRun = any(isSynth);
    if input.fixtureRun
        % A fixture run inside a study that also holds real subjects would
        % write its aggregates over theirs - same study-level file names, same
        % folder. The fixture belongs in its own study root.
        othersReal = setdiff({available.name}, names);
        for k = 1:numel(othersReal)
            if ~isSyntheticSubject(fullfile(searchRoot, othersReal{k}))
                error('NGL:fixtureInRealStudy', ...
                    ['Synthetic subject %s lives in a study that also holds ', ...
                     'real subject(s) (%s). Study-level outputs would collide. ', ...
                     'Give the fixture its own studyname and re-run.'], ...
                    strjoin(names(isSynth), ', '), strjoin(othersReal, ', '));
            end
        end
        fprintf(['\n*** SYNTHETIC DATA RUN: subject(s) %s are generated test ', ...
                 'data, not recordings. Results are for testing the pipeline ', ...
                 'only. ***\n\n'], strjoin(names(isSynth), ', '));
    end
end
