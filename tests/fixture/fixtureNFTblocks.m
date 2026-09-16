function blocks = fixtureNFTblocks(root, varargin)
%FIXTURENFTBLOCKS  opt.nft.blocks for the fixture, valid in every session.
%
% PURPOSE:
%   NGL08_NFT takes ONE opt.nft.blocks for the whole run, but each session's
%   tagging block starts at a slightly different second - trials are drawn per
%   session, so everything before the block drifts. A window taken from one
%   session can therefore fall outside another's block and epoch the silence
%   on either side of it.
%
%   This reads every requested session's block table and returns the
%   INTERSECTION: the window each stream occupies in all of them, shrunk by a
%   margin. Narrower than any single session's block, and correct in all.
%
% USAGE:
%   opt.nft.blocks = fixtureNFTblocks('D:\PIPELINE_TEST')
%   opt.nft.blocks = fixtureNFTblocks(root, 'sessions', {'19850214'})
%
% INPUTS:
%   root - the fixture's study folder.
%   Name/value pairs:
%     'subject'  default 'MRX'
%     'sessions' default every session found under the subject
%     'margin'   seconds trimmed from each end (default 1)
%
% OUTPUT:
%   blocks - struct array with .name, .base (Hz) and .window [t0 t1], the
%            shape opt.nft.blocks requires.
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).

    p = inputParser;
    p.addParameter('subject',  'MRX');
    p.addParameter('sessions', {});
    p.addParameter('margin',   1);
    p.parse(varargin{:});
    a = p.Results;

    base = fullfile(root, 'data', 'preprocessing', a.subject);
    assert(isfolder(base), 'fixtureNFTblocks:noSubject', ...
        'no fixture for %s under %s.', a.subject, root);
    sessions = a.sessions;
    if isempty(sessions)
        d = dir(base); d = d([d.isdir] & ~startsWith({d.name}, '.'));
        sessions = {d.name};
    end

    % One row per session, one column per stream: collect first, intersect
    % after, so the loop has nothing to grow.
    t0 = []; t1 = []; rates = [];
    for k = 1:numel(sessions)
        T = load(fullfile(base, sessions{k}, 'fixture_truth.mat'));
        b = T.blocks(strcmp({T.blocks.task}, 'nft'));
        if isempty(b), continue; end
        trials = b(1).trials;
        r = unique([trials.rate]);
        if isempty(rates)
            rates = r;
            t0 = nan(numel(sessions), numel(r));
            t1 = nan(numel(sessions), numel(r));
        end
        assert(isequal(r, rates), 'fixtureNFTblocks:rates', ...
            'session %s tags at different rates than the first.', sessions{k});
        for s = 1:numel(r)
            in = [trials.rate] == r(s);
            t0(k, s) = min([trials(in).tITI]); %#ok<AGROW> preallocated above
            t1(k, s) = max([trials(in).tEnd]); %#ok<AGROW>
        end
    end
    assert(~isempty(rates), 'fixtureNFTblocks:noBlock', ...
        'none of the sessions has a tagging block.');
    lo = max(t0, [], 1, 'omitnan');     % starts last in some session
    hi = min(t1, [], 1, 'omitnan');     % ends first in some session

    blocks = struct('name', {}, 'base', {}, 'window', {});
    for s = 1:numel(rates)
        w = [lo(s) + a.margin, hi(s) - a.margin];
        assert(diff(w) > 5 / rates(s), 'fixtureNFTblocks:tooNarrow', ...
            ['the %g Hz window common to all sessions is only %.1f s long; ', ...
             'reduce ''margin'' or use fewer sessions.'], rates(s), diff(w));
        blocks(s) = struct('name', sprintf('stream%.1fHz', rates(s)), ...
                           'base', rates(s), 'window', w);
    end
end
