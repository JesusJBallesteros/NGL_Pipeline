function blocks = sessionBlocks(EventRecord, opt, varargin)
%SESSIONBLOCKS  Recover block / stage boundaries from the session's event record.
%
% PURPOSE:
%   Some events do not belong to a trial: they mark where one stage of a
%   session ends and the next begins - acquisition vs extinction vs test, or
%   one task vs another in a session that runs several. They are read ONCE per
%   session, not per trial, and what they give is a set of time ranges: when to
%   draw a stage boundary on a plot, and which trials belong to which stage.
%
%   That second use matters as soon as a session holds more than one task, as
%   not every analysis applies to every task: a frequency-tagging spectrum is
%   meaningless on a delay-match block, and a correct-vs-incorrect contrast is
%   meaningless on a passive stream. The block table is what an analysis asks
%   before it runs.
%
% USAGE:
%   blocks = sessionBlocks(EventRecord, opt)
%   blocks = sessionBlocks(EventRecord, opt, 'labels', struct('dms', 18))
%   blocks = sessionBlocks(EventRecord, opt, 'onCode', 16, 'offCode', 17)
%
% INPUTS:
%   EventRecord - as written by NGL01 (.EventType, .TimeSecFromMidnight or
%                 .TimeMsFromMidnight).
%   opt         - options; reads opt.eventdef for the marker codes when they
%                 are not given explicitly (fields blockOn / blockOff).
%   Name/value pairs:
%     'onCode'   code that opens a block (default opt.eventdef.blockOn, or 16)
%     'offCode'  code that closes it (default opt.eventdef.blockOff, or 17).
%                A block with no closing code runs to the next block's start,
%                or to the last event of the session.
%     'labels'   struct mapping a block NAME to the code that follows the
%                opening marker and identifies it, e.g.
%                struct('dms', 18, 'arena', 19, 'nft', 20). Without it a block
%                is named by its identifying code, or 'block<N>' when the
%                opening marker carries none.
%     'trialCode' the trial-start code used to count trials per block
%                (default opt.eventdef.itiOn, or 0)
%     'relative' true (default) returns times on the RECORDING clock - the one
%                FT_data.time uses - instead of seconds-from-midnight.
%     't0'       the seconds-from-midnight value of recording sample 1. Only
%                consulted when 'relative' is true. The default is the first
%                event in the record, which is correct for INTAN (its first
%                event is emitted when recording starts) and wrong for any
%                record whose first event comes later - pass t0 explicitly
%                there, or the blocks come back shifted by that gap.
%
% OUTPUT (struct array, one entry per block, in time order):
%   .index      1..nBlocks
%   .label      block name
%   .code       the identifying code, NaN when the marker carried none
%   .tStart     block start, seconds
%   .tEnd       block end, seconds
%   .duration   tEnd - tStart
%   .nTrials    trial-start events inside the block
%   .trialIdx   indices of those trials, counted over the whole session, so
%               blocks(k).trialIdx indexes condition fields and trialdef rows
%
% NOTES:
%   * Trials are assigned to the block their START falls in. A trial that
%     straddles a boundary belongs to the block it began in, which is the only
%     assignment that keeps the trial count consistent.
%   * A session with no block markers returns an empty struct array, not an
%     error: most sessions are one block, and the caller decides whether the
%     markers were required.
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).

    p = inputParser;
    p.addParameter('onCode',    []);
    p.addParameter('offCode',   []);
    p.addParameter('labels',    struct());
    p.addParameter('trialCode', []);
    p.addParameter('relative',  true);
    p.addParameter('t0',        []);
    p.parse(varargin{:});
    a = p.Results;

    onCode    = localCode(a.onCode,    opt, 'blockOn',  16);
    offCode   = localCode(a.offCode,   opt, 'blockOff', 17);
    trialCode = localCode(a.trialCode, opt, 'itiOn',     0);

    codes = EventRecord.EventType(:);
    t = localTimes(EventRecord);
    if a.relative && ~isempty(t)
        if isempty(a.t0), t = t - t(1); else, t = t - a.t0; end
    end

    iOn = find(codes == onCode);
    blocks = struct('index', {}, 'label', {}, 'code', {}, 'tStart', {}, ...
                    'tEnd', {}, 'duration', {}, 'nTrials', {}, 'trialIdx', {});
    if isempty(iOn), return; end

    iOff = find(codes == offCode);
    trialStarts = t(codes == trialCode);
    names = fieldnames(a.labels);

    for k = 1:numel(iOn)
        b = struct();
        b.index  = k;
        b.tStart = t(iOn(k));
        % The identifying code is the next event, when it is not itself a
        % trial or block marker: the opening marker says "a block starts",
        % the code after it says which one.
        b.code = NaN;
        if iOn(k) < numel(codes)
            nxt = codes(iOn(k) + 1);
            if ~ismember(nxt, [onCode, offCode, trialCode])
                b.code = nxt;
            end
        end
        % End: the first closing marker after the start, else the next block's
        % start, else the end of the session.
        e = iOff(find(iOff > iOn(k), 1));
        if ~isempty(e)
            b.tEnd = t(e);
        elseif k < numel(iOn)
            b.tEnd = t(iOn(k + 1));
        else
            b.tEnd = t(end);
        end
        b.duration = b.tEnd - b.tStart;

        b.label = sprintf('block%d', k);
        if ~isnan(b.code)
            b.label = sprintf('code%d', b.code);
            hit = names(cellfun(@(n) isequal(a.labels.(n), b.code), names));
            if ~isempty(hit), b.label = hit{1}; end
        end

        in = trialStarts >= b.tStart & trialStarts < b.tEnd;
        b.trialIdx = find(in(:))';
        b.nTrials  = numel(b.trialIdx);
        blocks(k) = orderfields(b, {'index','label','code','tStart','tEnd', ...
                                    'duration','nTrials','trialIdx'}); %#ok<AGROW>
    end
end

function c = localCode(given, opt, name, dflt)
    if ~isempty(given), c = given; return; end
    if isstruct(opt) && isfield(opt, 'eventdef') && isfield(opt.eventdef, name)
        c = opt.eventdef.(name);
    else
        c = dflt;
    end
end

function t = localTimes(EventRecord)
    if isfield(EventRecord, 'TimeSecFromMidnight')
        t = EventRecord.TimeSecFromMidnight(:);
    elseif isfield(EventRecord, 'TimeMsFromMidnight')
        t = EventRecord.TimeMsFromMidnight(:) / 1000;
    else
        error('sessionBlocks:noTime', ...
            'EventRecord has neither TimeSecFromMidnight nor TimeMsFromMidnight.');
    end
end
