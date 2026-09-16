function FT_data = makeFixtureFT(lfp, geom, cfg, trialdef, alignName)
%MAKEFIXTUREFT  Package fixture LFP as FieldTrip data, continuous or per trial.
%
% USAGE:
%   FT = makeFixtureFT(lfp, geom, cfg)                       % continuous
%   FT = makeFixtureFT(lfp, geom, cfg, trialdef, 'stimOn1')  % trial-parsed
%
% INPUTS:
%   lfp       - [nChan x nSamples] continuous signal at cfg.fs
%   geom      - channel geometry (label, area)
%   cfg       - fixture configuration (fs)
%   trialdef  - the 2 x nAlign cell array, windows in ms
%   alignName - which alignment to cut
%
% OUTPUT:
%   FT_data - FieldTrip raw structure, with .chanArea as NGL01 writes it.
%
% NOTES:
%   A trial is dropped when its window falls outside the recording, or when
%   this alignment's event does not exist in it (a NaN offset - the tagging
%   stream has no stimOn2). Both happen in real sessions, so the fixture
%   reproduces them; .trialinfo then holds the ORIGINAL trial index of every
%   trial kept, which is how a per-trial condition vector is matched to the
%   trials that actually survived.
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).

    n = size(lfp, 2);
    if nargin < 4 || isempty(trialdef)
        FT_data = struct('label', {geom.label}, 'fsample', cfg.fs, ...
                         'trial', {{lfp}}, 'time', {{(0:n-1) / cfg.fs}}, ...
                         'chanArea', {geom.area}, 'sampleinfo', [1 n], ...
                         'cfg', struct('continuous', 'yes', 'fixture', true));
        return
    end

    col = find(strcmp(trialdef(1, :), alignName), 1);
    assert(~isempty(col), 'makeFixtureFT:align', ...
        'alignment ''%s'' is not in trialdef.', alignName);
    win = trialdef{2, col} / 1000;                      % ms -> s
    nTr = size(win, 1);
    trial = cell(1, nTr); time = cell(1, nTr);
    keep = false(nTr, 1);
    for k = 1:nTr
        if any(~isfinite(win(k, :))), continue; end     % no such event here
        i0 = round(win(k, 1) * cfg.fs) + 1;
        i1 = round(win(k, 2) * cfg.fs);
        if i0 < 1 || i1 > n || i1 <= i0, continue; end
        trial{k} = double(lfp(:, i0:i1));
        time{k}  = ((i0:i1) - 1) / cfg.fs - win(k, 3);
        keep(k)  = true;
    end
    FT_data = struct('label', {geom.label}, 'fsample', cfg.fs, ...
                     'trial', {trial(keep)}, 'time', {time(keep)}, ...
                     'chanArea', {geom.area}, 'trialinfo', find(keep), ...
                     'cfg', struct('fixture', true, 'alignment', alignName, ...
                                   'droppedTrials', sum(~keep)));
end
