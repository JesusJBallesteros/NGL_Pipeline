function bp = hilbertBandpass(FT_data, opt)
% hilbertBandpass  Per-band bandpass + Hilbert -> phase + envelope traces.
%
% PURPOSE:
%   Return continuous instantaneous phase and amplitude envelope per
%   band for every channel. Used by spike-phase relationships, phase-
%   amplitude coupling, event-triggered analyses. Complements
%   detectBursts (episode-level) - use hilbertBandpass when you need
%   the continuous time series.
%
% USAGE:
%   bp = hilbertBandpass(FT_data, opt);
%
% INPUTS:
%   FT_data - FieldTrip continuous struct with .label, .trial{1},
%             .time{1}, .fsample.
%   opt     - resolved options. Consumes:
%               .lfp.bands             cell of {name, [flo fhi]} rows
%                                       (default: theta/beta/gamma)
%               .lfp.tfrAreaFilter     (optional) restrict channels
%               .lfp.hilbert.storePhase   default true
%               .lfp.hilbert.storeEnv     default true
%               .lfp.hilbert.dtype        'single'|'double' (default single)
%
% OUTPUT (struct):
%   .bands       {nBands x 1} cell of band names
%   .bandRange   [nBands x 2] [flo fhi]
%   .channels    {nCh x 1} cell of channel labels used
%   .chanArea    {nCh x 1} cell of area labels
%   .fsample     sample rate (Hz)
%   .time        [1 x nSamples] time vector (s)
%   .phase       {nBands x 1} of [nCh x nSamples] instantaneous phase
%                 (radians, [-pi pi]); only when opt.lfp.hilbert.storePhase
%   .envelope    {nBands x 1} of [nCh x nSamples] amplitude envelope
%                 (uV); only when opt.lfp.hilbert.storeEnv
%   .provenance  buildLFPProvenance snapshot
%
% NOTES:
%   * Phase / envelope stored as single by default to keep session-size
%     .mat files bounded (a 30-min session at 1 kHz with 32 channels x
%     3 bands x double = ~1.4 GB; single halves that).
%   * ft_preproc_bandpassfilter runs Butterworth order 4, 'twopass' to
%     preserve phase. Matches detectBursts for consistency.
%
% SEE ALSO:
%   detectBursts, spikeFieldCoupling, NGL07_LFPanalysis.
%
% Last modified 26.06.2026 (Jesus) - new (LFP Pass 3).

    bands   = localOptField(opt, {'lfp','bands'}, ...
                {{'theta',[4 8]}, {'beta',[15 30]}, {'gamma',[30 90]}});
    doPhase = localOptField(opt, {'lfp','hilbert','storePhase'}, true);
    doEnv   = localOptField(opt, {'lfp','hilbert','storeEnv'},   true);
    dtype   = localOptField(opt, {'lfp','hilbert','dtype'},      'single');

    %% Optional per-area subselection (same idiom as detectBursts).
    keepChanIdx = 1:numel(FT_data.label);
    if isfield(opt,'lfp') && isfield(opt.lfp,'tfrAreaFilter') ...
            && ~isempty(opt.lfp.tfrAreaFilter) && isfield(FT_data,'chanArea')
        keepChanIdx = find(ismember(FT_data.chanArea, cellstr(opt.lfp.tfrAreaFilter)));
        if isempty(keepChanIdx)
            error('NGL:hilbertBandpass:noChans', ...
                'opt.lfp.tfrAreaFilter matched 0 channels.');
        end
    end

    x    = FT_data.trial{1}(keepChanIdx, :);
    fs   = FT_data.fsample;
    nCh  = size(x, 1);
    nS   = size(x, 2);
    nBands = numel(bands);

    bp.bands     = cellfun(@(b) b{1},  bands, 'uni', false)';
    bp.bandRange = cell2mat(cellfun(@(b) b{2}, bands, 'uni', false)');
    bp.channels  = FT_data.label(keepChanIdx);
    if isfield(FT_data,'chanArea')
        bp.chanArea = FT_data.chanArea(keepChanIdx);
    else
        bp.chanArea = repmat({'main'}, nCh, 1);
    end
    bp.fsample = fs;
    bp.time    = FT_data.time{1};
    if doPhase, bp.phase    = cell(nBands, 1); end
    if doEnv,   bp.envelope = cell(nBands, 1); end

    for b = 1:nBands
        flo = bp.bandRange(b, 1);  fhi = bp.bandRange(b, 2);
        phaseAcc = zeros(nCh, nS, dtype);
        envAcc   = zeros(nCh, nS, dtype);
        for c = 1:nCh
            xf  = ft_preproc_bandpassfilter(x(c, :), fs, [flo fhi], 4, 'but', 'twopass');
            ana = hilbert(xf);
            if doPhase, phaseAcc(c, :) = cast(angle(ana), dtype); end
            if doEnv,   envAcc(c,   :) = cast(abs(ana),   dtype); end
        end
        if doPhase, bp.phase{b}    = phaseAcc; end
        if doEnv,   bp.envelope{b} = envAcc;   end
        fprintf('hilbertBandpass: %s [%g-%g Hz] done (%d chans, dtype=%s).\n', ...
                bp.bands{b}, flo, fhi, nCh, dtype);
    end

    bp.provenance = buildLFPProvenance('', opt, struct( ...
        'analysis',   'hilbertBandpass', ...
        'nBands',     nBands, ...
        'nChans',     nCh, ...
        'nSamples',   nS, ...
        'dtype',      dtype, ...
        'storePhase', doPhase, ...
        'storeEnv',   doEnv));
end


% =======================================================================
function v = localOptField(opt, path, dflt)
    v = dflt;  cursor = opt;
    for k = 1:numel(path)
        if isstruct(cursor) && isfield(cursor, path{k})
            cursor = cursor.(path{k});
        else
            return
        end
    end
    v = cursor;
end
