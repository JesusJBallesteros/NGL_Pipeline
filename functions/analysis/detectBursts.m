function bursts = detectBursts(FT_data, opt)
% detectBursts  Per-band oscillatory burst detection (threshold-based first pass).
%
% PURPOSE:
%   Identifies oscillatory burst episodes in continuous LFP data. For
%   each band, bandpasses the signal, extracts the analytic envelope
%   (Hilbert), z-scores the envelope in log-space per channel, and
%   emits episodes where the z-scored envelope stays above
%   opt.lfp.burst.threshMult for at least opt.lfp.burst.minDurationMs
%   ms.
%
% METHOD (this is a FIRST PASS - see UPGRADE):
%   1. Bandpass filter with ft_preproc_bandpassfilter (Butterworth,
%      order 4, twopass) per band.
%   2. Analytic signal: hilbert() -> abs() = amplitude envelope.
%   3. Log-transform envelope + z-score PER CHANNEL over the whole
%      trace (log-normal is a decent stationary-Gaussian approximation
%      for envelope amplitudes and stabilizes the threshold across
%      1/f-slope differences).
%   4. Threshold at opt.lfp.burst.threshMult standard deviations.
%      Contiguous supra-threshold runs are candidate episodes.
%   5. Reject episodes shorter than opt.lfp.burst.minDurationMs.
%   6. For each surviving episode record start / end sample indices,
%      peak envelope, peak time, peak-instant phase (from the analytic
%      signal), and the channel.
%
% UPGRADE PATH:
%   This detector uses a fixed z-threshold; it does NOT separate
%   aperiodic 1/f activity from true rhythmic bursts. For studies where
%   that matters, install eBOSC (extended Better OSCillation, github/
%   BOSCbase/eBOSC) under toolboxes/eBOSC/ and re-implement this
%   function as a wrapper over BOSC's robust 1/f fit + power-time-of-
%   interest thresholding. Signature can stay the same.
%
% USAGE:
%   bursts = detectBursts(FT_data, opt);
%
% INPUTS:
%   FT_data - FieldTrip continuous struct with .label, .trial{1},
%             .time{1}, .fsample. Should carry .chanArea (backfilled by
%             ensureChanArea in NGL02_LFP / NGL07 load steps).
%   opt     - resolved options. Consumes:
%               .lfp.bands              cell of {name, [flo fhi]} rows
%               .lfp.burst.threshMult   z-sigma threshold (default 3)
%               .lfp.burst.minDurationMs (default 100)
%               .lfp.tfrAreaFilter      (optional) restrict to a subset
%                                        of channels by chanArea
%
% OUTPUT (struct):
%   .bands       - {nBands x 1} cell of band names.
%   .bandRange   - [nBands x 2] double, [flo fhi] per band.
%   .episodes    - table with columns:
%                    band          char (band name)
%                    channel       char (FT_data.label{i})
%                    chanArea      char (FT_data.chanArea{i} or 'main')
%                    startSample   int
%                    endSample     int
%                    startSec      double (s from trial start)
%                    endSec        double
%                    peakSample    int
%                    peakSec       double
%                    peakEnvAbs    double (raw envelope amplitude)
%                    peakEnvZ      double (z-scored envelope)
%                    peakPhase     double (radians, [-pi pi])
%                    durationMs    double
%   .params      - {threshMult, minDurationMs, filter order, source}
%   .provenance  - buildLFPProvenance snapshot.
%
% SEE ALSO:
%   hilbertBandpass (returns the full per-band phase / envelope traces
%       rather than discrete episodes; use it when you need continuous
%       signals for e.g. spike-phase relationships), NGL07_LFPanalysis.
%
% Last modified 26.06.2026 (Jesus) - new (LFP Pass 3).

    threshMult    = localOptField(opt, {'lfp','burst','threshMult'},    3);
    minDurationMs = localOptField(opt, {'lfp','burst','minDurationMs'}, 100);
    bands         = localOptField(opt, {'lfp','bands'}, ...
                        {{'theta',[4 8]}, {'beta',[15 30]}, {'gamma',[30 90]}});

    %% Optional per-area subselection.
    keepChanIdx = 1:numel(FT_data.label);
    if isfield(opt,'lfp') && isfield(opt.lfp,'tfrAreaFilter') ...
            && ~isempty(opt.lfp.tfrAreaFilter) && isfield(FT_data,'chanArea')
        keepAreas   = cellstr(opt.lfp.tfrAreaFilter);
        keepChanIdx = find(ismember(FT_data.chanArea, keepAreas));
        if isempty(keepChanIdx)
            error('NGL:detectBursts:noChans', ...
                'opt.lfp.tfrAreaFilter matched 0 channels.');
        end
    end

    x    = FT_data.trial{1}(keepChanIdx, :);   % [nChan x nSamples]
    fs   = FT_data.fsample;
    nCh  = size(x, 1);
    nS   = size(x, 2);
    minS = round(minDurationMs * fs / 1000);

    nBands    = numel(bands);
    bandName  = cell(nBands, 1);
    bandRange = zeros(nBands, 2);
    epCollect = struct('band',{}, 'channel',{}, 'chanArea',{}, ...
                       'startSample',{}, 'endSample',{}, ...
                       'startSec',{}, 'endSec',{}, ...
                       'peakSample',{}, 'peakSec',{}, ...
                       'peakEnvAbs',{}, 'peakEnvZ',{}, 'peakPhase',{}, ...
                       'durationMs',{});

    for b = 1:nBands
        bandName{b}     = bands{b}{1};
        bandRange(b, :) = bands{b}{2};
        flo = bandRange(b, 1);  fhi = bandRange(b, 2);

        for c = 1:nCh
            % Filter + Hilbert.
            xf   = ft_preproc_bandpassfilter(x(c, :), fs, [flo fhi], 4, 'but', 'twopass');
            ana  = hilbert(xf);
            env  = abs(ana);
            phi  = angle(ana);

            % Log-normal-friendly z-score per channel.
            logE  = log(env + eps);
            zLogE = (logE - mean(logE, 'omitnan')) ./ std(logE, 0, 'omitnan');

            % Threshold + segment.
            above = zLogE > threshMult;
            if ~any(above), continue, end
            [starts, ends] = localRunEndpoints(above);
            dur = ends - starts + 1;
            keep = dur >= minS;
            starts = starts(keep);  ends = ends(keep);

            for e = 1:numel(starts)
                seg           = starts(e):ends(e);
                [pkZ, pkOff]  = max(zLogE(seg));
                pk            = seg(pkOff);
                row           = struct();
                row.band       = bandName{b};
                row.channel    = FT_data.label{keepChanIdx(c)};
                if isfield(FT_data,'chanArea')
                    row.chanArea = FT_data.chanArea{keepChanIdx(c)};
                else
                    row.chanArea = 'main';
                end
                row.startSample = starts(e);
                row.endSample   = ends(e);
                row.startSec    = (starts(e) - 1) / fs;
                row.endSec      = (ends(e)   - 1) / fs;
                row.peakSample  = pk;
                row.peakSec     = (pk - 1) / fs;
                row.peakEnvAbs  = env(pk);
                row.peakEnvZ    = pkZ;
                row.peakPhase   = phi(pk);
                row.durationMs  = (ends(e) - starts(e) + 1) * 1000 / fs;
                epCollect(end + 1) = row; %#ok<AGROW>
            end
        end
        fprintf('detectBursts: %s [%g-%g Hz]: %d episodes (%d chans).\n', ...
                bandName{b}, flo, fhi, ...
                sum(strcmp({epCollect.band}, bandName{b})), nCh);
    end

    bursts.bands      = bandName;
    bursts.bandRange  = bandRange;
    if isempty(epCollect)
        bursts.episodes = table();
    else
        bursts.episodes = struct2table(epCollect);
    end
    bursts.params     = struct( ...
        'threshMult',    threshMult, ...
        'minDurationMs', minDurationMs, ...
        'filterOrder',   4, ...
        'filterType',    'butter twopass', ...
        'zMethod',       'log envelope, per-channel');
    bursts.provenance = buildLFPProvenance('', opt, struct( ...
        'analysis', 'detectBursts', ...
        'nBands',   nBands, ...
        'nChans',   nCh));
end


% =======================================================================
function [starts, ends] = localRunEndpoints(mask)
% Contiguous true-runs in a logical vector -> start/end sample indices.
    m = mask(:)';
    d = diff([false, m, false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
end

function v = localOptField(opt, path, dflt)
% Read opt.(path{1}).(path{2})... with default. Missing intermediate -> dflt.
    v = dflt;
    cursor = opt;
    for k = 1:numel(path)
        if isstruct(cursor) && isfield(cursor, path{k})
            cursor = cursor.(path{k});
        else
            return
        end
    end
    v = cursor;
end
