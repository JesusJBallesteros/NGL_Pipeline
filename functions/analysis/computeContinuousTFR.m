function TFR = computeContinuousTFR(FT_data, opt)
% computeContinuousTFR  Multitaper continuous-mode TFR (compute-only).
%
% PURPOSE:
%   Compute half of the pre-26.06.2026 continous_MTspectrogram. Runs
%   ft_freqanalysis with mtmconvol on the FIRST channel of the FT_data
%   struct (matches legacy behaviour; do NOT read this as "all channels
%   averaged"), applies a no-op ft_freqbaseline (baseline off), and
%   RETURNS + SAVES the TFR struct. Plotting is now a separate call
%   (plotContinuousTFR) so re-plotting doesn't re-compute the FFT.
%
% USAGE:
%   TFR = computeContinuousTFR(FT_data, opt);
%
% INPUTS:
%   FT_data - continuous FieldTrip data (the *_FTcont.mat payload).
%   opt     - resolved options struct. Used fields:
%               .analysis       output folder (TFR .mat lives there)
%               .SavFileName    session name (for filename)
%
% OUTPUT (struct):
%   .dB      - ft_freqbaseline-processed freq struct (baseline off, kept
%              for parity with the pre-split code path).
%   .cfg     - the ft_freqanalysis cfg used (so re-plotters can inspect
%              foi, method, tapers, etc.).
%   .provenance - {source_FTfile, source_mtime, opt_snapshot, savedAt}.
%
% FILE:
%   <opt.analysis>/<opt.SavFileName>_TFR_continuous.mat  (variable: TFR)
%
% NOTES:
%   * cfg is currently HARDCODED (foi = 1:1:40, t_ftimwin = 1 s,
%     tapsmofrq = 2, channel = FT_data.label(1)). Pass 2 lifts these
%     into opt.lfp.tfr.cont.* schema entries.
%   * Baseline is off (cfg2.baseline = 'no'). The pre-split code had a
%     commented-out anesthesia baseline; that hook lives at the plotter
%     side now.
%
% SEE ALSO:
%   plotContinuousTFR (draws the PNG from a computed TFR),
%   NGL02_LFP (caller), NGL07_LFPanalysis (planned research-grade stage).
%
% Last modified 26.06.2026 (Jesus) - Pass 1 split of continous_MTspectrogram.

    %% MTM cfg (legacy defaults; Pass 2 lifts to opt.lfp.tfr.cont.*).
    cfg              = [];
    cfg.channel      = FT_data.label(1);
    cfg.method       = 'mtmconvol';
    cfg.output       = 'pow';
    cfg.taper        = 'dpss';
    cfg.foi          = 1:1:40;
    cfg.t_ftimwin    = ones(length(cfg.foi), 1) .* 1;
    cfg.toi          = '50%';
    cfg.tapsmofrq    = 2;
    cfg.keeptrials   = 'no';
    cfg.polyremoval  = 0;
    cfg.pad          = 'nextpow2';

    %% 1 TFR analysis.
    abs = ft_freqanalysis(cfg, FT_data);

    %% 2 Baseline pass (currently off; retained for symmetry).
    cfg2              = [];
    cfg2.baselinetype = 'db';
    cfg2.baseline     = 'no';
    TFR.dB            = ft_freqbaseline(cfg2, abs);

    %% 3 Provenance + persist.
    TFR.cfg = cfg;
    TFR.provenance = buildLFPProvenance(localGuessSourceFT(opt), opt, struct( ...
        'mode',      'continuous', ...
        'channel',   cfg.channel, ...
        'foi_range', [min(cfg.foi), max(cfg.foi)]));

    if ~isfolder(opt.analysis), mkdir(opt.analysis); end
    outFile = fullfile(opt.analysis, [opt.SavFileName '_TFR_continuous.mat']);
    save(outFile, 'TFR', '-mat');
    fprintf('computeContinuousTFR: wrote %s\n', outFile);
end

% -----------------------------------------------------------------------
function p = localGuessSourceFT(opt)
% Best-effort source-file record for provenance. NGL02_LFP loads the
% *_FTcont.mat from opt.FolderProcDataMat; mirror that here without
% requiring FT_data to carry a `.cfg.previous.filename` breadcrumb.
    if isfield(opt, 'FolderProcDataMat') && isfield(opt, 'SavFileName')
        p = fullfile(opt.FolderProcDataMat, [opt.SavFileName '_FTcont.mat']);
    else
        p = '';
    end
end
