%% NGL_SetAndRunMe. User Configuration and Run Script (TEMPLATE)
%
% PURPOSE:
%   The single file a user must edit to drive the entire NGL pipeline.
%   Copy this template from the toolbox root into your project's
%   analysisCode\ folder and edit it there. Never run it from the toolbox
%   root itself.
%
% USAGE:
%   Run section-by-section (F9) or all at once (F5) from MATLAB.
%   All subsequent pipeline scripts (NGL01_Main, NGL02_postPhy, etc.) are
%   called from within this file — they should never be opened directly.
%
% SECTIONS:
%   1) PREPARE   - study metadata, data drive + project name, run NGL00_Prep.
%   2) SET       - choose subjects/sessions and configure opt / param:
%                    A)  subjects + dates
%                    A2) multi-area (optional)
%                    B)  opt struct (grouped B.1..B.15; commented lines show
%                        the schema default, uncomment to override)
%                    C)  param struct (analysis/plot tuning, semi-independent)
%                  SET 2 (optional) - merge two INTAN sessions
%                        (NGL_mergeSessionsINTAN); run on its own.
%   3) RUN       - call the pipeline stages in order (comment out any to skip):
%                    1  NGL01_Main          preprocessing + spike sorting
%                    2  NGL02_postPhy       spike side (requires Phy curation)
%                    3  NGL02_LFP           LFP quick-look (independent of Phy)
%                    4  NGL03_aggregate     cross-session / cross-subject
%                    5  NGL04_fireRate      cross-subject PSTH plots
%                    6  NGL04_PCA           cross-subject population PCA
%                    7  NGL05_plotting      group-level plots (TODO)
%                    8  NGL06_videoAnalysis DLC gaze pipeline
%                    9  NGL07_LFPanalysis   research-grade session-level LFP
%
% REQUIRED CONFIG FILES (place in analysisCode\ alongside this script):
%   NGL_machineConfig.m   - machine-specific Python env and toolbox paths
%   eventDefinitions.m    - project event code definitions
%   conditions_script.m   - condition grouping logic
%   chanMapXXX.mat        - custom Kilosort channel map (if applicable)
%
% NOTES:
%   - Every recognised option appears below, either UNCOMMENTED (the ones
%     most users override) or COMMENTED-OUT showing the canonical default.
%     Uncomment to override; a commented line falls back to its schema default.
%   - Options are SCHEMA-DRIVEN. The single source of truth is
%     functions/config/optSchema.m — each option is one optEntry(name,
%     default, validator, group, help) line. default_opt() builds the defaults
%     from it (generateDefaultsFromSchema); set_default() validates every field
%     against it (validateOptAgainstSchema + optPostChecks for cross-field
%     rules). To ADD an option, append one optEntry line to optSchema.m — that
%     is the only file to touch (default_opt.m / set_default.m are now thin
%     wrappers over the schema, not hand-maintained lists).
%   - NGL_machineConfig.m is machine-specific; do not commit it to git.
%
% LANDMINE: REDUNDANT FIELDS — DO NOT SET DIFFERENT VALUES FOR THE SAME THING.
%   * opt.binSize_ms / opt.stepSz_ms       canonical FR binning. ALL per-session
%                                          analyses (calculate_fireRate_general,
%                                          calculate_neural_pca, etc.) read
%                                          these. Override here if you want a
%                                          project-wide non-default.
%   * opt.fireRatePlot.binSize_ms / .stepSz_ms / .interval
%                                          ONLY used by NGL04_fireRate's
%                                          plotPSTH call. Independent from the
%                                          per-session bins above — change only
%                                          if you want the cross-subject PSTH
%                                          to use a different binning than
%                                          calcFireRate did per session.
%   * opt.pcaPlot.binSize_ms / .stepSz_ms / .interval
%                                          ONLY used by NGL04_PCA. Same logic.
%   * opt.popDyn.smoothSigma + opt.pcaPlot.smoothSigma_s
%                                          NOT redundant: the former drives
%                                          per-session calculate_neural_pca
%                                          (and jPCA/GPFA placeholders); the
%                                          latter drives NGL04_PCA cross-
%                                          subject. Keep them aligned unless
%                                          you have a reason to differ.
%
% Last modified 24.07.2026 

%% 1) PREPARE.
clear all

% A) README.TXT (written into <studyname>\readme.txt on first run).
readmecontent = ["Study name: DefaultName", ...
                 "Readme date: YYYY/MM/DD"                          , ...
                 "Person (1) responsible for data repository: Main researcher", ...
                 "Person(s) responsible for study: Main R and associates"    , ...
                 "Hardware used: ACQ SYS. Probe type. Arena/Box"                  , ...
                 "Related Publication(s): None."                    , ...
                 "Short description of study: A more elaborated description of project paradigm, goals, etc ." ];

% B) DATA LOCATION.
% ADD the toolbox folder to MATLAB folder system before running.
datadrive   = 'D';            % The LETTER of the drive where the data structure is/will be created.
studyname   = 'ProjectName';  % Name of the study to be used (main folder for the data).

% Creates the folder structure on first run; no-op on subsequent runs.
NGL00_Prep

%% 2) SET.
% A) SUBJECTS AND SESSIONS
subjects    = {'XXX'};                  % char 'all', or cell of subject IDs e.g. {'042'}, {'DOE','JNE'}.
dates       = {'YYYYMMDD','YYYYMMDD'};  % char 'all', or cell of session dates e.g. {'20260101'}.

% A2) MULTI-AREA (optional).
% Uncomment if your probe spans more than one brain area. One label per
% kcoords group in the chanMap (kcoords==1 -> Areas{1}, etc.). Repeated
% labels indicate shanks from the same area (processed together).
% Kilosort runs once per unique area; results go to
% preprocessing\<session>\<Area>\. Comment out for single-area behaviour.

% areas = {'NCL','NCL','STR'};   % e.g. 2 NCL shanks + 1 STR shank.

% B) OPTIONS (opt struct).
%    Fields not set here fall back to defaults in default_opt.m;
%    set_default.m validates every field. Commented entries below show
%    the default, uncomment ONLY to override.
opt = struct();

    % B.1  ACQUISITION
    opt.numChannels             = 32;          % electrode channel count.
    opt.bin                     = true;        % create Kilosort .bin file.
    opt.FieldTrip               = true;        % produce FieldTrip LFP .mat.
    opt.doNWB                   = true;        % NWB export (NeuroConv / matNWB).

    % B.1b LOCAL-PC PostPhy
    opt.regenFrom.preproc       = false;
    opt.regenFrom.system        = 'INTAN';  % or 'Deuteron'
    % When moving your preprocessed data to your desk PC, no RAW would in
    % principle need to move. This will cause problems in an standard run.
    % Making .preproc = true will flag the run for a 'regeneration' of
    % 'info' fields, based on the system of your choice. NGL01_Main then
    % forces kilosort / bombcell / phy / doNWB OFF and only re-runs
    % EventProcess to rebuild events.mat / trialdef.mat / condition.mat
    % from the EventRecord.mat already on disk.

    % B.2 EVENTS 
    opt.RetrieveEvents          = true;        % extract event log.
    opt.alignto                 = {'itiOn'};   % alignment events (cell of char).
    opt.addtime                 = 500;         % Adds (ms) around trial start/end. Gives a BL for itiOn, avoids edge artifacts in several calculations
    opt.trEvents                = {};          % inter-trial events (block changes, treatments).
    % opt.uselog                = false;       % Deuteron text-log fallback.
    % opt.timebreak             = false;       % expect a recording break (Deuteron battery change).

    % B.3 MOTION SENSORS 
    opt.GetMotionSensors        = false;       % head-direction sensor data.

    % B.4 PREPROCESSING FILTERS 
    % opt.lowpass               = 9000;        % .bin low-pass, Hz
    % opt.lowpassFT             = 300;         % FieldTrip LFP low-pass, Hz.
    % opt.linefilter            = 0;           % line-noise notch centre, Hz. 0 = off.
    % opt.noise                 = [];          % reserved for noise-rejection params.

    % B.5 SORTING & CURATION 
    opt.kilosort                = true;        % run Kilosort 4.
    opt.KSchanMapFile           = '';          % '' = linear array; or e.g. 'chanMap_ATLAS_E32-...mat'.
    opt.bombcell                = true;        % run Bombcell QC after sorting.
    opt.callBcGUI               = false;       % open Bombcell GUI after QC
    % opt.phy                   = false;       % open Phy right after QC and GUI (BLOCKS MATLAB).

    % B.6 NGL02 STAGE GATES 
    opt.doSpikething            = true;        % run NGL02_postPhy spike work.
    opt.doLFPthing              = true;        % run NGL02_LFP work.

    % B.7 SPIKE SIDE 
    % B.7.b  Waveform extraction (loadSpikes)
    opt.getwF                   = true;        % extract raw waveforms per cluster (useful).
    opt.gwfparams.nWf           = 1000;        % max waveforms per cluster.
    % opt.gwfparams.wfWin       = [-20 41];    % samples around spiketime.

    % B.7.c Cluster loading + ISI binning
    % opt.spparams.excludeNoise = true;        % skip Phy 'noise' clusters.
    % opt.spparams.loadPCs      = false;       % load principal components.
    % opt.isibins               = 0:0.5:200;   % ISI histogram bin edges, ms.

    % B.7.d Firing-rate binning (CANONICAL — overrides for the whole pipeline)
    opt.binSize_ms            = 200;         % FR sliding-bin width (ms).
    opt.stepSz_ms             = 20;          % FR sliding-bin step (ms).

    % B.8 POPULATION DYNAMICS (per-session, NGL02 path) 
    opt.popDyn.do             = false;         % master gate.
    % opt.popDyn.pca            = true;        % per-session PCA via calculate_pca_from_pool (single-trial overlay + CI tube).
    % opt.popDyn.jPCA           = false;       % rotational dynamics (PLACEHOLDER).
    % opt.popDyn.GPFA           = false;       % single-trial smooth trajectories (PLACEHOLDER).
    % opt.popDyn.trialEmbed     = false;       % legacy trial-similarity embedding.
    % opt.popDyn.smoothSigma    = 0.050;       % Gaussian smoothing sigma, SECONDS.
    % opt.popDyn.nComponents    = 3;           % output embedding dimensionality.
    % opt.popDyn.conditionVar   = '';          % per-condition grouping field name (empty = no grouping). Used by jPCA/GPFA/trialEmbed only.
    % opt.popDyn.dropAborted    = true;        % drop aborted trials before grouping.
    % opt.popDyn.alignIdx       = 1;           % 1-based index into opt.alignto. IGNORED by the new PCA path (it iterates all alignments). Kept for jPCA/GPFA placeholders.
    % opt.popDyn.trialEmbedMethod = 'tSNE';    % 'PCA' | 'tSNE' | 'UMAP'.
    % opt.popDyn.pcaConditions  = {'allInitiated'};  % cell of condition tokens for per-session PCA iteration. Each entry is one condition fieldname OR an 'X vs Y' comparison; produces one figure pair per (alignment, label, entry). E.g. {'allInitiated','correct vs incorrect'}.

    % B.9 LFP SIDE (NGL02_LFP - quick-look; research-grade analyses live in NGL07_LFPanalysis)
    opt.trialparsed           = false;       % load *_<align>.mat (trial-parsed) instead of *_FTcont.mat.
    opt.artifdet              = false;       % run LFP artifact detection / rejection.
    % opt.artZvalue             = 10;          % z-value cutoff for ft_artifact_zvalue. Used only if opt.artifdet=true.
    % opt.rejValue              = 'zero';      % how to fill rejected segments: 'zero' | 'nan' | numeric scalar.
    opt.spectrogram           = false;       % run multitaper TFR analysis (continuous or trial-parsed per opt.trialparsed).
    % TFR (trial-parsed) knobs consumed by trialparsed_MTspectrogram and plot_superletsTFR_*:
    opt.freqInterest          = {[1:1:29] [30:2:98] [100:5:250]};  % cell of freq vectors per band, e.g. {[4:1:30] [30:2:150]}. REQUIRED for trial-parsed TFR.
    % opt.TFRmethod             = 'wavelet';   % 'wavelet' | 'mtmconvol' |'superlet'. Automatically 'mtmconvol' for low frequencies
    % opt.superletOrder         = {};          % cell of superlet orders per band (only when TFRmethod = 'superlet').
    % opt.width                 = {};          % cell of superlet widths per band  (only when TFRmethod = 'superlet').
    % opt.combine               = 'additive';  % superlet combine mode.
    % opt.timeResol             = 0.2;         % TFR toi step (s). e.g. cfg.toi = opt.toi(1):timeResol:opt.toi(2).
    % opt.toi                   = [-0.5 4];      % TFR time-of-interest range [t0 t1] (seconds).
    % opt.blocks                = 'all';       % 'all' -> derive from condition.block, or integer N.
    % opt.chbych                = false;       % trigger plot_superletsTFR_extintion_chbych.
    % opt.trialbytrial          = false;       % trigger plot_superletsTFR_extintion_tbt.
    % opt.lfp.tfrAreaFilter     = '';          % restrict TFR to a subset of channels by chanArea; '' | 'NCL' | {'NCL','STR'}.
    % opt.lfp.tfrCacheDir       = '';          % override for the TFR cache folder; '' -> <analysis>/cache/lfp_tfr/.
    % opt.lfp.alignSubset       = {};          % restrict LFP path to a subset of opt.alignto (e.g. {'itiOn'}); {} -> all.
    % opt.lfp.downsampleFs      = [];          % explicit target Hz for ft_resampledata; overrides autoDownsample when set. [] -> per-band auto.
    % opt.lfp.autoDownsample    = true;        % per-band ft_resampledata to autoDownsampleFactor * max(band). ~5x speedup on low band.
    % opt.lfp.autoDownsampleFactor = 4;        % Nyquist safety (2x required; 4x = generous, 3x aggressive, 5x extra safe).
    % opt.lfp.autoMethod        = true;        % below autoMethodThresholdHz use mtmconvol+hanning (FT-recommended low-freq path; ~2x faster).
    % opt.lfp.autoMethodThresholdHz = 30;      % method-switch cutoff (Hz).
    % opt.lfp.autoFoi           = true;        % replace linear foi with quarter-octave log-spaced per band (halves nFreq).
    % opt.lfp.autoFoiStep       = 1/4;         % log-spacing step in octaves; 1/6 denser, 1/2 half-octave (coarser).
    % opt.lfp.parallel          = 'none';      % 'none' | 'bands' (parfor over freqInterest) | 'trials' (FT-native per-trial via cfg.parallel).
    
    % NGL07_LFPanalysis (research-grade session-level LFP; runs after NGL02_postPhy + NGL06_videoAnalysis)
    opt.lfp.session.do        = false;       % master gate for NGL07.
    % opt.lfp.session.tfr       = true;        % (a) trial-parsed TFR per opt.alignto.
    % opt.lfp.session.bursts    = false;       % (b) per-band burst detection (detectBursts).
    % opt.lfp.session.phase     = false;       % (c) continuous phase + envelope per band (hilbertBandpass).
    % opt.lfp.session.spikeField= false;       % (d) spike-field coupling PPC + coherence (needs NGL02_postPhy spike.mat).
    % opt.lfp.session.behReg    = false;       % (e) LFP x behaviour regression (STUB; needs NGL06 gaze.mat sidecar).
    % opt.lfp.session.flip      = false;       % (f) spectrolaminar (vFLIP) mapping.
    % opt.lfp.bands             = {{'theta',[4 8]},{'beta',[15 30]},{'gamma',[30 90]}}; % bands for (b) (c).
    % opt.lfp.burst.threshMult  = 3;           % detectBursts z-sigma threshold.
    % opt.lfp.burst.minDurationMs = 100;       % detectBursts minimum episode duration.
    % opt.lfp.spikeField.minSpikes = 50;       % spikeFieldCoupling: skip clusters below this.
    % opt.lfp.spikeField.ppcMethod = 'ppc2';   % 'ppc0'|'ppc1'|'ppc2' (Vinck et al.).
    % opt.lfp.spikeField.timwin = 0.5;         % STA window (s) around each spike.
    % opt.lfp.spikeField.foi    = 2:2:100;     % coherence frequency vector (Hz).
    % opt.lfp.hilbert.storePhase= true;        % save per-band phase matrix.
    % opt.lfp.hilbert.storeEnv  = true;        % save per-band envelope matrix.
    % opt.lfp.hilbert.dtype     = 'single';    % storage dtype (single|double).
    % opt.lfp.flip.laminaraxis  = 0:0.05:1.55; % probe laminar axis (mm).
    % opt.lfp.flip.freqaxis     = 1:150;       % vFLIP freq axis (Hz).
    % opt.lfp.flip.setfreqbool  = 0;           % 0 = vFLIP (auto), 1 = default fixed FLIP bands.
    % opt.lfp.plot.visible      = 'off';       % NGL02_LFP quick-look plots: figure visibility.
    % opt.lfp.plot.Resolution   = 300;         % exportgraphics DPI for the quick-look PNGs.
    % opt.lfp.plot.zlim         = [];          % colour-axis limits [zmin zmax]; [] -> auto.
    % opt.lfp.plot.colormap     = 'parula';    % colormap name or Nx3 matrix.
    % opt.lfp.plot.trialFilter  = 'correct';   % quick-look TFR trial mask: condition field to filter on ('' -> all trials).
    % opt.lfp.plot.baseline     = [-0.5 0];    % [t0 t1] seconds for dB baseline; [] -> raw power. Default matches new opt.toi = [-0.5 4].
    % opt.lfp.plot.interp       = 'bilinear';  % imagesc interpolation: 'bilinear' (smooth) | 'nearest' | 'none'.

    % B.10 NGL04_fireRate  (cross-subject PSTH plotter) 
    opt.fireRatePlot.interval        = [-1000 4000];      % ms window passed to plotPSTH.
    opt.fireRatePlot.binSize_ms      = opt.binSize_ms;    % FR sliding-bin width inside the plot.
    opt.fireRatePlot.stepSz_ms       = opt.stepSz_ms;     % FR sliding-bin step inside the plot.
    % opt.fireRatePlot.smoothPlot      = true;            % nanMeanSterrHistogram smoothing.
    % opt.fireRatePlot.errAlpha        = 0.4;             % error-shade alpha (0..1).
    % opt.fireRatePlot.labelPriority   = {'HumanLabel','KSLabel','bc_unitType'};   % resolution order for cluster-label tokens (also used by per-session calculate_neural_pca).
    % opt.fireRatePlot.busyWarnTraces  = 4;               % warn above N overlaid traces per subplot.
    % opt.fireRatePlot.outDir          = '';              % default: <input.analysis>/plots/NGL04_fireRate
    % opt.fireRatePlot.cacheDir        = '';              % default: <input.analysis>/cache/firepools (SHARED with NGL04_PCA; first run that pools a (align,cond,label) writes the cache, subsequent runs reuse it).

    % B.11 NGL04_PCA / per-session PCA  (state-space plots) 
    % NGL04_PCA reads these directly. calculate_neural_pca (per-session,
    % from NGL02) also reads .nBootstrap/.rngSeed/.sessionAlpha/.ciAlpha/
    % .ciStride/.variants/.outDir from here; for .interval it uses
    % opt.fireRatePlot.interval and for binning it uses opt.binSize_ms /
    % opt.stepSz_ms (see LANDMINE note at the top of this file).
    opt.pcaPlot.interval        = opt.fireRatePlot.interval;  % ms window around alignment (NGL04_PCA only).
    opt.pcaPlot.binSize_ms      = opt.binSize_ms;         % FR bin width (ms) inside NGL04_PCA only.
    opt.pcaPlot.stepSz_ms       = opt.stepSz_ms;          % FR bin step (ms) inside NGL04_PCA only.
    % opt.pcaPlot.smoothSigma_s  = 0.050;                 % Gaussian sigma for smoothing (s).
    % opt.pcaPlot.nComponents    = 3;                     % number of PCs to keep (>=2).
    % opt.pcaPlot.nBootstrap     = 100;                   % trial-bootstrap reps for CI tube; 0 disables CI.
    % opt.pcaPlot.rngSeed        = [];                    % integer seed for reproducible CI, or [] for random.
    % opt.pcaPlot.sessionAlpha   = 0.18;                  % alpha for grey single-trial / session-marginal traces.
    % opt.pcaPlot.ciAlpha        = 0.20;                  % alpha for 2D CI ribbon.
    % opt.pcaPlot.ciStride       = 10;                    % 3D CI crosshair every N bins.
    % opt.pcaPlot.variants       = {'singleTrials','ciTube'};  % which figure variants to render.
    % opt.pcaPlot.outDir         = '';                    % default for NGL04_PCA: <input.analysis>/plots/NGL04_PCA. Per-session PCA writes to <opt.analysis>/plots/population_dynamics/ regardless.

    % B.12 PROJECT-SPECIFIC GATES 
    % Off by default; turn on only for the matching paradigm. Code blocks
    % guarded by these flags live in the toolbox and stay dormant for
    % any other project.
    % opt.proj_chgDtctPCue       = false;      % Change-Detection P-Cue paradigm (trialdef /32).
    % opt.proj_socialLearning    = false;      % SocialLearning ASL (TFR testname, etc.).
    % opt.proj_extintion         = false;      % Extintion paradigm (fireRate_extintion, trialdef -1000 offset).
    % opt.proj_FLIP              = false;      % vFLIP laminar power analysis (was opt.FLIP).

    % B.12.a Video + social tracking
    % opt.offlineTrack          = false;       % run offline video blob detection (social-arena).
    % opt.useTrack              = false;       % spike-vs-social interaction indexing.

    % B.13 CROSS-SESSION AGGREGATION (NGL03_acrossSession) 
    % aggregateSubjects requires aggregateSessions=true; subjects can only
    % be stacked once sessions have been collapsed per subject.
    opt.aggregateSessions     = true;          % build <subject>_aggregated_<area>.mat per (subject, area).
    opt.aggregateSubjects     = true;          % build study-level aggregated_<area>.mat per area (requires aggregateSessions).

    % B.14 PER-AREA CONTEXT (set internally by NGL02) 
    % opt.area                  = 'all';       % single-area runs leave 'all'; multi-area mode iterates input.areaMap.uniqueAreas. Setting it here has no effect — NGL02 overwrites before each loadSpikes call.

    % B.15 GAZE/VIDEO PROCESSING (GazEstim) 
    % Wrapper around the GazEstim Python toolbox (toolboxes/GazEstim/{pose_clean,pose_render}.py
    % driven by configfiles/master_gaze.py). Pairs with B.12.a "Video +
    % social tracking" inside the future NGL06_VideoProcess stage.
    % Requires: copy configfiles/master_gaze.py (+ HexArena.png if used)
    % into your <project>\analysisCode\ and set GAZEpythonExe in
    % NGL_machineConfig.m.
    opt.gaze.do                 = false;       % master gate: run process_gaze on each discovered DLC csv.
    opt.gaze.features           = false;       % <output>_estimated_features.mat + a .png figure
    opt.gaze.video              = true;        % optional: features only, skip the mp4
    % opt.gaze.videoWidth       = 1250;        % source video width (px); must match DLC training.
    % opt.gaze.videoHeight      = 1160;        % source video height (px).
    % opt.gaze.parts            = {'beak','tail','wing_l','wing_r','back','HS_f'}; % can be a subset of the total
    % opt.gaze.roles.beak       = 'HS_f';      % which label will be the 'beak' if such does not exist
    % opt.gaze.roles.head       = 'HS_f';      %
    % opt.gaze.roles.left_wing  = 'wing_l';    % 
    % opt.gaze.roles.right_wing = 'wing_r';    %
    % opt.gaze.fps              = 59.94;       % input video fps (Hz).
    % opt.gaze.downsampleStep   = 2;           % keep 1 frame per N; out_fps = fps/N.
    % opt.gaze.pCut             = 0.5;         % DLC likelihood gate (0..1).
    % opt.gaze.devFac           = 0.6;         % jump threshold = devFac * bodyLength.
    % opt.gaze.smooth           = 5;           % temporal smoothing window (frames).
    % opt.gaze.wMed             = 9;           % rolling-median window (frames).
    % opt.gaze.boneTolFrac      = 0.4;         % bone-length tolerance fraction.
    % opt.gaze.boneTolMad       = 5.0;         % bone-length tolerance (xMAD).
    % opt.gaze.orderMargin      = 0.10;        % head-behind-wing clamp margin (xbody).
    % opt.gaze.drawCones        = true;        % render gaze cones on the output mp4.
    % opt.gaze.monoFOV          = 170;         % monocular field per eye (degrees).
    % opt.gaze.binoHalf         = 15;          % binocular half-angle (degrees).
    % opt.gaze.coneMult         = 2.5;         % cone length = coneMult * birdLength.
    % opt.gaze.eyeFwdFrac       = 1/3;         % eye base: fraction of head->beak from head.
    % opt.gaze.eyeLatFrac       = 1/5;         % eye lateral offset: fraction of back->wing.
    % opt.gaze.dpi              = 120;         % render DPI.
    % opt.gaze.crf              = 24;          % ffmpeg CRF (lower = better quality, larger files).
    % opt.gaze.preset           = 'veryfast';  % ffmpeg preset.
    % opt.gaze.background       = '';          % '' -> <analysisCode>/HexArena.png (NGL06 asserts it exists).
    % opt.gaze.masterScript     = '';          % '' -> <analysisCode>/master_gaze.py (NGL06 asserts it exists).
    % opt.gaze.pythonExe        = '';          % '' -> input.GAZEpythonExe (falls back to 'python' via NGL_machineConfig).
    % opt.gaze.previewFrame     = [];          % render one still PNG at this frame idx instead of the full mp4.
    % opt.gaze.targetFps        = [];          % alternative to downsampleStep; pick step so out_fps ~= targetFps.
    % opt.gaze.startTime        = [];          % [s] analyse from this time; [] -> from first frame.
    % opt.gaze.endTime          = [];          % [s] analyse up to this time; [] -> to last frame.
    % opt.gaze.maxFrames        = [];          % cap rendered output frames; [] -> unbounded.
    % opt.gaze.maxSeconds       = [];          % cap rendered output duration [s]; [] -> unbounded.
    % opt.gaze.overwrite        = false;       % true -> re-render even if output file already exists.

% C) PARAM (analysis/plot tuning, kept as inline-defaulted in functions).
%    `param` is a semi-independent struct that downstream analysis/plotting
%    functions fill with their own inline defaults if absent. Override here
%    only when you want non-default behaviour for THIS project.
param = struct();

    % Firing-rate analysis (calculate_fireRate_general)
    param.trial2plot      = 'allInitiated'; % 'allInitiated'|'correct'|'incorrect'|'omission'|...
    param.binSize         = opt.binSize_ms;           % ms; mirrors opt.binSize_ms.
    param.stepSz          = opt.stepSz_ms;            % ms; mirrors opt.stepSz_ms.
    param.interval        =  opt.fireRatePlot.interval; % ms window around alignment.
    param.baseline        = -param.interval(1);          % ms before alignment used for normalisation.
    param.smpRate         = 1000;          % samples/s used inside calcFireRate.
    param.blockchange     = [];            % trial indices marking block boundaries.

    % Plotting (plot_single_fireRate, plot_multi_fireRate)
    param.plot            = true;          % draw per-cluster + multi-cluster FR plots.   
    param.post            = param.interval(2);          % ms after event in raster plots.
    param.treatment       = false;         % colour by treatment level.
    param.plotevent       = 3;             % mark this many trial events in rasters.
    param.visible         = 'off';         % 'on'|'off' figure visibility.
    param.size            = 'adaptive';    % figure size policy.
    param.Resolution      = 300;           % exportgraphics resolution.

    % Block-aware FR (calculate_fireRate_byBlock)
    % param.blockBounds     = [];            % vector [t0 t1 ... tN] of trial indices defining N blocks.

%% SET 2 MERGE INTAN SESSIONS
% subject     = 'NNN';                   % subject ID (string or char)
% mergeDates  = {'' ''}; % {A, B} in chronological order
% mergeTag    = 'AABB_merged';          % user-supplied output folder name, preferably NOT a date (change it later)
% opt.merge.cleanupStaging = true;

% and RUN this ONLY
% NGL_mergeSessionsINTAN

%% 3) RUN.
%% 1  NGL01_Main — preprocessing.
%   Locates sessions, determines file format, extracts events and motion
%   data, converts to Kilosort and FieldTrip formats, runs Kilosort 4
%   spike sorting, optionally runs Bombcell QC, and (optionally) opens
%   Phy for manual curation.
%
%   PRODUCES (per session):
%     <session>.bin, kilosort\, EventRecord.mat, trialdef.mat, events.mat,
%     *_FTcont.mat / *_<event>.mat, preprocInfo.mat (per-session), and
%     analysisCode\preprocInfo_lastRun.mat (study-level snapshot).

NGL01_Main

%% 2  NGL02_postPhy — spike pipeline (REQUIRES Phy curation).
%   Loads curated KS/Phy clusters, sorts into trials, computes firing
%   rate, optionally runs population-dynamics analysis. Multi-area
%   aware: nested spike.<Area> / neurons.<Area> / fireRate.<Area> when
%   input.Areas is set.
%
%   PRODUCES (per session, when enabled):
%     spike.mat, neurons.mat, fireRate.mat, optional neuralDynamics.mat

NGL02_postPhy

%% 3  NGL02_LFP — LFP pipeline.
%   Independent of Phy curation: can be run any time after NGL01_Main
%   finishes (e.g., in parallel with manual curation). Loads
%   FieldTrip-formatted LFP, runs artifact rejection and time-frequency
%   analysis.
%
%   PRODUCES (per session, when enabled):
%     *_FT_data_NoArtif.mat (if opt.artifdet),
%     *_TFR_continuous.mat / *_TFR_trialparsed.mat (if opt.spectrogram)
%
%   If LFP work is not needed for this project, set opt.doLFPthing = false
%   or skip this section.

NGL02_LFP

%% 4  NGL03_acrossSession — cross-session and cross-subject aggregation.
%   Aggregates per-session outputs (spike / neurons / fireRate /
%   condition / events / trialdef; plus neuralDynamics and blob when
%   their gates are on) into cell arrays indexed by (subject, session).
%   Gated by opt.aggregateSessions (per-subject) and opt.aggregateSubjects
%   (study-level). Both default to false; set them in section 2B above.
%
%   PRODUCES (one file per area; see migrate_aggregated_to_perArea.m for
%   upgrading legacy nested files):
%     data\analysis\<subject>\<subject>_aggregated_<area>.mat   (per (subject, area))
%     data\analysis\aggregated_<area>.mat                       (per area, across subjects)
%   Single-area studies use area='main'.
%
%   LFP-side aggregation is planned but not yet wired here; the spike
%   side is the foundation.

NGL03_aggregate

%% 5  NGL04_fireRate — cross-subject FR PSTH plots.
%   Consumes the NGL03_acrossSession output. Set `request` as a 1x3
%   cell, each entry a single value or a 'X vs Y' comparison. Categories
%   (condition, cluster label, alignment) are inferred from content;
%   slot order is irrelevant. Zero to three entries may be 'vs'.
%
%   Layout rule:
%     ALIGNMENT controls subplot layout (side-by-side, one per align).
%     CONDITION x LABEL overlay within each subplot.
%
%   Examples:
%     request = {'correct vs incorrect', 'good', 'stim2'};
%       % 1 subplot, 2 overlaid traces (cond varies)
%     request = {'correct', 'good', 'itiOn vs stim2'};
%       % 2 subplots side-by-side, 1 trace each (align varies)
%     request = {'correct vs incorrect', 'good', 'itiOn vs stim2'};
%       % 2 subplots side-by-side, 2 overlaid traces each
%

request = {'correct', 'good', 'stim2'};
NGL04_fireRate

%% 6  NGL04_PCA — cross-subject population PCA state-space plots.
%   Same 3-cell `request` semantics as NGL04_fireRate. Reuses the shared
%   firepools cache (opt.fireRatePlot.cacheDir), so if NGL04_fireRate was
%   run first for this request, the pools are loaded from disk and only
%   the PCA fit + plot run here. Produces two PNG variants per
%   (alignment, label):
%     <encoded-request>_a<A>_l<L>_singleTrials.png   (mean + per-session
%                                                     grey marginals)
%     <encoded-request>_a<A>_l<L>_ciTube.png         (mean + bootstrap CI)
%   plus an <encoded-request>_pca.mat with all PCA results.
%
%   Examples:
%     request = {'correct vs incorrect', 'good', 'stim2'};
%       % one PC space, 2 trajectories (correct vs incorrect overlaid)
%     request = {'allInitiated', 'good', 'itiOn vs stim2'};
%       % two figures (one PC space per alignment), 1 trajectory each

request = {'correct vs incorrect', 'good', 'stim2'};
NGL04_PCA

%% 7  NGL03_plotting — group-level visualisation (TODO).
%   Comprehensive plots across sessions and conditions. Building on the
%   per-session plots already produced by NGL02_postPhy.

% NGL03_plotting

%% 8  NGL06_videoAnalysis — video-based gaze pipeline.
%   Fans out the GazEstim Python pipeline (pose_clean + pose_render,
%   driven by configfiles/master_gaze.py) across every discovered
%   (subject, session). One DLC csv per session is expected at
%   <input.bhvfolder>/<subject>/<session>/*.csv; NGL06 writes the
%   gaze-cone-overlaid mp4 (or single PNG when opt.gaze.previewFrame is
%   set) next to the csv.
%
%   PROJECT-LOCAL COPIES REQUIRED at <input.analysisCode>/:
%     master_gaze.py   copy once from <toolbox>/configfiles/master_gaze.py
%     HexArena.png     copy once from <toolbox>/configfiles/HexArena.png
%
%   Toggle opt.gaze.do = true in B.16 to enable. See docs/gaze_pipeline.md.

% NGL06_videoAnalysis

%% 9  NGL07_LFPanalysis - research-grade session-level LFP analyses.
%   Consumes NGL01 FTcont + NGL02_postPhy spikes + (optional) NGL06 gaze
%   sidecar, and runs the analyses enabled under opt.lfp.session.*:
%     (a) trial-parsed TFR  (b) burst detection  (c) phase + envelope
%     (d) spike-field PPC / coherence  (e) LFP x behaviour regression
%
%   Multi-area handled via FT_data.chanArea + opt.lfp.tfrAreaFilter -
%   the SAME FT file feeds each per-area subselect (Q3 (b) design).
%
%   See docs/lfp_pipeline.md for the full opt.lfp.* reference.

% NGL07_LFPanalysis

%% More custom stages...
% NGLXX_something
