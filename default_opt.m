function opts = default_opt()
% default_opt  Return the canonical default options struct for the NGL toolbox.
%
% PURPOSE:
%   Single source of truth for every configurable option. Every recognised
%   option name MUST appear here with a safe default value. Downstream
%   functions must never hard-code their own defaults — they rely on this
%   struct already being complete and validated before they are called.
%
% USAGE:
%   opts = default_opt();
%   Called internally by set_default. Users should not call this directly;
%   instead, set options in the opt struct inside NGL_SetAndRunMe.m.
%
% OUTPUT:
%   opts - struct with all NGL option fields pre-filled to safe defaults.
%          Field names are the canonical names recognised by set_default.
%
% ADDING A NEW OPTION:
%   1. Add the field here with its default value and a comment.
%   2. Add validation logic in set_default (Section 2) if needed.
%   3. Document it in wiki_NGL01_pipeline.md (Section 6).
%
% Last modified 07.05.2026 (Jesus)

    % Data format
    opts.numChannels     = 32;      % Expected channel count (override for 32-ch Deuteron)
    opts.bin             = true;    % Normally, we always check if the .bin file exists
    opts.FieldTrip       = true;    % Produce a FieldTrip-ready .mat file
    opts.doNWB           = true;    % NWB export: INTAN via NeuroConv (Python), Deuteron via matNWB (MATLAB)

    % Events
    opts.RetrieveEvents  = true;     % Extract event log from session
    opts.alignto         = {'itiOn'};% Alignment events; cell array of char vectors
    opts.trEvents        = {};       % 'Special' ITI events (treatments, tutors, etc.)
    opts.addtime         = 0;        % Padding around trial start/end in ms
    opts.uselog          = false;    % By default, use Deuteron data files to extract events. 
                                      % When true, uses the text log. For cases when the events 
                                      % were not properly transmitted to the system but logged.

    % Motion sensors
    opts.GetMotionSensors = false;   % Extract head-direction sensor data

    % Data preprocessing
    opts.noise           = [];       % Reserved for noise-rejection parameters
    opts.lowpass         = 9000;     % High boundary for low-pass (Hz). [] = off.
    opts.lowpassFT       = 200;      % Low-pass for FieldTrip LFP stream (Hz)
    opts.highpass        = [];       % Low boundary for high-pass (Hz). [] = off.
    opts.linefilter      = 0;        % Line-noise notch centre frequency. 0 = off.
    opts.CAR             = 0;        % Common-average re-referencing. 0 = off.
    opts.dwnsmplRate     = [];       % LFP downsample target (Hz). [] = auto (937.5 Hz).
    opts.timebreak       = false;    % If a break in the recording is expected (e.g. Deuteron battery change)

    % Sorting & curation
    opts.kilosort        = true;     % Default to Kilosort4
    opts.KSchanMapFile   = '';       % Empty = linear array; set to 'chanMapXXX.mat' for custom
    opts.bombcell        = true;     % Run Bombcell QC on Kilosort output
    opts.callBcGUI       = false;    % GUI after BC metrics
    opts.phy             = false;    % Open Phy after sorting (blocks MATLAB)

    % NGL02_postPhy options
    opts.doSpikething    = true;     % Process single-unit/spike data
    opts.doLFPthing      = true;     % Process LFP data
    opts.offlineTrack    = false;    % Run offline video blob detection
    opts.useTrack        = false;    % Index spiking against social-tracking events
    opts.trialparsed     = false;    % Load trial-parsed FT file (vs continuous)
    opts.artifdet        = false;    % Run LFP artifact detection and rejection
        opts.artZvalue   = 10;       % z-value cutoff for ft_artifact_zvalue
        opts.rejValue    = 'zero';   % value to insert into rejected segments ('zero'|'nan'|numeric)
    opts.spectrogram     = false;    % Run multitaper time-frequency analysis
    opts.neurDyn.do      = false;    % LEGACY trial-state embedding (kept for back-compat;
                                     % retired in favour of opt.popDyn below).

    % Population-dynamics family (NGL02 -> calculate_population_dynamics wrapper).
    % Each .<method> flag opts that method in/out independently; the wrapper
    % aggregates results into a single neuralDynamics struct.
    opts.popDyn          = struct( ...
        'do',           false,  ...  % master gate: run any population-dynamics step
        'pca',          true,   ...  % trial-averaged smoothed-rate PCA (real time-trajectories)
        'jPCA',         false,  ...  % rotational dynamics (PLACEHOLDER: not yet implemented)
        'GPFA',         false,  ...  % single-trial smooth trajectories (PLACEHOLDER)
        'trialEmbed',   false,  ...  % legacy trial-similarity embedding (PCA/tSNE/UMAP per trial)
        'smoothSigma',  0.050,  ...  % Gaussian smoothing kernel sigma, SECONDS (~50 ms)
        'nComponents',  3,      ...  % output embedding dimensionality
        'conditionVar', '',     ...  % field name on `condition` for per-condition grouping (empty = no grouping)
        'dropAborted',  true,   ...  % drop aborted trials before grouping (matches legacy fireRate filter)
        'alignIdx',     1,      ...  % which opt.alignto entry to analyse (popDyn methods operate on one alignment at a time)
        'trialEmbedMethod','tSNE');  % method used by the legacy trialEmbed view: 'PCA'|'tSNE'|'UMAP'

    % Waveform extraction (loadSpikes; consumed during NGL02)
    opts.getwF           = false;    % Extract raw waveforms per cluster (slow)
    opts.gwfparams       = struct( ...
        'wfWin',    [-20 41], ...    % samples around spiketime (negative=before)
        'nWf',      2000,     ...    % max waveforms per cluster
        'dataType', 'int16',  ...    % .bin sample type (overridden by loadSpikes)
        'nCh',      []);             % set inside loadSpikes from spikes.n_channels_dat
    % Notes on gwfparams:
    %   - wfWin and nWf are the fields users typically tune.
    %   - dataType and nCh are kept for completeness but loadSpikes derives
    %     them at runtime (hard-codes int16 and reads channel count from
    %     the loaded KS struct). Setting them in NGL_SetAndRunMe has no
    %     effect today; left here for compatibility.

    % Cluster loading and ISI binning (loadSpikes / calc_isihist).
    opts.spparams        = struct( ...
        'excludeNoise', true, ...    % skip Phy 'noise'-labelled clusters
        'loadPCs',      false);      % load principal components (rarely needed)
    opts.isibins         = 0:0.5:200;% ISI histogram bin edges, ms

    % Firing-rate binning (canonical values, shared by calcFireRate's
    % `param` defaults and consumed by the population-dynamics family
    % to size its smoothing kernel). Keep these in sync with whatever
    % the user overrides in `param.binSize` / `param.stepSz` at call
    % time — if they diverge, set_default does not (and cannot) check
    % per-call param values, so the popDyn smoothing will use the
    % canonical opt values regardless.
    opts.binSize_ms      = 200;      % FR sliding-bin width (ms)
    opts.stepSz_ms       = 20;       % FR sliding-bin step (ms)

    % Per-area context, set by NGL02 around each loadSpikes call.
    % Single-area runs leave this as 'all'; multi-area runs set it to the
    % current area name (one of input.areaMap.uniqueAreas) before calling
    % loadSpikes, which uses it to tag every cluster's spike.roi.
    opts.area            = 'all';

    % Project-specific gates (opt-in code paths for specific paradigms).
    % Each flag guards code blocks that would otherwise be commented out
    % or hard-wired for one paradigm. Default false everywhere; users
    % opt in from their NGL_SetAndRunMe. Renaming convention: opt.proj_*.
    opts.proj_chgDtctPCue    = false;  % Change-Detection P-Cue paradigm fixes (trialdef /32).
    opts.proj_socialLearning = false;  % SocialLearning ASL: TFR testname, etc.
    opts.proj_extintion      = false;  % Extintion paradigm: fireRate_extintion, trialdef -1000 offset.
    opts.proj_FLIP           = false;  % vFLIP laminar power analysis (renamed from opts.FLIP).

    % Cross-subject FR PSTH plotter (NGL04_fireRate). All defaults are
    % overridable from NGL_SetAndRunMe via opt.fireRatePlot.<field>.
    opts.fireRatePlot = struct( ...
        'interval',       [-2000 4000],                              ...  % ms window around alignment
        'binSize_ms',     200,                                        ...  % overrides inside the plot
        'stepSz_ms',      20,                                         ...  % overrides inside the plot
        'smoothPlot',     true,                                       ...  % nanMeanSterrHistogram smoothing
        'errAlpha',       0.4,                                        ...  % error-shade alpha (0..1)
        'labelPriority',  {{'HumanLabel','KSLabel','bc_unitType'}},   ...  % resolution order for cluster-label tokens
        'busyWarnTraces', 4,                                          ...  % warn above N overlaid traces per subplot
        'outDir',         '');                                              % '' -> default <input.analysis>/plots/NGL04_fireRate

    % Cross-session aggregation (NGL03_acrossSession).
    % aggregateSubjects requires aggregateSessions: subjects can only be
    % stacked after each subject's sessions have been collapsed into a
    % per-subject .mat. Both default to false so the aggregation stage
    % is opt-in.
    opts.aggregateSessions = false;  % build <subject>_aggregated.mat per subject
    opts.aggregateSubjects = false;  % build study-level aggregated.mat across subjects

end
