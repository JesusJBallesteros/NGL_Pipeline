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
%   1) PREPARE   - set study metadata, drive and project name, run NGL00_Prep
%   2) SET       - choose subjects/sessions and configure opt / param
%                    2A) subjects + dates
%                    2A2) multi-area (optional)
%                    2B) opt struct (required + commonly-defaulted)
%                    2C) param struct (analysis/plot tuning, kept semi-independent)
%   3) RUN       - call the pipeline stages in order:
%                    3.1  NGL01_Main           preprocessing
%                    3.2  NGL02_postPhy        spike (requires Phy curation)
%                    3.3  NGL02_LFP            LFP (independent of Phy)
%                    3.4  NGL03_plotting       (TODO)
%                    3.5  NGL03_acrossSession  cross-session aggregation
%
% REQUIRED CONFIG FILES (place in analysisCode\ alongside this script):
%   NGL_machineConfig.m   - machine-specific Python env and toolbox paths
%   eventDefinitions.m    - project event code definitions
%   conditions_script.m   - condition grouping logic
%   chanMapXXX.mat        - custom Kilosort channel map (if applicable)
%
% NOTES:
%   - Fields shown UNCOMMENTED are the ones most users override; fields
%     COMMENTED-OUT show their default value for reference — uncomment
%     only to override the default.
%   - The full canonical defaults live in default_opt.m and are validated
%     by set_default.m.
%   - NGL_machineConfig.m is machine-specific; do not commit it to git.
%
% Last modified 29.05.2026 (Jesus) - comprehensive opt/param reference

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

% input.Areas = {'NCL','NCL','STR'};   % e.g. 2 NCL shanks + 1 STR shank.

% B) OPTIONS (opt struct).
%    Fields not set here fall back to defaults in default_opt.m;
%    set_default.m validates every field. Commented entries below show
%    the default, uncomment ONLY to override.
opt = struct();

    % Acquisition 
    opt.numChannels         = 32;            % electrode channel count.
    % opt.bin               = true;          % create Kilosort .bin file.
    % opt.FieldTrip         = true;          % produce FieldTrip LFP .mat.
    % opt.doNWB             = true;          % NWB export (NeuroConv/matNWB).

    % Events 
    opt.RetrieveEvents      = true;          % extract event log.
    opt.alignto             = {'itiOn'};     % alignment events (cell of char).
    % opt.trEvents          = {};            % inter-trial events (block changes, treatments).
    % opt.addtime           = 0;             % ms padding around trial start/end.
    % opt.uselog            = false;         % Deuteron text-log fallback.

    % Motion sensors 
    % opt.GetMotionSensors  = false;         % head-direction sensor data.

    % Filtering / preprocessing
    % opt.lowpass           = 9000;          % bin file low-pass, Hz. [] = off.
    % opt.lowpassFT         = 200;           % FieldTrip LFP low-pass, Hz.
    % opt.highpass          = [];            % bin file high-pass, Hz. [] = off.
    % opt.linefilter        = 0;             % line-noise notch centre, Hz. 0 = off.
    % opt.CAR               = 0;             % common-average rereferencing. 0 = off.
    % opt.dwnsmplRate       = [];            % LFP downsample target, Hz. [] = auto, 1 KHz
    % opt.timebreak         = false;         % expect a recording break (Deuteron battery change).
    % opt.noise             = [];            % reserved for noise-rejection params.

    % Sorting & curation 
    % opt.kilosort          = true;          % run Kilosort 4.
    opt.KSchanMapFile       = '';            % '' = linear array; or e.g. 'chanMap_ATLAS_E32-...mat'.
    % opt.bombcell          = true;          % run Bombcell QC after sorting.
    % opt.callBcGUI         = false;         % open Bombcell GUI.
    % opt.phy               = false;         % open Phy right after each session sorting (blocks MATLAB).

    % NGL02 stage gates 
    % opt.doSpikething      = true;          % run NGL02_postPhy spike work for each session.
    % opt.doLFPthing        = true;          % run NGL02_LFP work for each session.

    % Spike side: video + social tracking
    % opt.offlineTrack      = false;         % run offline video blob detection (social-arena projects).
    % opt.useTrack          = false;         % spike-vs-social interaction indexing.

    % Spike side: waveform extraction (loadSpikes)
    % opt.getwF             = false;         % extract raw waveforms per cluster (slow).
    % opt.gwfparams.wfWin   = [-20 41];      % samples around spiketime (negative = before).
    % opt.gwfparams.nWf     = 2000;          % max waveforms per cluster.
    % opt.gwfparams.dataType = 'int16';      % overridden internally by loadSpikes.
    % opt.gwfparams.nCh     = [];            % overridden internally by loadSpikes.

    % Spike side: cluster loading + ISI binning
    % opt.spparams.excludeNoise = true;      % skip Phy 'noise' clusters.
    % opt.spparams.loadPCs      = false;     % load principal components.
    % opt.isibins           = 0:0.5:200;     % ISI histogram bin edges, ms.

    % Firing-rate binning
    % opt.binSize_ms        = 200;           % FR sliding-bin width (ms).
    % opt.stepSz_ms         = 20;            % FR sliding-bin step (ms).

    % Population dynamics
    % opt.popDyn.do          = false;        % master gate.
    % opt.popDyn.pca         = true;         % trial-averaged smoothed-rate PCA (real time-trajectories).
    % opt.popDyn.jPCA        = false;        % rotational dynamics.
    % opt.popDyn.GPFA        = false;        % single-trial smooth trajectories.
    % opt.popDyn.trialEmbed  = false;        % legacy trial-similarity embedding.
    % opt.popDyn.smoothSigma = 0.050;        % Gaussian smoothing sigma, SECONDS.
    % opt.popDyn.nComponents = 3;            % output embedding dimensionality.
    % opt.popDyn.conditionVar = '';          % field name on `condition` for per-condition grouping (empty = no grouping).
    % opt.popDyn.trialEmbedMethod = 'tSNE';  % 'PCA' | 'tSNE' | 'UMAP'.

    % LFP side (NGL02_LFP)
    % opt.trialparsed       = false;         % load *_stimOn2.mat (trial-parsed) instead of *_FTcont.mat.
    % opt.artifdet          = false;         % run LFP artifact detection / rejection.
    % opt.spectrogram       = false;         % run multitaper TFR analysis.
    % opt.FLIP              = false;         % vFLIP laminar power analysis (PLACEHOLDER, project-specific).

    % Legacy (kept for backward compatibility)
    % opt.neurDyn.do        = false;         % LEGACY trial-state embedding; superseded by opt.popDyn.

% C) PARAM (analysis/plot tuning, kept as inline-defaulted in functions).
%    `param` is a semi-independent struct that downstream analysis/plotting
%    functions fill with their own inline defaults if absent. Override here
%    only when you want non-default behaviour for THIS project.
param = struct();

    % Firing-rate analysis (calculate_fireRate_general)
    % param.trial2plot      = 'allInitiated'; % 'allInitiated'|'correct'|'incorrect'|'omission'|...
    % param.binSize         = 200;           % ms; mirrors opt.binSize_ms.
    % param.stepSz          = 20;            % ms; mirrors opt.stepSz_ms.
    % param.interval        = [-2000 10000]; % ms window around alignment.
    % param.smpRate         = 1000;          % samples/s used inside calcFireRate.
    % param.baseline        = 2000;          % ms before alignment used for normalisation.
    % param.plot            = true;          % draw per-cluster + multi-cluster FR plots.
    % param.blockchange     = [];            % trial indices marking block boundaries.

    % Plotting (plot_single_fireRate, plot_multi_fireRate)
    % param.visible         = 'off';         % 'on'|'off' figure visibility.
    % param.Resolution      = 300;           % exportgraphics resolution.
    % param.treatment       = false;         % colour by treatment level.
    % param.post            = 2000;          % ms after event in raster plots.
    % param.plotevent       = 3;             % mark this many trial events in rasters.
    % param.size            = 'adaptive';    % figure size policy.

    % Block-aware FR (calculate_fireRate_byBlock)
    % param.blockBounds     = [];            % vector [t0 t1 ... tN] of trial indices defining N blocks.


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

%% 4  NGL03_acrossSession — cross subject and cross-session aggregation (TODO)
%   Aggregate per-session outputs (spike / neurons / fireRate / condition
%   / TFR) into study-level structures.

% NGL03_acrossSession

%% 4  NGL03_plotting — group-level visualisation (TODO).
%   Comprehensive plots across sessions and conditions. Building on the
%   per-session plots already produced by NGL02_postPhy.

% NGL04_plotting

%% More custom stages...
% NGLXX_something
