# NGL Pipeline — Operator Guide
**Running `NGL01_Main`: from raw recordings to sorted spikes and LFP**

*Last updated: May 2026 · Maintainer: Jesus Ballesteros · Repository: `ephys-data-pipeline`*

---

## Table of Contents

1. [Overview](#1-overview)
2. [Prerequisites](#2-prerequisites)
3. [Folder Structure (IKN Standard)](#3-folder-structure-ikn-standard)
4. [Machine Configuration (`NGL_machineConfig`)](#4-machine-configuration-ngl_machineconfig)
5. [User Configuration (`NGL_SetAndRunMe`)](#5-user-configuration-ngl_setandrunme)
6. [All Options Reference (`opt` fields)](#6-all-options-reference-opt-fields)
7. [Event Coding System](#7-event-coding-system)
8. [Customising `eventDefinitions.m`](#8-customising-eventdefinitionsm)
9. [Pipeline Stages — What Happens When You Run](#9-pipeline-stages--what-happens-when-you-run)
10. [Outputs](#10-outputs)
11. [Running Post-Phy Processing (`NGL02_postPhy`)](#11-running-post-phy-processing-ngl02_postphy)
12. [Troubleshooting](#12-troubleshooting)

---

## 1. Overview

The NGL pipeline converts raw electrophysiology recordings — acquired with **INTAN RHX** or **Deuteron** hardware — into analysis-ready formats:

- **`.bin` file** for spike sorting with [Kilosort 4](https://github.com/MouseLand/Kilosort)
- **FieldTrip `.mat` structure** for LFP analysis
- **`EventRecord.mat` / `trialdef.mat` / `events.mat`** for trial-parsed alignment
- **Bombcell QC report** for semi-automatic cluster curation

The pipeline is driven by a single user-facing script, **`NGL_SetAndRunMe.m`**, which lives in your project's `analysisCode` folder. You only ever need to edit that one file.

```
NGL_SetAndRunMe  →  NGL00_Prep  →  NGL01_Main  →  [Phy]  →  NGL02_postPhy
```

---

## 2. Prerequisites

### MATLAB toolboxes (bundled under `toolboxes/`)

| Toolbox | Purpose |
|---|---|
| `fieldtrip_light` | LFP preprocessing and FieldTrip format |
| `Intan` | INTAN RHD/RHS header reader |
| `Deuteron` | Deuteron binary reader + event DLL |
| `bombcell` | Automatic spike-cluster QC |
| `npy-matlab` | NumPy `.npy` reader (KS4 output) |
| `prettify_matlab` | Plot formatting |
| `spikes` | Spike analysis utilities |
| `BDPAT_NGL` | NGL-specific analysis helpers |
| `Viewer` | Data viewer |

All paths are added automatically by `set_default` on each run — **do not add them manually** to MATLAB's persistent path.

### Python environments

Two separate Conda environments are required:

| Environment | Purpose | Config key |
|---|---|---|
| KS4 Python env | Kilosort 4 spike sorting | `cfg.KSpythonExe` |
| PHY Python env | Phy manual curation | `cfg.PHYpythonExe` |

Kilosort 4 must be installed (`pip install kilosort`) inside its environment. NeuroConv is only needed if `opt.doNWB = true`.

### Project files (in `analysisCode/`)

| File | Required? | Notes |
|---|---|---|
| `NGL_SetAndRunMe.m` | Yes | Copy template from toolbox root, edit for your project |
| `NGL_machineConfig.m` | Yes | Machine-specific paths |
| `eventDefinitions.m` | Yes | Copy template from `configfiles/`, add project events |
| `conditions_script.m` | Yes | Project-specific condition grouping |
| `chanMapXXX.mat` | If custom probe | Custom channel map for Kilosort |

---

## 3. Folder Structure (IKN Standard)

`NGL00_Prep` creates the full folder tree on first run. The expected layout is:

```
<datadrive>:\<studyname>\
│
├── analysisCode\               ← NGL_SetAndRunMe, eventDefinitions, chanMaps, etc.
│
├── data\
│   ├── raw\
│   │   └── <subject>\
│   │       └── <YYYYMMDD>\    ← raw INTAN .dat / Deuteron .DT2 / .DF1 files
│   │
│   ├── preprocessing\
│   │   └── <subject>\
│   │       └── <YYYYMMDD>\    ← .bin, kilosort\4\, EventRecord.mat
│   │
│   ├── spikeSorted\
│   │   └── <subject>\
│   │       └── <YYYYMMDD>\    ← post-Phy spike data
│   │
│   ├── trialSorted\
│   │   └── <subject>\
│   │       └── <YYYYMMDD>\    ← trialdef.mat, events.mat, FT .mat files
│   │
│   ├── behaviour\
│   │   └── <subject>\
│   │       └── <YYYYMMDD>\    ← behavioural files
│   │
│   └── analysis\
│       └── <subject>\
│           └── <YYYYMMDD>\    ← analysis outputs
│
├── manuscript\
├── paradigmCode\
└── training\
```

The `readme.txt` in the study root is written once from `readmecontent` in `NGL_SetAndRunMe`.

> **IKN standard reference:** `gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure`

---

## 4. Machine Configuration (`NGL_machineConfig`)

`NGL_machineConfig.m` must exist in `analysisCode/` and return a struct with machine-specific absolute paths. Copy the template from `configfiles/` and fill in your paths:

```matlab
function cfg = NGL_machineConfig()
    cfg.toolbox      = 'C:\Code\ephys-data-pipeline';   % toolbox root
    cfg.KSpythonExe  = 'C:\envs\kilosort4\';            % KS4 conda env (contains python.exe)
    cfg.PHYpythonExe = 'C:\envs\phy2\';                 % Phy conda env
    % cfg.NCpythonExe  = 'C:\envs\neuroconv\';           % NeuroConv env (only if doNWB=true)
end
```

> `NGL_machineConfig` is **per machine, not per project**. Do not commit machine-specific paths to the shared repository. Add `NGL_machineConfig.m` to `.gitignore` or keep a generic template.

---

## 5. User Configuration (`NGL_SetAndRunMe`)

`NGL_SetAndRunMe.m` is the **only file you need to edit** for each project or recording batch. It lives in `analysisCode/` alongside your data-specific config files.

### Step-by-step walkthrough

**Section 1 — Prepare**

```matlab
clear all

readmecontent = ["Study name: MyStudy", ...
                 "Readme date: 01/01/2026", ...
                 "Person responsible: Your Name", ...
                 "Hardware used: INTAN RHX + Cambridge NeuroTech P32", ...
                 "Short description: ..."];

datadrive = 'D';           % drive letter where the study folder lives
studyname = 'MyStudy';     % must match the folder name on disk

NGL00_Prep                 % creates folder structure if needed
```

**Section 2 — Set**

```matlab
subjects = {'ABC'};                          % or 'all' to process every subject
dates    = {'20260101', '20260102'};         % or 'all' for every session

opt = struct();
opt.numChannels   = 32;                      % electrode count
opt.KSchanMapFile = 'chanMap_P32.mat';       % '' for linear array
opt.alignto       = {'itiOn', 'stimOn1'};    % events to align LFP trials to
opt.FieldTrip     = true;                    % produce LFP .mat file
opt.bombcell      = true;                    % run QC after sorting

NGL01_Main
```

**Section 3 — Run post-Phy**

After manual Phy curation:

```matlab
cd(input.analysisCode)
postPhy_param();
cd(toolbox)
NGL02_postPhy
```

---

## 6. All Options Reference (`opt` fields)

Every field must be set in `default_opt.m`. If you set a field in `NGL_SetAndRunMe` that does not appear in `default_opt`, `set_default` will issue a warning and ignore it.

### Data format

| Field | Default | Description |
|---|---|---|
| `numChannels` | `32` | Expected electrode channel count. Override to `64` for SpikeLog-64C Deuteron. For INTAN, the actual count is read from the header — set this to match. |
| `bin` | `true` | Create the `.bin` file for Kilosort. Set `false` to skip if the `.bin` already exists. |
| `FieldTrip` | `false` | Produce a FieldTrip-ready `.mat` file from the LFP stream. |
| `doNWB` | `false` | Run INTAN→NWB conversion via NeuroConv (testing; requires Python NeuroConv env). |

### Events

| Field | Default | Description |
|---|---|---|
| `RetrieveEvents` | `true` | Extract event log from the session. |
| `alignto` | `{'itiOn'}` | Event(s) to use as trial time-zero. Single char or cell array of chars matching names in `eventDefinitions`. At minimum keep `'itiOn'`. |
| `trEvents` | `{}` | "Special" ITI events (e.g. drug treatments, tutor calls). Cell array of event names. |
| `addtime` | `0` | Extra padding in ms added around each trial's start and end. |
| `uselog` | `false` | **Deuteron only.** `false` = use `Event_File_Reader_9_0.exe` (default). `true` = parse `logevents.txt` directly; use when events were logged but not transmitted to the recording system. |

### Motion sensors

| Field | Default | Description |
|---|---|---|
| `GetMotionSensors` | `false` | Extract head-direction / accelerometer data from Deuteron AUX channels or INTAN AUX*.dat files. |

### Data preprocessing

| Field | Default | Description |
|---|---|---|
| `lowpass` | `9000` | High-cut frequency (Hz) for the spike-band low-pass filter. `[]` = off. Must be < 9500. |
| `lowpassFT` | `200` | Low-pass cutoff (Hz) for the FieldTrip LFP stream (Butterworth 4th order, twopass). |
| `highpass` | `[]` | High-pass cutoff (Hz). `[]` = off. **For Deuteron wideband (Kilosort input), set to 300 Hz** — no high-pass is applied to the raw DF1 data by default. |
| `linefilter` | `0` | Line-noise notch centre (Hz). `0` = off. Applies a ±2 Hz band-stop filter. |
| `CAR` | `0` | Common-average re-referencing. `0` = off. When on, uses median reference. For > 32 channels, applies per 32-channel bank. |
| `dwnsmplRate` | `[]` | LFP downsample target (Hz). `[]` = auto → 937.5 Hz (= 30000/32, integer factor from standard INTAN rate). Must be an integer divisor of the raw sample rate. |
| `noise` | `[]` | Reserved for future noise-rejection parameters. |

### Sorting and curation

| Field | Default | Description |
|---|---|---|
| `kilosort` | `1` | Boolean flag: `1` = run Kilosort 4; `0` = skip spike sorting. Only Kilosort 4 is supported. |
| `KSchanMapFile` | `''` | Channel map filename (`.mat`). `''` = linear array (no custom map). Otherwise, place the file in `analysisCode/` and give its name here (e.g. `'chanMap_P32.mat'`). |
| `bombcell` | `true` | Run Bombcell automatic QC on Kilosort output. |
| `phy` | `false` | Open Phy after sorting. **Blocks MATLAB until Phy is closed.** |

### NGL02 options (post-Phy)

| Field | Default | Description |
|---|---|---|
| `doSpikething` | `true` | Process single-unit / spike data. |
| `doLFPthing` | `true` | Process LFP data. |
| `offlineTrack` | `false` | Run offline video blob detection. |
| `FLIP` | `false` | Run vFLIP laminar power analysis. |
| `useTrack` | `false` | Index spiking against social-tracking events. |
| `trialparsed` | `false` | Load trial-parsed FT file instead of continuous. |
| `artifdet` | `false` | LFP artifact detection and rejection. |
| `spectrogram` | `false` | Run multitaper time-frequency analysis. |
| `neurDyn.do` | `false` | Run neural-dynamics analysis. |

---

## 7. Event Coding System

### Concept

Behavioural events (stimulus onset, reward delivery, animal responses) are timestamped within the electrophysiology time series by toggling digital output pins. The NGL pipeline reads these pin states and converts them to a list of (time, event-code) pairs in `EventRecord`.

### Pin encoding

Events are encoded as **4-bit binary words** (INTAN uses 16-pin words, giving codes 0–65535). The 4 least-significant pins carry the "reserved" event vocabulary (codes 0–15). Additional project-specific events use higher-order pins (codes ≥ 16).

The binary-to-decimal conversion uses **LSB-first** ordering (`binvec2dec`), consistent with both Deuteron and INTAN hardware.

### Reserved event codes (do not change)

| Name | Decimal | Binary (LSB first) | Meaning |
|---|---|---|---|
| `itiOn` | 0 | `[0 0 0 0]` | Trial start / inter-trial interval onset |
| `stimOn1` | 1 | `[1 0 0 0]` | Stimulus 1 onset (sample, INI, etc.) |
| `stimOn2` | 2 | `[0 1 0 0]` | Stimulus 2 onset (match, choice, etc.) |
| `bhv` | 3 | `[1 1 0 0]` | Behavioural response detected |
| `end1` | 4 | `[0 0 1 0]` | Trial end — omission |
| `oms1` | 5 | `[1 0 1 0]` | Omission to Stim 1 |
| `oms2` | 6 | `[0 1 1 0]` | Omission to Stim 2 |
| `rwd` | 7 | `[1 1 1 0]` | Reward delivered |
| `preIni` | 8 | `[0 0 0 1]` | Transition (meaningless to INTAN) |
| `tr1` | 9 | `[1 0 0 1]` | Treatment / block / phase 1 start |
| `end2` | 10 | `[0 1 0 1]` | Trial end — punishment |
| `pun` | 11 | `[1 1 0 1]` | Punishment delivered |
| `na1` | 12 | `[0 0 1 1]` | Transition (meaningless to INTAN) |
| `tr2` | 13 | `[1 0 1 1]` | Treatment / block / phase 2 end |
| `na2` | 14 | `[0 1 1 1]` | Transition (meaningless to INTAN) |
| `end3` | 15 | `[1 1 1 1]` | Trial end — reward |

> **The reserved block (0–15) must never be changed.** These codes are hard-wired in the hardware delivery scripts and in `reservedEvents()` inside `eventDefinitions.m`.

### Event sequence example (correct trial)

```
itiOn    [0 0 0 0]  → trial starts
stimOn1  [1 0 0 0]  → stimulus appears
bhv      [1 1 0 0]  → animal responds
rwd      [1 1 1 0]  → reward
end3     [1 1 1 1]  → trial ends
         [1 1 1 0]  ← step back toward zero (transition)
         [1 1 0 0]  ← step back
         [1 0 0 0]  ← step back
itiOn    [0 0 0 0]  → next trial starts
```

### Debouncing

`INTAN_ExtractEvents` applies a debounce window of **28 samples** (≈ 0.93 ms at 30 kHz) to avoid registering the same edge twice due to pin noise.

---

## 8. Customising `eventDefinitions.m`

Copy `configfiles/eventDefinitions.m` into your project's `analysisCode/` folder and edit only the **`ONLY MODIFY THIS TWO BLOCKS`** section.

### Adding project-specific INTAN events

```matlab
elseif strcmpi(format,'fileperch')
    % Codes 16–65535 are available. Use a new integer for each event.
    eventdef.MyStim   = 16;   % [0000100000000000]
    eventdef.MyChoice = 17;   % [1000100000000000]
end
```

The comment shows the 16-bit binary pattern for reference. The decimal value is what the pipeline actually uses; the binary is only for your hardware delivery script.

### Rules

- **Never** modify `reservedEvents()` or any code above the `ONLY MODIFY` comment.
- Use unique integer values starting from **16** for your own events.
- Name fields exactly as they are referenced in `opt.alignto` and `conditions_script.m`.
- Deuteron-format events go in the `strcmpi(format,'DF1')` block.

---

## 9. Pipeline Stages — What Happens When You Run

### Stage 00 — Validation (`set_default`)

`set_default(input, opt)` is called once and does everything needed before the session loop:

1. Loads all defaults from `default_opt()`.
2. Overlays user `opt` fields, warning on unknown names.
3. Validates types (e.g. `opt.alignto` must be a cell of chars, `opt.kilosort` must be 2 or 4).
4. Checks cross-option consistency (e.g. `doNWB` + H5 DLL conflict, `phy` without `bombcell`).
5. Builds the full IKN path structure from `datadrive` and `studyname`.
6. Reads `NGL_machineConfig` for Python env paths, toolbox root, and Deuteron EXE paths.
7. Resolves the subject list (`'all'` or a cell of names).
8. Adds all dependency folders to MATLAB path and calls `ft_defaults`.

### Stage 01 — Find sessions (`findSessions`)

Scans `data/raw/<subject>/<date>/` for each requested subject and date, building `input.sessions(s).list{}` of valid session folders.

### Stage 02 — Prepare session (`prepforsession`)

For each `(subject, session)` pair:

1. Changes directory to the raw data folder.
2. Calls `chckV()` to detect the recording format (`fileperch`, `filepertype`, `tradFormat`, `DT2`, `DF1`, or `FieldTrip`).
3. Sets all session-specific path fields in `opt` (e.g. `opt.FolderProcDataMat`, `opt.KSfolder`, `opt.trialSorted`).
4. Creates missing session-specific output folders.
5. For INTAN formats, calls `findSetting` to read the RHD header and `settings.xml` for sample rate and channel count.
6. Checks actual `.dat` file count against `opt.numChannels`; runs `reduceChanMap` if fewer channels are present than expected.

### Stage 03 — Format-specific pipeline

The switch dispatches on `info.fileformat`:

#### INTAN (`fileperch`, `filepertype`, `tradFormat`) → `INTAN_PipelineWrapper`

1. **Event extraction** (`EventProcess` → `INTAN_ExtractEvents`): reads all `board-DIGITAL-IN*.dat` files, detects pin transitions, applies debounce, converts bit vectors to decimal codes, calls `trialdefGen` to produce trial boundaries in ms, runs `conditions_script`.
2. **NWB export** (`intan2NWB_neuroconv`): only if `opt.doNWB = true`. Requires NeuroConv Python env.
3. **Kilosort `.bin`** (`Intan2Kilosort_wrapper`): reads amp\*.dat files, scales to µV (× 0.195), optionally detrends + high-pass + low-pass (for the spike band), writes a flat binary interleaved `.bin` in chunks of 300 s.
4. **FieldTrip LFP** (`intan2MAT_wrapper` + `MAT2FieldTrip`): only if `opt.FieldTrip = true`. Reads amp\*.dat, scales, optionally CAR + low-pass + line-filter, downsamples to 937.5 Hz, builds a pseudo-FieldTrip struct, calls `ft_redefinetrial` to produce continuous and/or trial-parsed `.mat` files.
5. **Motion sensors** (`GetMotionSensors`): only if `opt.GetMotionSensors = true`. Reads AUX\*.dat accelerometer channels.

#### Deuteron (`DT2`, `DF1`) → `Deuteron_PipelineWrapper`

1. **Metadata** (`Deuteron_GetMetaData`): reads hardware constants (32 kHz sample rate, 16-bit ADC, 1.95 × 10⁻⁷ V/bit for DF1). Channel count is **not** available in the file header — it is inferred from the event log in the next step.

2. **Event extraction** (`EventProcess` → `Deuteron_ExtractEvents`): by default (`opt.uselog = false`) invokes `Event_File_Reader_9_0.exe` on all `NEUR*.DF1` files to produce an `EventRecord.CSV`. The CSV is parsed into the `EventRecord` struct. The channel-mapping entry in the CSV is used to set `opt.channelOrder` and `opt.numChannels` — these are critical for correct data reshaping downstream. If the EXE path is unavailable or events were only logged (not transmitted), set `opt.uselog = true` to parse `logevents.txt` instead. After extraction, `trialdefGen` and `conditions_script` run identically to the INTAN path.

3. **Wideband `.bin`** (`Deuteron2Kilosort`): only if `opt.bin = true`. Reads each `NEUR*.DF1` file block by block via `Deuteron_extractData` (stream = 1), converts raw ADC samples to int16 µV using `int16((voltageResolution × (data − offset)) × 1e6)`, optionally applies CAR, detrends each channel, applies a high-pass filter if `opt.highpass` is set (recommended: 300 Hz for Kilosort 4), and appends to a flat interleaved `.bin` file in channels × samples layout. Skips if a non-empty `.bin` already exists.

4. **FieldTrip LFP** (`Deuteron2Fieldtrip` + `MAT2FieldTrip`): only if `opt.FieldTrip = true`. Reads the same `NEUR*.DF1` files, converts to µV, optionally applies CAR, then filters channel by channel (detrend → low-pass at `opt.lowpassFT` = 200 Hz → optional line-noise band-stop), and downsamples to `sampleRate / 32` = **1000 Hz** (Deuteron's integer factor, analogous to INTAN's 937.5 Hz). Intermediate output is cached as `<session>_filt_dwn.mat` to avoid reprocessing on re-runs. The pseudo-FieldTrip struct is then passed to `MAT2FieldTrip`, which produces `_FTcont.mat` and per-alignment-event `.mat` files identically to the INTAN path.

5. **Motion sensors** (`GetMotionSensors` → `getfrom_Deuteron`): only if `opt.GetMotionSensors = true`. Reads the motion-sensor stream (stream = 2) from each `NEUR*.DF1` file via `Deuteron_extractData`, extracting Accelerometer, Gyroscope, and Magnetometer data from the embedded MPU-9250 blocks at 1000 Hz. Applies `magcal` for hard/soft-iron magnetometer correction, runs an AHRS filter (`Deuteron_estimateheading`) to estimate heading and dead-reckoning position, plots raw sensor timeseries (`Deuteron_PlotMotionSensors`), and saves `MotionData.mat`.

**Key differences from INTAN:**

| Property | INTAN | Deuteron |
|---|---|---|
| Raw sample rate | 30 000 Hz | 32 000 Hz |
| LFP downsample target | 937.5 Hz (30000/32) | 1000 Hz (32000/32) |
| Channel count source | RHD header / `.dat` file count | Event log (`opt.channelOrder`) |
| Event source | `board-DIGITAL-IN*.dat` pin transitions | `NEUR*.DF1` via EXE or `logevents.txt` |
| Pin state convention | starts at `[0 0 0 0]` | starts at `[1 1 0 0]`, cumulative tracking |
| File layout | one `.dat` per channel (flat) | block format, all channels per file |
| Wideband creator | `Intan2Kilosort_wrapper` | `Deuteron2Kilosort` |
| LFP creator | `intan2MAT_wrapper` | `Deuteron2Fieldtrip` |
| Motion sensors | AUX*.dat (3-axis accelerometer) | MPU-9250 (accel + gyro + magnetometer) |

### Stage 04 — Spike sorting (`master_kilosort4`)

1. Activates the KS4 Python environment via `pyenv`.
2. Copies `parameters.py` and `master_kilosort4.py` from the toolbox into the KS Python environment.
3. Calls `pyrunfile` with the `.bin` path, channel count, and channel map path.
4. Kilosort writes its output to `preprocessing/<subject>/<session>/kilosort/4/`.

### Stage 05 — Quality metrics (`Bombcell_Main`)

Only if `opt.bombcell = true`. Loads Kilosort output from the KS folder and runs `bc.qm.runAllQualityMetrics`. Parameters can be overridden in `bombcellConfig.m` (stored in `analysisCode/`).

### Stage 06 — Manual curation (`Phy`)

Only if `opt.phy = true`. Changes directory to `opt.FolderProcDataMat` and calls `system('phy template-gui params.py')`. **This blocks MATLAB until Phy is closed.** Manual curation is typically done as a separate step, after reviewing Bombcell results.

---

## 10. Outputs

| File | Location | Created by | Notes |
|---|---|---|---|
| `<session>.bin` | `preprocessing/<subj>/<session>/` | `Intan2Kilosort_wrapper` (INTAN) · `Deuteron2Kilosort` (Deuteron) | Flat int16 interleaved binary for KS4 |
| `kilosort/4/` folder | `preprocessing/<subj>/<session>/` | `master_kilosort4` | KS4 templates, spike times, amplitudes |
| `EventRecord.CSV` | `preprocessing/<subj>/<session>/` | `Deuteron_ExtractEvents` (EXE path) | Raw Deuteron event log generated by Event_File_Reader_9_0.exe |
| `EventRecord.mat` | `preprocessing/<subj>/<session>/` | `EventProcess` | Parsed event struct (both INTAN and Deuteron) |
| `trialdef.mat` | `trialSorted/<subj>/<session>/` | `trialdefGen` | Trial boundaries in ms per alignment event |
| `events.mat` | `trialSorted/<subj>/<session>/` | `trialdefGen` | Trial-aligned event struct |
| `condition.mat` | `trialSorted/<subj>/<session>/` | `conditions_script` | Condition grouping |
| `<session>_FTcont.mat` | `trialSorted/<subj>/<session>/` | `MAT2FieldTrip` | Continuous FieldTrip LFP data |
| `<session>_<event>.mat` | `trialSorted/<subj>/<session>/` | `MAT2FieldTrip` | Trial-parsed FieldTrip LFP per alignment eve