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
13. [Multi-Area Mode](#13-multi-area-mode)
14. [NWB Export](#14-nwb-export)
15. [Contributors](#15-contributors)
16. [External Tools & Citations](#16-external-tools--citations)

---

## 1. Overview

The NGL pipeline converts raw electrophysiology recordings — acquired with **INTAN RHX** or **Deuteron** hardware — into analysis-ready formats:

- **`.bin` file** for spike sorting with [Kilosort 4](https://github.com/MouseLand/Kilosort)
- **FieldTrip `.mat` structure** for LFP analysis
- **`EventRecord.mat` / `trialdef.mat` / `events.mat`** for trial-parsed alignment
- **Bombcell QC report** for semi-automatic cluster curation
- **NWB file** for data sharing and archiving

The pipeline is driven by a single user-facing script, **`NGL_SetAndRunMe.m`**, which lives in your project's `analysisCode` folder. You only ever need to edit that one file.

```
NGL_SetAndRunMe  -  NGL00_Prep  -  NGL01_Main  -  [Phy]  -  NGL02_postPhy
```

---

## 2. Prerequisites

### MATLAB toolboxes (bundled under `toolboxes\`)

| Toolbox | Purpose |
|---|---|
| `fieldtrip_light` | LFP preprocessing and FieldTrip format |
| `Intan` | INTAN RHD/RHS header reader |
| `Deuteron` | Deuteron binary reader + event DLL |
| `bombcell` | Automatic spike-cluster QC |
| `npy-matlab` | NumPy `.npy` reader (KS4 output) |
| `prettify_matlab` | Plot formatting |
| `spikes` | Spike analysis utilities |
| `CADopti` | Cell assembly detection |
| `matnwb` | NWB MATLAB interface |
| `BDPAT_NGL` | NGL-specific analysis helpers |
| `Viewer` | Data viewer |

All paths are added automatically by `set_default` on each run — **do not add them manually** to MATLAB's persistent path. The `functions\` folder is partitioned into logical subfolders (`pipeline\`, `intan\`, `deuteron\`, `sorting\`, `events\`, `analysis\`, `video\`, `plotting\`, `ethology\`, `utils\`); each is added by name, deliberately excluding `_deprecated\`.

### Python environments

Three local Conda environments are required for full functionality:

| Environment | Purpose | Config key |
|---|---|---|
| KS4 Python env | Kilosort 4 spike sorting | `cfg.KSpythonExe` |
| PHY Python env | Phy manual curation | `cfg.PHYpythonExe` |
| NeuroConv env | NWB export (`opt.doNWB = true`) | `cfg.NCpythonExe` |

Look up each project's documentation to find out how to install them on your own machine.

### Project config files (in `analysisCode\`)

`set_default` calls `checkAnalysisCode(input, opt)` automatically and raises a single error listing every missing required file before the session loop begins.

| File | Required when | Notes |
|---|---|---|
| `chanMap*.mat` | Always | Channel map for Kilosort |
| `eventDefinitions.m` | `opt.RetrieveEvents` | Copy from `configfiles\`, add project events |
| `master_kilosort4.py` | `opt.kilosort` | Auto-copied from toolbox at runtime |
| `parameters.py` | `opt.kilosort` | KS4 parameters - tune per project/probe |
| `bombcellConfig.m` | `opt.bombcell` | Bombcell QC thresholds |
| `master_neuroconv.py` | `opt.doNWB` | Auto-copied from toolbox at runtime |
| `nwb_metadata.yaml` | `opt.doNWB` | Copy `nwb_metadata_template.yaml`, fill in and change name|
| `conditions_script.m` | Recommended | Condition grouping; warning if absent |
| `postPhy_param.m` | Recommended | NGL02 parameters; warning if absent |

---

## 3. Folder Structure (IKN Standard)

`NGL00_Prep` creates the full folder tree on first run. The expected layout is:

```
<datadrive>:\<studyname>\
│
├── analysisCode\               - NGL_SetAndRunMe, eventDefinitions, chanMaps, etc.
│
├── data\
│   ├── raw\
│   │   └── <subject>\
│   │       └── <DDMMYYYY>\    - raw INTAN .dat / Deuteron .DT2 / .DF1 files
│   │
│   ├── preprocessing\
│   │   └── <subject>\
│   │       └── <DDMMYYYY>\    - .bin, kilosort output, EventRecord.mat, .nwb
│   │
│   ├── spikeSorted\
│   │   └── <subject>\
│   │       └── <DDMMYYYY>\    - post-Phy spike data
│   │
│   ├── trialSorted\
│   │   └── <subject>\
│   │       └── <DDMMYYYY>\    - trialdef.mat, events.mat, FT .mat files
│   │
│   ├── behaviour\
│   │   └── <subject>\
│   │       └── <DDMMYYYY>\    - behavioural files
│   │
│   └── analysis\
│       └── <subject>\
│           └── <DDMMYYYY>\    - analysis outputs
│
├── manuscript\				   - publications related to the data
├── paradigmCode\			   - the final paradigm code used during the project
└── training\ 				   - relevant data related to animal training
```

> **Session folder naming:** Raw data folders are expected to use the format `DDMMYYYY` (day-month-year). This convention is used by the NWB fallback timestamp parser when session start time is absent from the recording file header.

> **IKN standard reference:** `gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure`

---

## 4. Machine Configuration (`NGL_machineConfig`)

`NGL_machineConfig.m` must exist at the toolbox root. Returns a struct with machine-specific absolute paths. This is a template, fill in your own paths:

```matlab
function cfg = NGL_machineConfig() 						% located toolbox root
    cfg.toolbox      = 'C:\Code\ephys-data-pipeline';   % toolbox root
    cfg.KSpythonExe  = 'C:\...\envs\kilosort4\';            % KS4 conda env
    cfg.PHYpythonExe = 'C:\...\envs\phy2\';                 % Phy conda env
    cfg.NCpythonExe  = 'C:\...\envs\neuroconv\';            % NeuroConv env
end
```

> `NGL_machineConfig` is a single file that exists **per computer**.

---

## 5. User Configuration (`NGL_SetAndRunMe`)

`NGL_SetAndRunMe.m` is the **only file you need to edit** for each project or recording batch.

### Step-by-step walkthrough

> All are **examples**, not necessarily your optimal choices

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
dates    = {'YYYYMMDD', 'YYYYMMDD'};         % or 'all' for every session

opt = struct();
opt.numChannels   = 32;                      % electrode count
opt.KSchanMapFile = 'chanMap_*.mat';         % '' (empty) for a simple linear array
opt.alignto       = {'itiOn', 'stimOn1'};    % events to align LFP trials to
opt.FieldTrip     = false;                   % don't produce LFP .mat file
opt.bombcell      = true;                    % run QC after sorting

NGL01_Main
```

**Multi-area** — additionally set:

```matlab
input.Areas = {'NCL', 'NCL', 'STR'};        % one label per kcoords group in chanMap. In this case, a 2-shank probe is in NCL and a second 1-shank probe in STR
```

**NWB export** — additionally set:

```matlab
opt.doNWB = true;                            % requires NeuroConv env + nwb_metadata.yaml
```

**Section 3 — Run post-Phy**

After manual Phy curation:

```matlab
postPhy_param();
NGL02_postPhy
```

---

## 6. All Options Reference (`opt` fields)

Every field must be set in `default_opt.m`. If you set a field in `NGL_SetAndRunMe` that does not appear in `default_opt`, `set_default` will issue a warning and ignore it.

### Data format

| Field | Default | Description |
|---|---|---|
| `numChannels` | `32` | Expected electrode channel count. For INTAN, the actual count is read from the header, but set this to match.|
| `bin` | `true` | Create the `.bin` file for Kilosort. Set `false` to skip for any reason. |
| `FieldTrip` | `true` | Produce a FieldTrip-ready `.mat` file from the LFP stream. |
| `doNWB` | `true` | Run INTAN→NWB conversion via NeuroConv. Requires `cfg.NCpythonExe` and `nwb_metadata.yaml` in `analysisCode\`. |

### Events

| Field | Default | Description |
|---|---|---|
| `RetrieveEvents` | `true` | Extract event log from the session. |
| `alignto` | `{'itiOn'}` | Event(s) to use as trial time-zero. Single char or cell array of chars matching names in `eventDefinitions`. At minimum keep `'itiOn'`. |
| `trEvents` | `{}` | "Special" ITI events (e.g. drug treatments, tutor calls). Cell array of event names. |
| `addtime` | `0` | Extra padding in ms added around each trial's start and end. |
| `uselog` | `false` | **Deuteron only.** `false` = use `Event_File_Reader_9_0.exe` (default). `true` = parse `logevents.txt` directly. |

### Motion sensors

| Field | Default | Description |
|---|---|---|
| `GetMotionSensors` | `false` | Extract head-direction / accelerometer data from Deuteron or INTAN AUX dedicated channels. |

### Data preprocessing

| Field | Default | Description |
|---|---|---|
| `lowpass` | `9000` | High-cut frequency (Hz) for the spike-band low-pass filter. `[]` = off. Must be < 9500. |
| `lowpassFT` | `200` | Low-pass cutoff (Hz) for the FieldTrip LFP stream (Butterworth 4th order, twopass). |
| `highpass` | `[]` | High-pass cutoff (Hz). `[]` = off. **For Deuteron wideband (Kilosort input), set to 300 Hz.** |
| `linefilter` | `0` | Line-noise notch centre (Hz). `0` = off. Applies a ±2 Hz band-stop filter. |
| `CAR` | `0` | Common-average re-referencing. `0` = off. When on, uses median reference. For > 32 channels, applies per 32-channel bank. |
| `dwnsmplRate` | `[]` | LFP downsample target (Hz). `[]` = auto -> 937.5 Hz (INTAN) or 1000 Hz (Deuteron). |
| `noise` | `[]` | Reserved for future noise-rejection parameters. |

### Sorting and curation

| Field | Default | Description |
|---|---|---|
| `kilosort` | `true` | logic flag: `true` = run Kilosort 4; `false` = skip. Only Kilosort 4 is supported. |
| `KSchanMapFile` | `''` | Channel map filename (`.mat`). `''` = linear array. Place the file in `analysisCode\`. |
| `bombcell` | `true` | Run Bombcell automatic QC on Kilosort output. |
| `phy` | `false` | Open Phy after sorting. **Blocks MATLAB until Phy is closed.** |

### NGL02 options (post-Phy)

| Field | Default | Description |
|---|---|---|
| `doSpikething` | `true` | Process single-unit / spike data. |
| `doLFPthing` | `true` | Process LFP data. |
| `offlineTrack` | `false` | Run offline video blob detection. |
| `FLIP` | `false` | Run vFLIP laminar power analysis. TESTING |
| `useTrack` | `false` | Index spiking against social-tracking events. |
| `trialparsed` | `false` | Load trial-parsed FT file instead of continuous. |
| `artifdet` | `false` | LFP artifact detection and rejection. |
| `spectrogram` | `false` | Run multitaper time-frequency analysis. |
| `neurDyn.do` | `false` | Run neural-dynamics analysis. TESTING|

---

## 7. Event Coding System

### Concept

Behavioural events (stimulus onset, reward delivery, animal responses) are timestamped within the electrophysiology time series by toggling digital output pins. The NGL pipeline reads these pin states and converts them to a list of (time, event-code) pairs in `EventRecord`.

### Pin encoding

Events are encoded as **N-bit binary words** (INTAN uses 16-pin words, Deuteron uses 4-pin for now). The 4 least-significant pins carry the "reserved" event vocabulary (codes 0–15). Additional project-specific events use higher-order pins (codes ≥ 16).

The binary-to-decimal conversion uses **LSB-first** ordering (`binvec2dec`), consistent with both Deuteron and INTAN hardware.

### Reserved event codes (NOT to be changed)

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

### Event sequence example (for a 'correct' trial)

```
itiOn    [0 0 0 0]  - trial starts
stimOn1  [1 0 0 0]  - stimulus appears
bhv      [1 1 0 0]  - animal responds
rwd      [1 1 1 0]  - reward
end3     [1 1 1 1]  - trial ends
         [1 1 1 0]  - step back toward zero (transition)
         [1 1 0 0]  - step back
         [1 0 0 0]  - step back
itiOn    [0 0 0 0]  - next trial starts
```

### Debouncing

`INTAN_ExtractEvents` applies a debounce window of **28 samples** (≈ 0.93 ms at 30 kHz) to avoid registering the same edge twice due to pin noise.

---

## 8. Customising `eventDefinitions.m`

Copy `configfiles\eventDefinitions.m` into your project's `analysisCode\` folder and edit only the **`ONLY MODIFY THIS TWO BLOCKS`** section.

### Adding project-specific INTAN events

```matlab
elseif strcmpi(format,'fileperch')
    % Codes 16–65535 are available. Use a new integer for each event.
    eventdef.MyStim   = 16;   % [0000100000000000]
    eventdef.MyChoice = 17;   % [1000100000000000]
end
```

The comment shows the 16-bit binary pattern for reference. The decimal value is what the pipeline actually uses.

### Rules

- **Never** modify `reservedEvents()` or any code above the `ONLY MODIFY` comment.
- Use unique integer values starting from **16** for your own events.
- Name fields exactly as they are referenced in `opt.alignto` and `conditions_script.m`.
- Deuteron-format events go in the `strcmpi(format,'DF1')` block.

---

## 9. Pipeline Stages

### Stage 00 - Validation and initialisation (`set_default`)

`set_default(input, opt)` is called once before the session loop and does everything needed for a safe run:

1. Loads all defaults from `default_opt()`.
2. Overlays user `opt` fields; issues a warning for any unrecognised field name.
3. Validates field types (e.g. `opt.alignto` must be a cell of chars, `opt.numChannels` must be a positive scalar).
4. Checks cross-option consistency (e.g. `opt.phy = true` without `opt.bombcell = true`).
5. Builds the full IKN path structure from `datadrive` and `studyname`.
6. Reads `NGL_machineConfig` for Python env paths, toolbox root, and Deuteron EXE paths.
7. Resolves the subject list (`'all'` or a cell of names).
8. Adds all dependency folders to MATLAB path (explicit per-subfolder `addpath`, `_deprecated\` excluded) and calls `ft_defaults`.
9. Calls `checkAnalysisCode(input, opt)` — validates all required config files exist in `analysisCode\`; raises a descriptive error listing every missing file if any are absent.
10. If `input.Areas` is set, calls `buildAreaMap` to generate per-area channel masks and write per-area channel map `.mat` files to `analysisCode\`.

### Stage 01 - Find sessions (`findSessions`)

Scans `data\raw\<subject>\<date>\` for each requested subject and date, building `input.sessions(s).list{}` of valid session folders.

### Stage 02 - Prepare session (`prepforsession`)

For each `(subject, session)` pair:

1. Changes directory to the raw data folder.
2. Calls `chckV()` to detect the recording format (`fileperch`, `filepertype`, `tradFormat`, `DT2`, `DF1`, or `FieldTrip`).
3. Sets all session-specific path fields in `opt` (e.g. `opt.FolderProcDataMat`, `opt.KSfolder`, `opt.trialSorted`).
4. Creates missing session-specific output folders.
5. For INTAN formats, calls `findSetting` to read the RHD header and `settings.xml` for sample rate and channel count.
6. Checks actual `.dat` file count against `opt.numChannels`; runs `reduceChanMap` if fewer channels are present than expected.

### Stage 03 - Format-specific pipeline

The switch dispatches on `info.fileformat`:

#### INTAN (`fileperch`, `filepertype`, `tradFormat`) - `INTAN_PipelineWrapper`

1. **Event extraction** (`EventProcess` -> `INTAN_ExtractEvents`): reads all `board-DIGITAL-IN*.dat` files, detects pin transitions, applies debounce, converts bit vectors to decimal codes, calls `trialdefGen` to produce trial boundaries in ms, runs `conditions_script`.
2. **NWB export** (`intan2NWB_neuroconv`): only if `opt.doNWB = true`. See [§14](#14-nwb-export).
3. **Kilosort `.bin`** (`Intan2Kilosort_wrapper`): reads amp\*.dat files, scales to µV (× 0.195), optionally detrends + high-pass + low-pass, writes a flat binary interleaved `.bin` in chunks of 300 s.
4. **FieldTrip LFP** (`intan2MAT_wrapper` + `MAT2FieldTrip`): only if `opt.FieldTrip = true`. Reads amp\*.dat, scales, optionally CAR + low-pass + line-filter, downsamples to 937.5 Hz, builds a pseudo-FieldTrip struct, calls `ft_redefinetrial` to produce continuous and/or trial-parsed `.mat` files.
5. **Motion sensors** (`GetMotionSensors`): only if `opt.GetMotionSensors = true`. Reads AUX\*.dat accelerometer channels.

#### Deuteron (`DT2`, `DF1`) - `Deuteron_PipelineWrapper`

1. **Metadata** (`Deuteron_GetMetaData`): reads hardware constants (32 kHz, 16-bit ADC, 1.95 × 10⁻⁷ V/bit for DF1). Channel count is **not** in the file header — it is inferred from the event log.

2. **Event extraction** (`EventProcess` -> `Deuteron_ExtractEvents`): by default invokes `Event_File_Reader_9_0.exe` on all `NEUR*.DF1` files to produce an `EventRecord.CSV`. The channel-mapping entry in the CSV is used to set `opt.channelOrder` and `opt.numChannels`. If the EXE is unavailable or events were only logged, set `opt.uselog = true` to parse `logevents.txt` instead.

3. **Wideband `.bin`** (`Deuteron2Kilosort`): reads each `NEUR*.DF1` file block by block, converts raw ADC to int16 µV, optionally applies CAR, detrends each channel, applies high-pass filter if `opt.highpass` is set (recommended: 300 Hz for KS4), appends to a flat interleaved `.bin`.

4. **FieldTrip LFP** (`Deuteron2Fieldtrip` + `MAT2FieldTrip`): reads `.DF1` files, converts to µV, optionally applies CAR, filters channel by channel (detrend → low-pass at `opt.lowpassFT` → optional line notch), downsamples to 1000 Hz (32000/32). Intermediate output cached as `<session>_filt_dwn.mat`.

5. **Motion sensors** (`GetMotionSensors` -> `getfrom_Deuteron`): reads MPU-9250 motion stream from `.DF1` files at 1000 Hz. Applies `magcal` for magnetometer correction, runs AHRS filter (`Deuteron_estimateheading`) for heading + dead-reckoning, saves `MotionData.mat`.

**Key differences between INTAN and Deuteron:**

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

### Stage 04 - Spike sorting (`master_kilosort4`)

1. Activates the KS4 Python environment via `pyenv`.
2. Copies `parameters.py` and `master_kilosort4.py` from `analysisCode\` to the KS Python environment folder.
3. Calls `pyrunfile` with the `.bin` path, channel count, and channel map path.
4. Kilosort writes its output to `preprocessing\<subject>\<session>\kilosort4\` (single-area) or `preprocessing\<subject>\<session>\<areaLabel>\` (multi-area).

In multi-area mode, Stage 04 loops over `opt.KSfolders` (one entry per area) and runs a separate Kilosort job per area.

### Stage 05 - Quality metrics (`Bombcell_Main`)

Only if `opt.bombcell = true`. Locates the Kilosort output via `opt.KSfolder` and runs `bc.qm.runAllQualityMetrics`. Parameters are read from `bombcellConfig.m` in `analysisCode\`.

**Area tagging** (multi-area mode only): after Bombcell completes, if `opt.KSfolders` is present (the multi-area guard), the area label is written to three complementary locations:

- **`qMetric.area`** - cell array of area labels added to `qMetrics.mat`. Enables `strcmp(qMetric.area, 'NCL')` filtering in downstream MATLAB analysis.
- **`cluster_area.tsv`** - written to the KS output folder. Phy automatically loads any `cluster_*.tsv` and displays it as an extra column in the cluster table. No Phy configuration is needed.
- **`area_label.txt`** - plain-text file alongside the Bombcell `.npy` output files, for Python analysis scripts.

### Stage 06 - Manual curation (`Phy`)

Only if `opt.phy = true`. Changes directory to `opt.FolderProcDataMat` and calls `system('phy template-gui params.py')`. **This blocks MATLAB until Phy is closed.** Manual curation is typically done as a separate step after reviewing Bombcell results. In multi-area mode, Phy is opened once per area.

---

## 10. Outputs

| File | Location | Created by | Notes |
|---|---|---|---|
| `<session>.bin` | `preprocessing\<subj>\<session>\` | `Intan2Kilosort_wrapper` · `Deuteron2Kilosort` | Flat int16 interleaved binary for KS4 |
| `kilosort4\` (single-area) | `preprocessing\<subj>\<session>\` | `master_kilosort4` | KS4 templates, spike times, amplitudes |
| `\<area>\` (multi-area) | `preprocessing\<subj>\<session>\` | `master_kilosort4` | One KS4 output folder per area |
| `cluster_area.tsv` | `\<area>\bombcell\` | `Bombcell_Main` | Area label per cluster; auto-loaded by Phy |
| `qMetrics.mat` | `preprocessing\<subj>\<session>\kilosort\bombcell\qMetrics\` | `Bombcell_Main` | Bombcell quality metrics + `qMetric.area` field |
| `area_label.txt` | `preprocessing\<subj>\<session>\kilosort\bombcell\` | `Bombcell_Main` | Plain-text area label for Python scripts |
| `<session>.nwb` | `preprocessing\<subj>\<session>\` | `intan2NWB_neuroconv` | NWB file with raw ephys; only if `doNWB=true` |
| `EventRecord.CSV` | `preprocessing\<subj>\<session>\` | `Deuteron_ExtractEvents` (EXE) | Raw Deuteron event log |
| `EventRecord.mat` | `preprocessing\<subj>\<session>\` | `EventProcess` | Parsed event struct (INTAN and Deuteron) |
| `trialdef.mat` | `trialSorted\<subj>\<session>\` | `trialdefGen` | Trial boundaries in ms per alignment event |
| `events.mat` | `trialSorted\<subj>\<session>\` | `trialdefGen` | Trial-aligned event struct |
| `condition.mat` | `trialSorted\<subj>\<session>\` | `conditions_script` | Condition grouping |
| `<session>_FTcont.mat` | `preprocessing\<subj>\<session>\` | `MAT2FieldTrip` | Continuous FieldTrip LFP data |
| `<session>_<event>.mat` | `trialSorted\<subj>\<session>\` | `MAT2FieldTrip` | Trial-parsed FieldTrip LFP per alignment event |
| `MotionData.mat` | `trialSorted\<subj>\<session>\` | `GetMotionSensors` | Heading, accelerometer, dead-reckoning position |

---

## 11. Running Post-Phy Processing (`NGL02_postPhy`)

After completing manual curation in Phy, run `NGL02_postPhy` to build analysis-ready spike and LFP variables:

1. Copy `postPhy_param.m` from `configfiles/` to `analysisCode/` if you haven't already.
2. Edit `postPhy_param.m` for your project-specific parameters.
3. In `NGL_SetAndRunMe.m`, Section 3:

```matlab
postPhy_param();           % loads NGL02 parameters into workspace
NGL02_postPhy              % runs post-Phy processing
```

`NGL02_postPhy` reads the curated Phy output (`.tsv` cluster labels), loads `qMetrics.mat`, and builds per-unit spike time vectors, spike waveforms, and trial-parsed spike data into `spikeSorted\`. It also produces trial-parsed LFP structures in `trialSorted/` if `opt.doLFPthing = true`.

---

## 12. Troubleshooting

**`checkAnalysisCode` raises an error listing missing files.**
All listed files must exist in `analysisCode\` before `NGL01_Main` proceeds. Copy the relevant templates from `configfiles/` and fill in the project-specific values.

**Kilosort fails to find Python.**
Check `cfg.KSpythonExe` in `NGL_machineConfig.m`. The path must point to the Conda environment root (the folder containing `python.exe`), not to `python.exe` itself.

**`session_start_time` warning in NeuroConv output.**
This is expected for INTAN recordings — the `.rhd` format does not store an absolute timestamp. The pipeline parses the date from the session folder name (DDMMYYYY convention). Check the folder name matches the convention if the warning shows a wrong date.

**Deuteron channel count is wrong.**
`opt.numChannels` is set from the event log CSV, not from the `.DF1` header. If the channel count looks wrong, check the event log file for the channel-mapping line, or confirm that `Event_File_Reader_9_0.exe` ran successfully.

**Phy does not show the `area` column.**
Confirm that `cluster_area.tsv` exists in the KS output folder for the session and that it is tab-delimited. Phy loads any `cluster_*.tsv` automatically; no Phy configuration is required.

**`opt.lowpass` validation error (`must be < 9500`).**
The pipeline enforces a hard upper limit of 9499 Hz to stay below the Nyquist frequency for 30 kHz recordings. Set a value below that or use `[]` to disable.

**Files processed by a previous run are being overwritten.**
Set `opt.bin = false` to skip `.bin` re-creation if the file already exists. Kilosort, Bombcell, and FieldTrip always re-run unless you comment out the relevant blocks in `NGL01_Main`.

---

## 13. Multi-Area Mode

Multi-area mode activates when `input.Areas` is defined in `NGL_SetAndRunMe`:

```matlab
input.Areas = {'NCL', 'NCL', 'STR'};
```

Each entry in `input.Areas` labels the `kcoords` group at the same index in the channel map. Entries can repeat (when multiple shanks record from the same area). The pipeline uses these labels to:

1. **Build per-area channel maps** - `buildAreaMap` reads `opt.KSchanMapFile` from `analysisCode/`, extracts the set of channels belonging to each unique area, and writes a separate `chanMap_<area>.mat` per area back to `analysisCode\`. It also returns `input.areaMap`, a struct with fields `uniqueAreas`, `channelMasks`, and `areaLabels`.

2. **Sort each area independently** - Kilosort is called once per unique area, using the per-area channel map. Output goes to `kilosort\<areaLabel>\` inside the session preprocessing folder.

3. **Run Bombcell per area** - `opt.KSfolder` is set to the per-area Kilosort output folder before each Bombcell run. After `runAllQualityMetrics` completes, three area-tagging outputs are written (see Stage 05).

4. **Tag clusters for downstream analysis** - `qMetric.area` allows simple MATLAB filtering, e.g.:

```matlab
load('qMetrics.mat')
nclUnits = find(strcmp(qMetric.area, 'NCL') & qMetric.label == 1);
```

In Phy, the `area` column appears automatically in the cluster table from `cluster_area.tsv`. Python scripts can read `area_label.txt` from the Bombcell output folder.

### Requirements for multi-area mode

- `opt.KSchanMapFile` must be set (a linear map without `kcoords` groups is not supported).
- The `kcoords` array in the channel map must have one integer per channel, with values matching the order of unique areas in `input.Areas`.
- `opt.bombcell = true` is recommended; area labels are only written when Bombcell runs.

---

## 14. NWB Export

NWB (Neurodata Without Borders) export converts the raw INTAN recording to a standardised HDF5 file for data sharing and archiving. It is enabled by setting `opt.doNWB = true`.

### How it works

`intan2NWB_neuroconv.m` is called from `INTAN_PipelineWrapper` after event extraction. It:

1. Checks whether an NWB file already exists for the session (skip if it does and is larger than 1 MB, indicating a completed conversion).
2. Locates the `info.rhd` header file in the raw data folder.
3. Calls `pyrunfile` with `master_neuroconv.py` (which is auto-copied from the toolbox to the NeuroConv Python environment folder at runtime).
4. `master_neuroconv.py` creates an `IntanRecordingInterface`, loads `nwb_metadata.yaml` from `analysisCode/`, merges it with auto-extracted metadata, resolves the session start time (from the `.rhd` header if present; otherwise parsed from the DDMMYYYY session folder name), and calls `interface.run_conversion`.
5. The `.nwb` file is written to `preprocessing/<subj>/<session>/<session>.nwb` with `overwrite=True`.

Future pipeline stages (spike sorting results, behavioural data) will append to the same NWB file using `overwrite=False`.

### Project metadata YAML

Copy `configfiles/nwb_metadata_template.yaml` to `analysisCode/nwb_metadata.yaml` and fill in the project-stable fields:

```yaml
NWBFile:
  session_description: "Spatial working memory task in freely-moving birds"
  lab: "Neural Basis of Learning"
  institution: "Ruhr-University Bochum"
  experimenter:
    - "Fulano, Mengano"
  experiment_description: "In vivo extracellular recordings during spatial working memory"
  keywords:
    - "electrophysiology"
    - "working memory"

Subject:
  species: "Columba livia"
  sex: "U"
  age: "P5Y"       # ISO 8601 duration; approximate age at study start
```

Fields left absent from the YAML (e.g. `session_start_time`, `session_id`) are either extracted automatically from the recording file or set programmatically by `master_neuroconv.py`. Do not set `session_start_time` in the YAML — it is always resolved at conversion time.

### Timestamp handling

Intan's `.rhd` format does **not** store an absolute recording timestamp. The pipeline uses a two-step fallback:

1. Try `get_metadata()` from NeuroConv (should be `None` for standard Intan files).
2. Parse DDMMYYYY from the immediate parent folder of `info.rhd`. Midnight of the session date is used as the approximation.

If the session folder name does not conform to DDMMYYYY, `master_neuroconv.py` raises a `ValueError` with a descriptive message. Rename the session folder or add an explicit `session_start_time` to the YAML as a fallback.

---

## 15. Contributors

| Name | Role |
|---|---|
| Jesus Ballesteros | Lead developer, pipeline architecture |
| Jonas Rose | Principal investigator |
| Aylin Apostel | Pipeline development |
| Lukas Hahn | Pipeline development |
| Sara Santos | Pipeline development |
| Juan Peschken | Pipeline development |
| Farina Lingstädt | Pipeline development |
| Winston Seah | Pipeline development |

---

## 16. External Tools & Citations

The NGL pipeline relies on the following open-source software and hardware. Please cite them in publications that use data processed with this toolbox.

**Hardware**

- **INTAN Technologies** — RHD2000 series multi-channel amplifier chips and RHX recording software. [intantech.com](https://intantech.com/)
- **Deuteron Technologies** — Neurolog miniature wireless logger. [deuterontech.com](https://www.deuterontech.com/)

**Spike sorting & QC**

- **Kilosort 4**: Pachitariu M, Sridhar S, Pennington J & Stringer C (2024). Spike sorting with Kilosort4. *Nature Methods*, 21, 914–921. https://doi.org/10.1038/s41592-024-02232-7
- **Bombcell**: Fabre JMJ, van Beest EH, Peters AJ, Carandini M & Harris KD (2023). Bombcell: automated curation and cell classification of spike-sorted electrophysiology data. *Zenodo*. https://doi.org/10.5281/zenodo.8172821
- **Phy**: Rossant C et al. — cortex-lab/phy. https://github.com/cortex-lab/phy

**LFP & signal processing**

- **FieldTrip**: Oostenveld R, Fries P, Maris E & Schoffelen JM (2011). FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. *Computational Intelligence and Neuroscience*, 2011, 156869. https://doi.org/10.1155/2011/156869

**NWB data standard**

- **Neurodata Without Borders (NWB)**: Rübel O, Tritt A, Ly R, Dichter BK et al. (2022). The Neurodata Without Borders ecosystem for neurophysiological data science. *eLife*, 11:e78362. https://doi.org/10.7554/eLife.78362
- **MatNWB**: NeurodataWithoutBorders/matnwb. https://github.com/NeurodataWithoutBorders/matnwb
- **NeuroConv**: Mayorquin H, Baker C et al. (2025). NeuroConv: Streamlining Neurophysiology Data Conversion to the NWB Standard. *Proceedings of the 24th Python in Science Conference*. https://doi.org/10.25080/cehj4257

**Analysis toolboxes**

- **CADopti** (cell assembly detection): Russo E & Durstewitz D (2017). Cell assemblies at multiple time scales with arbitrary lag constellations. *eLife*, 6:e19428. https://doi.org/10.7554/eLife.19428; Oettl LL et al. (2020).
- **npy-matlab**: Rossant C & colleagues — kwikteam/npy-matlab. https://github.com/kwikteam/npy-matlab
- **spikes**: Steinmetz N, Okun M & colleagues — cortex-lab/spikes. https://github.com/cortex-lab/spikes
- **prettify_matlab**: Julie Fabre — Julie-Fabre/prettify_matlab. https://github.com/Julie-Fabre/prettify_matlab

**Event coding**

Event-coding system based on the OTBR Toolbox:
`OTBR-Toolbox@ruhr-uni-bochum.de` · `gitlab.ruhr-uni-bochum.de/ikn/OTBR`
