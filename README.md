# ephys-data-pipeline

MATLAB toolbox for processing INTAN RHX and Deuteron electrophysiology recordings at NGL.

Converts raw multi-channel recordings to Kilosort-ready `.bin` files, FieldTrip LFP structures, and trial-parsed event data. Runs Kilosort 4 spike sorting and Bombcell QC automatically.

📖 **Full operator guide → [wiki_NGL01_pipeline.md](wiki_NGL01_pipeline.md)**

---

## Quick Start

1. **Clone** this repository to a stable local path (e.g. `C:\Code\ephys-data-pipeline`).
2. **Copy** `NGL_SetAndRunMe.m` and `configfiles/eventDefinitions.m` into your project's `analysisCode/` folder.
3. **Create** `NGL_machineConfig.m` in `analysisCode/` with paths to your Python environments and toolbox root (see [Wiki §4](wiki_NGL01_pipeline.md#4-machine-configuration-ngl_machineconfig)).
4. **Edit** `NGL_SetAndRunMe.m` for your study name, subjects, dates, and recording options.
5. **Open** `NGL_SetAndRunMe.m` in MATLAB and run it (F5 or section by section with F9).

---

## Requirements

**MATLAB:** R2021b or later

**Python environments (separate Conda envs):**
- Kilosort 4 env: `pip install kilosort`
- Phy 2 env: `pip install phy`

**Bundled toolboxes** (present under `toolboxes/`): FieldTrip light, Intan RHD reader, Deuteron reader, Bombcell, npy-matlab, and others. All paths are added automatically at runtime by `set_default` — do not add them manually to MATLAB's permanent path.

---

## Pipeline Flow

```
NGL_SetAndRunMe  →  NGL00_Prep  →  NGL01_Main  →  [Phy curation]  →  NGL02_postPhy
```

| Script | What it does |
|---|---|
| `NGL00_Prep` | Creates the IKN standard folder structure on disk (run once per project) |
| `NGL01_Main` | Detects format, extracts events, creates `.bin` + LFP `.mat`, runs KS4 + Bombcell |
| `NGL02_postPhy` | Post-curation: builds spike matrices, trial-parsed LFP, spike-field variables |

---

## Repository Layout

```
ephys-data-pipeline/
│
├── NGL01_Main.m              ← top-level pipeline (do not edit)
├── NGL02_postPhy.m           ← post-curation processing
├── NGL_SetAndRunMe.m         ← user config TEMPLATE (copy to analysisCode/)
├── NGL00_Prep.m              ← IKN folder structure creator
├── default_opt.m             ← single source of truth for all option defaults
├── set_default.m             ← option validation, path builder, dependency loader
│
├── functions/                ← all pipeline functions
│   ├── findSessions.m            session directory discovery
│   ├── prepforsession.m          per-session setup and path assignment
│   ├── chckV.m                   recording format detection
│   ├── findSetting.m             INTAN header reader
│   ├── reduceChanMap.m           reduced channel map generator
│   ├── EventProcess.m            event extraction dispatcher + caching
│   ├── INTAN_ExtractEvents.m     INTAN digital-input pin reader
│   ├── INTAN_PipelineWrapper.m   INTAN stage orchestrator
│   ├── Intan2Kilosort_wrapper.m  .bin file creator (dispatch)
│   ├── Intan2Kilosort_fileperch.m  fileperch → .bin
│   ├── Intan2Kilosort_filepertype.m  filepertype → .bin
│   ├── intan2MAT_wrapper.m       INTAN → pseudo-FieldTrip LFP
│   ├── MAT2FieldTrip.m           continuous + trial-parsed FT .mat creator
│   ├── trialdefGen.m             trial boundary builder from EventRecord
│   ├── events2align.m            event-name to decimal-value resolver
│   ├── Deuteron_PipelineWrapper.m  Deuteron stage orchestrator
│   ├── Deuteron_ExtractEvents.m    Deuteron event log reader (EXE or text log)
│   ├── Deuteron2Kilosort.m         DF1 → Kilosort .bin file creator
│   ├── Deuteron2Fieldtrip.m        DF1 → pseudo-FieldTrip LFP struct
│   ├── Deuteron_GetMetaData.m      Deuteron hardware constants lookup
│   ├── Deuteron_extractData.m      Low-level DF1 block reader (neural/motion/audio)
│   ├── GetMotionSensors.m          Deuteron/INTAN motion sensor dispatcher
│   ├── Deuteron_estimateheading.m  AHRS heading + dead-reckoning from MPU-9250
│   ├── Deuteron_PlotMotionSensors.m  Motion sensor / orientation visualiser
│   ├── master_kilosort4.m          MATLAB → Python KS4 caller
│   ├── Bombcell_Main.m           Bombcell QC wrapper
│   ├── downsampleVolt.m          integer-factor voltage downsampler
│   ├── binvec2dec.m              LSB-first binary vector to decimal
│   └── ...
│
├── configfiles/              ← TEMPLATES — copy to analysisCode/ before use
│   ├── eventDefinitions.m    ← event code definitions (copy + customise per project)
│   ├── NGL_machineConfig.m   ← machine-specific paths (copy + fill in, do not commit)
│   └── ...
│
├── toolboxes/                ← bundled third-party toolboxes (read-only)
│
├── wiki_NGL01_pipeline.md    ← full operator guide
└── README.md                 ← this file
```

---

## Folder Structure on Disk (IKN Standard)

```
<datadrive>:\<studyname>\
├── analysisCode\       NGL_SetAndRunMe, eventDefinitions, chanMaps, ...
└── data\
    ├── raw\            original INTAN .dat / Deuteron .DT2/.DF1 files
    ├── preprocessing\  .bin, KS4 output, EventRecord.mat
    ├── spikeSorted\    post-Phy spike variables
    ├── trialSorted\    trialdef.mat, events.mat, FieldTrip .mat files
    ├── behaviour\
    └── analysis\
```

Reference: `gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure`

---

## Key Options (`opt` fields)

Set these in `NGL_SetAndRunMe.m`. All unset fields receive safe defaults from `default_opt.m`.

| Field | Default | Notes |
|---|---|---|
| `numChannels` | `32` | Expected electrode count |
| `KSchanMapFile` | `''` | `''` = linear array; or `'chanMapXXX.mat'` |
| `alignto` | `{'itiOn'}` | Alignment event(s) for LFP trial parsing |
| `FieldTrip` | `false` | Produce FieldTrip LFP `.mat` file |
| `bombcell` | `true` | Run Bombcell QC after sorting |
| `phy` | `false` | Open Phy after KS4 (blocks MATLAB) |
| `lowpass` | `9000` | Spike-band low-pass Hz (`[]` = off) |
| `lowpassFT` | `200` | LFP low-pass Hz (Butterworth 4th order) |
| `highpass` | `[]` | High-pass Hz (`[]` = off; set to 300 for Deuteron Kilosort input) |
| `CAR` | `0` | Common-average re-referencing |
| `uselog` | `false` | **Deuteron only** — `true` = parse `logevents.txt` instead of using EXE |
| `GetMotionSensors` | `false` | Extract head-direction / accelerometer data |

→ Full reference in [Wiki §6](wiki_NGL01_pipeline.md#6-all-options-reference-opt-fields)

---

## Branching Convention

| Branch | Purpose |
|---|---|
| `master` | Stable, tested code |
| `maintenance/cleanup` | Active maintenance and bug fixes |
| `feature/<name>` | New features |

---

## Citation

Event-coding system based on the OTBR Toolbox:
`OTBR-Toolbox@ruhr-uni-bochum.de` · `gitlab.ruhr-uni-bochum.de/ikn/OTBR`
