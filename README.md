# ephys-data-pipeline

MATLAB toolbox for processing INTAN RHX and Deuteron electrophysiology recordings at NGL.

Converts raw multi-channel recordings to Kilosort-ready `.bin` files, FieldTrip LFP structures, and trial-parsed event data. Runs Kilosort 4 spike sorting and Bombcell QC automatically. Supports multi-probe / multi-area recordings, NWB export, and head-direction motion-sensor extraction.

---

## Quick Start

1. Open `NGL_SetAndRunMe.m` and fill in your project paths and options.
2. Run `NGL_SetAndRunMe` — the pipeline dispatches automatically based on recording format.
3. Results land in `preprocessing/<subject>/<session>/` and `trialSorted/<subject>/<session>/`.

> For a full walkthrough see the **WIKI**.

---

## Pipeline Stages

| Script | Stage | Description | Depends on |
|---|---|---|---|
| `NGL00_Prep.m` | Preparation | Creates folder structure; validates `analysisCode/` | — |
| `NGL01_Main.m` | Preprocessing | Event extraction → `.bin` → Kilosort 4 → Bombcell QC | NGL00 |
| `NGL02_postPhy.m` | Spike analysis | Loads curated KS/Phy units; firing rates; population dynamics | NGL01 + **Phy curation** |
| `NGL02_LFP.m` | LFP analysis | Loads FieldTrip data; artifact rejection; time-frequency analysis | NGL01 (no curation needed) |
| `NGL03_acrossSession.m` | Aggregation | Cross-session and cross-subject pooling into cell arrays sized (subject × session). Gated by `opt.aggregateSessions` and `opt.aggregateSubjects` | NGL02_postPhy (NGL02_LFP later) |
| `NGL03_plotting.m` (TODO) | Visualisation | Group-level plots across sessions and conditions | NGL02_postPhy / NGL02_LFP / NGL03_acrossSession |

`NGL02_postPhy` and `NGL02_LFP` are siblings. The LFP path has no dependency on Phy curation, so it can be run as soon as `NGL01_Main` finishes — in parallel with manual curation if desired. All stages are launched via `NGL_SetAndRunMe.m`, which sets options and calls them in sequence.

---

## Supported Recording Platforms

| Platform | Format | Neural | Events | Motion |
|---|---|---|---|---|
| Intan RHX | `fileperch` (one file per channel) | `amp-*.dat` → `.bin` | Digital TTL via `readEvents` | ADXL335 on `AUX*.dat` channels |
| Deuteron DF1 | `DF1` | `NEUR*.DF1` → `.bin` | Deuteron event log (via `Event_File_Reader_9_0.exe`) | MPU-9250 (acc + gyro + mag) in stream 2 |

---

## Repository Layout

```
ephys-data-pipeline/
├── NGL_SetAndRunMe.m          Entry point — set options and run
├── NGL_machineConfig.m        Machine-specific paths (Python, toolboxes)
├── NGL00_Prep.m               Folder preparation
├── NGL01_Main.m               Main preprocessing loop
├── NGL02_postPhy.m            Post-sorting analysis
├── NGL03_plotting.m           Group-level plotting
├── set_default.m              Validates and merges opt with defaults
├── default_opt.m              Canonical defaults for every option field
├── Pipeline_LiveScript.m      Interactive live-script alternative to SetAndRunMe
│
├── functions/
│   ├── analysis/              MAT2FieldTrip, artifact detection, spectrograms, firing rate
│   ├── deuteron/              Deuteron2Kilosort, Deuteron2NWB, Deuteron2Fieldtrip,
│   │                          GetMotionSensors, Deuteron_PipelineWrapper, …
│   ├── ethology/              estimate_pecking, detect_jerkEvents, getSocialEvents
│   ├── events/                EventProcess, trialdefGen, sort2trials, …
│   ├── intan/                 INTAN_PipelineWrapper, Intan2Kilosort_*, intan2NWB_*, …
│   ├── pipeline/              findSessions, checkAnalysisCode, buildAreaMap, …
│   ├── plotting/              plot_single_fireRate, plot_multi_fireRate, densityScatterChart, …
│   ├── sorting/               master_kilosort4, Bombcell_Main, plot_KSresults, …
│   ├── utils/                 bandFilter, readyaml, downsampleVolt, …
│   └── video/                 vFLIP_NGL, processAndTrack_video, makeOrientationVideoFromMotionData
│
├── configfiles/               Templates to copy into your project's analysisCode\ folder
│   ├── bombcellConfig.m       Bombcell QC thresholds
│   ├── conditions_script.m    Condition grouping logic
│   ├── eventDefinitions.m     Event code → name mapping
│   ├── master_kilosort4.py    Kilosort 4 Python entry point
│   ├── master_neuroconv.py    NeuroConv NWB conversion script (INTAN)
│   ├── nwb_metadata_template.yaml  Project-level NWB metadata (YAML)
│   ├── parameters.py          Kilosort 4 parameters
│   └── postPhy_param.m        Post-Phy thresholds and flags
│
├── channelmaps/               Kilosort channel map .mat files for supported probes
├── toolboxes/                 Bundled dependencies (see below)
└── functions/_deprecated/     Retired code kept for reference
```

---

## Project `analysisCode\` Folder

Every project needs an `analysisCode\` folder containing:

| File | Required | Description |
|---|---|---|
| `eventDefinitions.m` | **Yes** | Maps event codes to names |
| `conditions_script.m` | **Yes** | Defines trial condition grouping |
| `bombcellConfig.m` | **Yes** | Bombcell QC thresholds |
| `postPhy_param.m` | **Yes** | Post-Phy analysis parameters |
| `chanMapXXX.mat` | **Yes** | Kilosort channel map (from `channelmaps/`) |
| `master_kilosort4.py` | **Yes** | Kilosort 4 entry point |
| `master_neuroconv.py` | **Yes** | NeuroConv conversion script |
| `nwb_metadata.yaml` | **Yes** | Project-level metadata |
| `parameters.py` | **Yes** | Kilosort 4 parameters |

Copy templates from `configfiles/` and fill in your project-specific values.

---

## Key Options (`opt` struct)

Set these in `NGL_SetAndRunMe.m` before running. All fields have safe defaults in `default_opt.m`.

| Option | Default | Description |
|---|---|---|
| `opt.numChannels` | `32` | Electrode count |
| `opt.bin` | `true` | Create Kilosort `.bin` file |
| `opt.FieldTrip` | `true` | Create FieldTrip LFP `.mat` file |
| `opt.doNWB` | `false` | Export to NWB format |
| `opt.RetrieveEvents` | `true` | Extract event log |
| `opt.alignto` | `{'itiOn'}` | Alignment event(s) for trial parsing |
| `opt.highpass` | `[]` | High-pass cutoff (Hz); `[]` = off |
| `opt.CAR` | `0` | Common-average re-referencing |
| `opt.kilosort` | `1` (= KS4) | Kilosort version |
| `opt.KSchanMapFile` | `''` | Channel map filename in `analysisCode\` |
| `opt.bombcell` | `true` | Run Bombcell QC after sorting |
| `opt.GetMotionSensors` | `false` | Extract accelerometer / IMU data |

> Full reference in [Wiki §6](wiki_NGL01_pipeline.md#6-all-options-reference)

---

## Multi-Area Mode

For recordings spanning more than one brain region, set `input.Areas` in `NGL_SetAndRunMe.m`:

```matlab
input.Areas = {'NCL', 'HP'};   % one label per kcoords group in chanMap
```

When `input.Areas` is set, `set_default` calls `buildAreaMap` to split channels by shank (`kcoords`), and Kilosort runs once per unique area. Results go to `preprocessing/<subj>/<session>/<Area>/`. `opt.KSchanMapFile` is required in multi-area mode.

---

## NWB Export

**INTAN:** set `opt.doNWB = true`. Requires a NeuroConv Python environment (path set in `NGL_machineConfig.m`). Calls `intan2NWB_neuroconv.m` → `master_neuroconv.py` via `pyrunfile`.

**Deuteron:** set `opt.doNWB = true`. Uses the bundled `matnwb` MATLAB toolbox (no Python needed). `Deuteron2Kilosort` saves a temporary `_raw.mat`; `Deuteron2NWB` reads it, builds an `NwbFile` with electrode geometry from the channel map, and writes `<session>.nwb`. YAML metadata is read by `readyaml` (`functions/utils/`).

Both paths read project-level metadata from `analysisCode/nwb_metadata.yaml`.

---

## Bundled Toolboxes

| Toolbox | Purpose |
|---|---|
| `Intan/` | Intan RHD2000 low-level file readers |
| `Deuteron/` | Deuteron DF1 extractor + motion-sensor constants |
| `fieldtrip_light/` | Lightweight FieldTrip preprocessing functions |
| `matnwb/` | MATLAB NWB read/write (used by Deuteron NWB path) |
| `npy-matlab/` | NumPy `.npy` reader (Kilosort / Phy output) |
| `bombcell/` | Automated spike-sorting QC |
| `spikes/` | Spike waveform utilities |
| `CADopti/` | CAD optimisation utilities |
| `prettify_matlab/` | Figure formatting helpers |
| `BDPAT_NGL/` | NGL-specific batch processing templates |
| `Viewer/` | Spike viewer |

---

## Branching Convention

| Branch | Purpose |
|---|---|
| `master` | Stable, tested releases |
| `maintenance/cleanup` | Active maintenance and bug fixes |
| `feature/<name>` | New feature development |

---

## Contributors

| Name | Affiliation | Contribution |
|---|---|---|
| Jesus Ballesteros | IKN, Ruhr-Universität Bochum | Pipeline architecture, maintenance |
| Jonas Rose | IKN, Ruhr-Universität Bochum | Lab PI; project direction |
| Aylin Apostel | IKN, Ruhr-Universität Bochum | Collaborator |
| Lukas Hahn | IKN, Ruhr-Universität Bochum | Collaborator |
| Sara Santos | IKN, Ruhr-Universität Bochum | Collaborator |
| Juan Peschken | IKN, Ruhr-Universität Bochum | Collaborator |
| Farina Lingstädt | IKN, Ruhr-Universität Bochum | Collaborator |
| Winston Seah | IKN, Ruhr-Universität Bochum | Collaborator |

---

## Citations

**Kilosort 4**
Pachitariu M, Sridhar S, Pennington J, Stringer C (2024). *Spike sorting with Kilosort4.* Nature Methods.

**Bombcell**
Bhagat J et al. (2024). *Bombcell: automated spike sorting quality control.* eLife.

**FieldTrip**
Oostenveld R, Fries P, Maris E, Schoffelen JM (2011). *FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data.* Computational Intelligence and Neuroscience.

**matNWB**
Teeuwen J et al. MatNWB: MATLAB interface for NWB files. Zenodo. `https://doi.org/10.5281/zenodo.6982050`

**readyaml**
Jongeneel MJ (2023). *readyaml — Read YAML files.* MATLAB Central File Exchange. `https://www.mathworks.com/matlabcentral/fileexchange/136369`

**OTBR Toolbox** (event-coding system)
`OTBR-Toolbox@ruhr-uni-bochum.de` · `gitlab.ruhr-uni-bochum.de/ikn/OTBR`
