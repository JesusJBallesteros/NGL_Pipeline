# Ephys Data Pipeline (NGL)

Scripts, functions, and toolboxes to process electrophysiology recordings at NGL. The pipeline converts raw INTAN/Deuteron recordings into spike-sorting and FieldTrip-compatible LFP outputs, supports automated Kilosort/Bombcell runs, and provides post-curation analysis utilities.

## Why this repo exists
This repository standardizes the workflow from **raw recordings** to **curated spikes and LFP analysis**, while enforcing a consistent on-disk data structure. It is designed to be the lab’s reproducible “one-stop” pipeline for typical NGL electrophysiology datasets.

## Quick start
1. **Set your MATLAB working directory** to this repository’s root (where `NGL01_Main.m` lives). Adding the repo to MATLAB’s permanent path is recommended.
2. **Edit `NGL_SetAndRunMe.m`** to configure your:
   - `datadrive`, `studyname`, and `toolbox` paths
   - subject/session selection (`subjects`, `dates`)
   - pipeline options (`opt.*`)
3. **Run the script**:
   - `NGL00_Prep` creates the folder structure for a new project.
   - `NGL01_Main` preprocesses data, exports Kilosort/FieldTrip files, and optionally runs Kilosort/Bombcell/Phy.
   - `NGL02_postPhy` performs post-curation processing (spikes, trial sorting, LFP workflows).

> **Note:** Raw data is expected to follow the IKN hard-disk structure. See:
> `gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure`

---

## Repository structure (high-level)
```
.
├── NGL_SetAndRunMe.m     # Main “driver” script you configure per project
├── NGL00_Prep.m          # Initializes project folder structure and README
├── NGL01_Main.m          # Main preprocessing pipeline (INTAN/Deuteron)
├── NGL02_postPhy.m       # Post-curation processing of spikes + LFP
├── default_opt.m         # Default pipeline options
├── set_default.m         # Input validation + dependency setup
├── functions/            # MATLAB helpers used by the pipeline
├── toolboxes/            # External toolboxes (FieldTrip, Kilosort utils, etc.)
├── channelmaps/          # Channel map files for probes
├── configfiles/          # Configs for specific processing steps
└── Instructions/         # Files to copy into a project’s analysisCode folder
```

---

## Pipeline overview

### 1) NGL00_Prep — Project initialization
- Creates a standard folder tree under `datadrive:\studyname\`.
- Writes a project `readme.txt` with basic metadata.

### 2) NGL01_Main — Preprocessing + spike sorting
- Detects session formats (INTAN vs. Deuteron vs. preprocessed FT).
- Converts raw data into:
  - `.bin` for Kilosort
  - FieldTrip-ready `.mat` for LFP
  - event/motion logs where applicable
- Runs Kilosort automatically (version 2 or 4).
- Optionally runs Bombcell and/or launches Phy for manual curation.

### 3) NGL02_postPhy — Post-curation and analysis
- Loads curated spike output.
- Builds trial-sorted neuron structures.
- Computes firing rates and optional tracking-related outputs.
- Loads LFP data and optionally performs trial parsing and time-frequency analysis.

---

## Configuration and defaults

### Main configuration (recommended)
Edit **`NGL_SetAndRunMe.m`** to set your project paths, subject/session selection, and pipeline options. This is intended to be the only file you routinely edit.

### Default options
`default_opt.m` defines a baseline `opt` structure (filters, Kilosort version, plotting parameters, etc.). You typically override these in `NGL_SetAndRunMe.m`.

### Dependency setup
`set_default.m` validates input settings, establishes default paths for analysis/raw/trial folders, and adds all required toolboxes to the MATLAB path.

---

## Important notes for newcomers
- **Manual curation is required.** Kilosort and Bombcell run automatically, but you still need to curate in Phy before post-processing spikes.
- **Stick to the expected data structure.** The pipeline assumes the IKN standard hard-disk layout and will break if data is placed elsewhere.
- **Keep toolboxes accessible.** `set_default.m` adds MATLAB paths to all required toolboxes; ensure those directories are present.

---

## What to learn next
If you’re new to this pipeline, these are great next steps:
- Understand **FieldTrip data structures** (LFP format and trial parsing).
- Learn the **Kilosort → Phy → postPhy** workflow end-to-end.
- Explore helper functions in `functions/` to see how events, trial definitions, and conversions are implemented.
- Review the `Instructions/` folder for project-specific configs (Kilosort configs, channel maps).

---

## Troubleshooting tips
- If NWB export fails after a previous run, you may need to restart MATLAB due to DLL conflicts (see `set_default.m`).
- If Phy launches but no data appears, verify that `.bin` and `params.py` are in the same session folder.

---

## License / attribution
This is an internal NGL electrophysiology pipeline. If you plan to reuse or publish parts of it, check with the lab maintainers.
