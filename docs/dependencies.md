# Toolbox Dependencies

**ephys-data-pipeline** — complete dependency manifest  
Last updated: 2026-05-18

---

## MATLAB Official Toolboxes

The table below lists every MathWorks toolbox required to run the full NGL01 pipeline (INTAN + Deuteron paths). Toolboxes marked **Required** will produce errors if absent; those marked **Optional** degrade gracefully (the calling function checks availability or the code path is not triggered by default).

| Toolbox | Status | Functions used | Where |
|---|---|---|---|
| **Signal Processing Toolbox** | Required | `butter`, `filtfilt`, `decimate`, `bandpass`, `highpass`, `lowpass`, `fft`, `ifft`, `fir1`, `spectrogram`, `findpeaks`, `downsample` | Preprocessing (LFP, motion sensors), `bandFilter.m`, `GetMotionSensors.m` |
| **Statistics and Machine Learning Toolbox** | Required | `kmeans`, `ksdensity` | `kmeans` in analysis functions; `ksdensity` in `densityScatterChart.m` (optional density method) |
| **Sensor Fusion and Navigation Toolbox** | Required (Deuteron 9-DoF) | `ahrsfilter`, `quaternion`, `rotvec`, `eulerd`, `rotatepoint`, `slerp` | `Deuteron_estimateheading.m`, `Deuteron_PlotMotionSensors.m`, `makeOrientationVideoFromMotionData.m`; future `ProcessMotionSensors.m` (NGL02) |
| **Image Processing Toolbox** | Required (video path) | `rgb2gray`, `imadjust`, `imshow`, `imcrop` | `TreatVideo.m` |
| **Parallel Computing Toolbox** | Legacy / optional | `gather` (GPU array retrieval) | `rez2Phy.m` — legacy KS2/3 output converter only; not in active KS4 path |

**Notes:**
- `magcal` is **not** from the Sensor Fusion Toolbox. A standalone local implementation lives at `functions/ethology/magcal.m` and shadows the MathWorks version.
- No Control Systems, Curve Fitting, Symbolic Math, DSP System, or Optimization Toolboxes are required.
- No explicit MATLAB version guards exist in the codebase. The pipeline has been developed with R2023b–R2024b. `readyaml.m` replaces the R2023b-only `yamlread` for NWB metadata; R2021b+ should work for the rest.

---

## Bundled Third-Party Toolboxes

All toolboxes are vendored under `toolboxes/`. They must be on the MATLAB path before running the pipeline (`set_default.m` adds them automatically).

| Folder | Origin / Author | Purpose | Actively used in NGL01? |
|---|---|---|---|
| `BDPAT_NGL` | Lab-internal (NGL) | Spike-rate analysis utilities: `calcFireRate`, `plotRaster`, `plotPSTH`, `nanste`, `plotHeatmap`, `generateSpikeTrain` | Yes — analysis and plotting functions |
| `CADopti` | Russo & Durstewitz 2017 (*eLife* 6:e19428) | Cell assembly detection at optimal time scale | NGL02 level; not called in NGL01 |
| `Deuteron` | Deuteron Technologies (official SDK) | DF1/DT2 block-file parsing: `Deuteron_extractData`, `MotionSensorConstants`, `ScaleMotionSensorData`, `DataTypeEnum` | Yes — Deuteron recording path |
| `Intan` | Custom FieldTrip-compatible wrappers | Intan `.rhd`/`.rhs` header and data reading: `ft_read_data_INTAN`, `ft_read_header_INTAN` | Yes — LFP path |
| `Viewer` | MathWorks example (MPU-9250 orientation demo) | `HelperBox`, `HelperOrientationViewer`, `displayMessage` | `makeOrientationVideoFromMotionData.m` only |
| `bombcell` | Fabre et al. 2023, UCL ([Zenodo 10.5281/zenodo.8172821](https://doi.org/10.5281/zenodo.8172821)) | Automated spike quality metrics: `bc.qm.runAllQualityMetrics` | Yes — `Bombcell_Main.m` (NGL01 stage 05) |
| `fieldtrip_light` | Oostenveld et al., Donders Institute ([fieldtriptoolbox.org](https://www.fieldtriptoolbox.org)) | LFP preprocessing and TFR: `ft_preproc_*`, `ft_freqanalysis`, `ft_redefinetrial`, `ft_rejectartifact`, `ft_artifact_zvalue` | Yes — LFP/FieldTrip path |
| `matnwb` | NeurodataWithoutBorders ([github.com/NeurodataWithoutBorders/matnwb](https://github.com/NeurodataWithoutBorders/matnwb)) | NWB file creation and export: `NwbFile`, `types.core.*`, `types.hdmf_common.*`, `nwbExport`, `util.table2nwb` | Yes — `Deuteron2NWB.m` |
| `npy-matlab` | cortex-lab, UCL ([github.com/kwikteam/npy-matlab](https://github.com/kwikteam/npy-matlab)) | NumPy `.npy` file I/O: `writeNPY`, `readNPY` | Yes — `rez2Phy.m` (legacy KS2/3) and KS4 output reading |
| `prettify_matlab` | File Exchange [#154567](https://uk.mathworks.com/matlabcentral/fileexchange/154567) | Publication-quality plot styling | Not called in pipeline; interactive/manual use only |
| `spikes` | Harris & Carandini lab, UCL ([github.com/cortex-lab/spikes](https://github.com/cortex-lab/spikes)) | Post-KS visualization: `computeWFampsOverDepth`, `templatePositionsAmplitudes` | Yes — `plot_genstuff.m`, `plot_KSresults.m` |

### Standalone utility (not a full toolbox)

| File | Origin | Purpose | Used by |
|---|---|---|---|
| `functions/utils/readyaml.m` | Maarten J. Jongeneel, TU/e (Oct 2023) | Pure-MATLAB YAML parser. No R2023b+ dependency. Does **not** support `>` folded or `\|` literal block scalars — those values silently return empty string. | `Deuteron2NWB.m` |

---

## External Executables

| Executable | Purpose | Invoked from |
|---|---|---|
| **Python interpreter** (3.9 for KS4; 3.10+ for NeuroConv) | Runs `.py` scripts via MATLAB `pyrunfile`/`pyenv`. Two separate Conda environments are needed — KS4 and NeuroConv have conflicting dependencies and must not share an environment. | `master_kilosort4.m`, `intan2NWB_neuroconv.m` |
| **phy** >= 2.0 | Manual spike sorting curation GUI. Launched via `system('phy template-gui params.py')` from the KS output folder. | `NGL01_Main.m` (interactive stage, `opt.openPhy = true`) |
| **Conda / Miniconda** | Python environment manager. Required to create and activate the KS4 and NeuroConv environments. | Setup only (not called at runtime) |

---

## Python Packages

### Kilosort 4 environment (Python 3.9)

```bash
conda create --name kilosort python=3.9
conda activate kilosort
pip install kilosort[gui]
pip uninstall torch
conda install pytorch pytorch-cuda=11.7 -c pytorch -c nvidia
```

| Package | Purpose |
|---|---|
| `kilosort` >= 4.0 | Spike sorting engine. `run_kilosort`, `DEFAULT_SETTINGS`. |
| `torch` + CUDA 11.7 | GPU backend. CUDA-capable NVIDIA GPU strongly recommended; CPU is extremely slow. |
| `numpy` | Used in `parameters.py` for parameter definitions. |

MATLAB passes five arguments to `master_kilosort4.py`: `ksdir`, `binfile`, `nchan`, `chanmap`, `results_dir`.

### NeuroConv environment (INTAN → NWB path)

```bash
conda create --name neuroconv python=3.10
conda activate neuroconv
pip install neuroconv
```

| Package | Purpose |
|---|---|
| `neuroconv` >= 0.4 | INTAN recording conversion. `IntanRecordingInterface`, `load_dict_from_file`, `dict_deep_update`. Pulls in `spikeinterface`, `pynwb`, `hdmf`, `neo`, `probeinterface` as transitive dependencies. |
| `pynwb` (via neuroconv) | NWB Python API. |
| `pyyaml` (via neuroconv) | YAML metadata loading (`load_dict_from_file`). |
| `sys`, `datetime`, `zoneinfo`, `pathlib` | Standard library — no separate install needed. |

MATLAB passes three arguments to `master_neuroconv.py`: path to `info.rhd`, path for the output `.nwb` file, path to `nwb_metadata.yaml`.

---

## Dependency Tree (summary)

```
NGL01 pipeline
├── MATLAB base  (R2021b+)
│
├── Official toolboxes
│   ├── Signal Processing Toolbox          [always required]
│   ├── Statistics & ML Toolbox            [always required]
│   ├── Sensor Fusion & Navigation Toolbox [Deuteron 9-DoF path]
│   ├── Image Processing Toolbox           [video processing path]
│   └── Parallel Computing Toolbox         [legacy rez2Phy only]
│
├── Bundled toolboxes (toolboxes/)
│   ├── Deuteron       — Deuteron recording path
│   ├── Intan          — INTAN LFP path
│   ├── fieldtrip_light — LFP/TFR analysis
│   ├── matnwb         — Deuteron NWB export
│   ├── npy-matlab     — KS output I/O
│   ├── bombcell       — QC stage (NGL01 stage 05)
│   ├── spikes         — post-KS visualization
│   ├── BDPAT_NGL      — analysis and plotting
│   ├── Viewer         — orientation video
│   ├── CADopti        — NGL02 (not NGL01)
│   └── prettify_matlab — interactive only
│
├── Python: kilosort env  (Python 3.9)
│   ├── kilosort >= 4.0
│   ├── torch + CUDA 11.7
│   └── numpy
│
├── Python: neuroconv env  (Python 3.10+)
│   └── neuroconv >= 0.4  (+pynwb, spikeinterface, neo, ...)
│
└── CLI tools
    ├── phy >= 2.0      — manual spike sorting curation
    └── conda/miniconda — Python environment management
```
