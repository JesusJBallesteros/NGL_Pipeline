# NGL LFP pipeline (NGL02_LFP + NGL07_LFPanalysis)

Two LFP stages sit alongside the spike-side chain. They serve different roles.

- **NGL02_LFP** — per-session quick-look. Runs after NGL01. Verifies the FT file loads, optionally rejects artifacts, computes a lightweight continuous or trial-parsed TFR. Purpose: catch bad sessions early.
- **NGL07_LFPanalysis** — research-grade session-level LFP. Runs after NGL02_postPhy (needs spike sort) and, optionally, NGL06_videoAnalysis (behaviour sidecar). Does the trial-parsed TFR per condition, oscillation / burst detection, phase + envelope, spike-field coupling, LFP × behaviour regression, and spectrolaminar (vFLIP) mapping.

## Data flow

```
NGL01_Main
   |     writes <SavFileName>_FTcont.mat and <SavFileName>_<align>.mat
   |     with FT_data.chanArea populated per channel (Pass 2, 26.06.2026).
   v
NGL02_postPhy          NGL02_LFP (quick-look)
   |                        |
   |  spike.mat             |  continuous / trial-parsed TFR .mat
   v                        v
NGL06_videoAnalysis   (research-grade LFP not here)
   |
   |  DLC-driven gaze/position/HD (mp4/png today; .mat sidecar planned)
   v
NGL07_LFPanalysis      <-- Pass 3 stage, this doc
```

## Multi-area design (Q3 (b), Pass 2)

One FT file per session with `FT_data.chanArea{i}` per channel — NOT one FT file per area. Analyses that want per-area outputs restrict channels via `opt.lfp.tfrAreaFilter` and `ft_selectdata`. This keeps cross-area work (coherence, PPC across areas) feasible without materialising a separate file per area.

Legacy (pre-Pass-2) FT files without `chanArea` are auto-backfilled at load time by `ensureChanArea`.

## NGL07 analyses

Each analysis is opt-gated (`opt.lfp.session.<name>`). NGL07 skips silently when its gate is off.

### (a) Trial-parsed TFR — `computeTrialparsedTFR`

Per opt.alignto: `ft_freqanalysis` with `cfg.keeptrials='yes'`. Per-band via `opt.freqInterest`; supports `wavelet` / `mtmconvol` / `superlet`. Fingerprinted cache (`tfrCacheKey` / `loadTFRcache` / `saveTFRcache`) invalidated by mtime of the source FT file. Output: `<SavFileName>_LFP_TFR_<align>.mat`.

### (b) Oscillation / burst detection — `detectBursts`

First-pass threshold detector. Per band: bandpass → Hilbert envelope → log → z-score per channel → threshold at `opt.lfp.burst.threshMult` σ → require ≥ `opt.lfp.burst.minDurationMs` ms duration. Output: `<SavFileName>_LFP_bursts.mat` with a MATLAB `table` of episodes.

Limitation: does not separate 1/f aperiodic activity from true rhythmic bursts. Upgrade path: install **eBOSC** (extended Better OSCillation, github/BOSCbase/eBOSC) under `toolboxes/eBOSC/` and rewrite `detectBursts` as a wrapper over eBOSC's robust 1/f fit; keep the same signature/output shape.

### (c) Continuous phase + envelope — `hilbertBandpass`

Per band: `ft_preproc_bandpassfilter` (Butterworth-4 twopass) → `hilbert()` → `angle()` / `abs()`. Output: `<SavFileName>_LFP_phase.mat` with `bp.phase{b}` and `bp.envelope{b}` as `[nCh × nSamples]`. `dtype` defaults to `single` to keep files bounded (30 min × 1 kHz × 32 chan × 3 bands × single ≈ 700 MB).

### (d) Spike-field coupling — `spikeFieldCoupling`

For every cluster with ≥ `opt.lfp.spikeField.minSpikes`:
- **PPC**: `ft_spiketriggeredspectrum` → `ft_spiketriggeredspectrum_stat` with `cfg.method = opt.lfp.spikeField.ppcMethod` (default `'ppc2'`, robust to spike-rate dependencies per Vinck 2011).
- **Coherence**: `ft_spike_convert2fieldtrip` (spikes as LFP-rate 0/1 series) → `ft_appenddata` → `ft_freqanalysis` (`mtmfft`, fourier) → `ft_connectivityanalysis` with `cfg.method='coh'`, `cfg.complex='absimag'` (FT-recommended for volume-conduction bleed suppression).

Output: `<SavFileName>_LFP_spikeField.mat` with `sfc.ppc` `[nClust × nFreq]` and `sfc.coherence` `[nClust × nChanLFP × nFreq]`.

Cross-area falls out for free: index `sfc.coherence` on `sfc.chanArea` dimension to get e.g. NCL-cluster ↔ STR-LFP coherence.

### (e) LFP × behaviour regression — `lfpBehaviorRegression` (STUB)

**Not yet usable.** Blocked on NGL06 emitting a `<csvName>_gaze.mat` sidecar with per-frame `t / x / y / head_dir / gaze_az`. Add the sidecar-save hook on the Python side (`configfiles/master_gaze.py`); once available, this function will do:
- resample covariate onto TFR trial-time grid
- `ft_freqstatistics` with `cfg.statistic = 'ft_statfun_depsamplesregrT'`, cluster-permutation
- output `LFP_behReg_<covariate>.mat` per covariate

### (f) Spectrolaminar mapping — `spectrolaminarFLIP`

Thin wrapper around `vFLIP_NGL`. Identifies superficial / deep channels via low- vs high-freq power crossover along the shank (Mendoza-Halliday et al. 2024). Output: `<SavFileName>_LFP_FLIP.mat`.

## `opt.lfp.*` reference (Pass 2 + Pass 3)

### NGL02_LFP quick-look (Pass 1)
Flat legacy names still supported:
- `opt.doLFPthing`, `opt.trialparsed`, `opt.spectrogram`
- `opt.artifdet`, `opt.artZvalue`, `opt.rejValue`
- `opt.freqInterest`, `opt.TFRmethod`, `opt.superletOrder`, `opt.width`, `opt.combine`, `opt.timeResol`, `opt.blocks`, `opt.chbych`, `opt.trialbytrial`

### Nested (Pass 2 / Pass 3)
- `opt.lfp.tfrAreaFilter` — restrict computation to a subset of channels by area (char or cellstr).
- `opt.lfp.tfrCacheDir` — override for the TFR cache folder.
- `opt.lfp.bands` — cell of `{name, [flo fhi]}` per band. Default: theta / beta / gamma.

### NGL07 gates
- `opt.lfp.session.do` — master gate.
- `opt.lfp.session.tfr` / `.bursts` / `.phase` / `.spikeField` / `.behReg` / `.flip` — per-analysis gates.

### Burst detector
- `opt.lfp.burst.threshMult` (default 3)
- `opt.lfp.burst.minDurationMs` (default 100)

### Spike-field
- `opt.lfp.spikeField.minSpikes` (default 50)
- `opt.lfp.spikeField.ppcMethod` (`'ppc0'` | `'ppc1'` | `'ppc2'`; default `'ppc2'`)
- `opt.lfp.spikeField.timwin` (default 0.5 s)
- `opt.lfp.spikeField.foi` (default 2:2:100)

### Hilbert
- `opt.lfp.hilbert.storePhase` / `.storeEnv` (default true)
- `opt.lfp.hilbert.dtype` (`'single'` | `'double'`, default `'single'`)

### vFLIP
- `opt.lfp.flip.laminaraxis` (default 0:0.05:1.55)
- `opt.lfp.flip.freqaxis` (default 1:150)
- `opt.lfp.flip.setfreqbool` (0 = vFLIP, 1 = default fixed FLIP bands)

### Behaviour regression (stub)
- `opt.lfp.behReg.covariates` (default `{'x','y','head_dir'}`)
- `opt.lfp.behReg.numrand` (default 1000)
- `opt.lfp.behReg.alpha` (default 0.05)

## Provenance (Pass 2)

Every LFP output file carries a `provenance` struct via `buildLFPProvenance`:
- `source_FTfile`, `source_mtime`
- `opt_snapshot`
- `matlabVer`, `ftVer`
- `savedAt`, `host`, `callingFunction`
- plus per-analysis extras (area, alignment, band, cluster id, …).

Use it as first-line diagnosis when a downstream plot looks off — the exact opt struct that produced the file is right there.

## What ships with Pass 3 vs. what's still deferred

Shipped:
- (a) TFR ✓
- (b) burst detection (threshold-based first pass) ✓
- (c) phase + envelope ✓
- (d) spike-field coupling ✓
- (f) spectrolaminar FLIP ✓
- Full NGL07 driver ✓

Deferred:
- Aggregation across sessions (`NGL03` LFP entries) — Q4 answer, TODO.
- Cross-subject `NGL04_TFR` plotter — Q4 answer, TODO.
- eBOSC upgrade of `detectBursts` — install-time decision.
- (e) LFP × behaviour regression — waiting on NGL06 gaze.mat sidecar.

## Pointers

- `functions/analysis/ensureChanArea.m` — one source of truth for per-channel area tagging.
- `functions/analysis/compareByBlock.m` — generic per-trial contrast primitive; use to build project-specific analyses instead of embedding early/late/block-vs-block loops.
- `functions/analysis/projects/socialLearning/computeTrialparsedTFR_ASL.m` — example of a project-specific wrapper that could be rebuilt on top of `computeTrialparsedTFR` + `compareByBlock`.
- `functions/_deprecated/LFP_Fieldtrip.m` / `continous_MTspectrogram.m` — for reference only; do not call.
