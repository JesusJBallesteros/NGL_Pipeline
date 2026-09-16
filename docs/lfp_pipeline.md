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

### (g) Event-centered power contrasts — `computeTFRcontrast` (Phase 1)

Gated on `opt.lfp.session.contrast`. For every alignment × area × entry in
`opt.lfp.contrast.pairs`:

1. `parseTrialContrast` turns `'correct vs error'` into two trial masks.
   Each side is a condition field (non-zero = in), `~field` for its
   complement, or `all`. A struct with `.A` / `.B` masks covers contrasts a
   condition field cannot express. A length mismatch against the TFR's trials
   raises — artifact rejection after `condition.mat` was written is the usual
   cause, and silently comparing the wrong trials is worse than stopping.
   Overlapping sides warn: the independent-samples test assumes they don't.
2. `computeTFRcontrast` normalises each side (maps) and tests the contrast on
   **raw** power (statistics), via the Phase 0 helpers. Averaging is trials
   first, then channels.
3. `plotTFRcontrast` draws condition A, condition B (shared colour scale) and
   their difference with the cluster outlined, one row per band.

Output: `<SavFileName>_LFP_TFRcontrast_<area>_<align>_<A>-vs-<B>.{mat,png}`.
The TFR itself comes from `computeTrialparsedTFR`, whose cache means turning
this on alongside (a) costs one TFR, not two.

Bands keep separate colour scales here, unlike the quick-look plot: power
falls steeply with frequency, and one scale across a 4 Hz and a 60 Hz band
flattens the faster one. The question here is each band's shape, not which
band is larger.

### (h) Current source density — `computeCSD` (Phase 2)

Gated on `opt.lfp.session.csd`. Per alignment × shank, from the trial-parsed
data:

    CSD = -sigma * d²(phi)/dz²

**Sign convention: a sink is negative, a source positive.** A sink is current
entering cells — the extracellular negativity under an active input — and the
figure draws it blue.

`lfpChannelGeometry` reads the channel map (`xcoords`, `ycoords` in µm,
`kcoords` per shank) and returns one entry per shank with contacts in depth
order. Channels are matched to map rows **by position** — the same assumption
Kilosort makes reading the `.bin` — and the function refuses when the counts
differ rather than pairing the first N and hoping (a reduced map from
`reduceChanMap` is the usual cause; point `'chanMap'` at the reduced one).

Things that matter:

- **Uniform spacing is required and enforced.** The second difference weights
  unequal gaps wrongly; an irregular array needs an inverse method (iCSD),
  which is not implemented.
- **Smoothing across contacts** (`opt.lfp.csd.smoothPasses`, default 1) is the
  main knob. Differentiating amplifies noise, so some smoothing is standard —
  too little gives a speckled map, too much merges neighbouring sinks.
- **Vaknin's extension** (default on) duplicates the end contacts so CSD is
  defined at every depth. The two edge depths then rest on the assumption that
  the potential is flat beyond the probe, not on a measurement; the figure
  labels them.
- The CSD is computed on the **trial average**. `'keeptrials'` returns
  per-trial CSD, which is far noisier and mainly useful as input to a
  statistic.
- Colour limits use wide percentiles (0.2–99.8) here, unlike the TFR panels:
  an evoked CSD occupies tens of milliseconds of the epoch, and the default
  robust limits would scale to the baseline and saturate the response.

Verified by round trip on the real ATLAS 2-shank map (16 contacts, 50 µm): a
known CSD integrated twice into a potential, then recovered — sink returned
within half a contact spacing of where it was planted, at the right latency,
with the other shank flat.

Output: `<SavFileName>_LFP_CSD_<area>_<align>_shank<N>.{mat,png}`.

### (f) Spectrolaminar mapping — `spectrolaminarFLIP`

Thin wrapper around `vFLIP_NGL`. Identifies superficial / deep channels via low- vs high-freq power crossover along the shank (Mendoza-Halliday et al. 2024). Output: `<SavFileName>_LFP_FLIP.mat`.

## Shared foundation (Phase 0, Sep 2026)

Five helpers every power analysis is built on. They exist so that a power map,
a contrast, a comodulogram and a CSD from the same session are normalised,
tested, drawn and named the same way — and so each of those decisions travels
with the file that resulted from it.

### `normalizeTFR` — power relative to what

`[freq, info] = normalizeTFR(freq, opt)`. Methods: `db` (default), `relchange`,
`percent`, `z`, `absolute`, `none`. Window from `opt.lfp.norm.baseline`,
falling back to `opt.lfp.plot.baseline`.

Unlike `ft_freqbaseline` it keeps trials (statistics need them), adds the
z-score against baseline *variability*, and refuses input that would silently
produce `Inf`: a zero/negative baseline yields NaN with one warning naming the
channels — normalising an already-normalised TFR is the usual cause. Normalise
first, average second: with dB the two orders differ.

`info` (method, window, samples, single-trial, units) goes into provenance.

### `lfpClusterStats` — one corrected test

`[stat, info] = lfpClusterStats(A, B, opt)`, wrapping `ft_freqstatistics`
cluster permutation (Maris & Oostenveld 2007). Designs:

- `trials` (default) — two conditions' trials, independent samples.
- `paired` — same units in two conditions (sessions/subjects at group level).
- `baseline` — a condition against its own baseline, dependent samples. Run it
  on **raw** power: normalising first makes the comparison circular.

What it returns is a cluster-level claim — where a cluster is significant the
data differ *somewhere in it*; the bins inside are not individually
significant and the cluster's edges are not a confidence interval. `info`
carries counts, p-values and `pFloor = 1/numrand`, the smallest p obtainable.

Set `opt.lfp.stats.latency` to the post-event window when that is the question;
leaving the baseline in weakens a contrast and makes a `baseline` design partly
circular.

FieldTrip's *"Not all replications are used for the computation of the
statistic"* is expected on wavelet TFRs — the edge cone is NaN, so those bins
hold fewer trials than the design lists. FieldTrip drops them per bin and the
test stays valid; a `latency` inside the cone-free window silences it.

### `lfpStyle` / `plotTFRpanel` — one panel, reused

`lfpStyle(opt)` resolves colour map, limits, fonts and size once;
`lfpStyle(opt, 'diverging')` gives a blue-white-red map with limits symmetric
about zero, for signed quantities (contrasts, CSD) where zero must be legible.

`plotTFRpanel(ax, t, f, M, st, ...)` draws one time × frequency panel: robust
colour limits (2nd–98th percentile, so one artifact bin cannot flatten the
map), the event line, and significance as an **outline** rather than a blank —
a non-significant trend stays visible instead of being hidden by the test. NaN
bins (a wavelet's edge cone) are transparent, so they cannot be read as a
strong effect.

### `lfpResultName` / `saveLFPresult` / `saveLFPfigure` — where things land

`<SavFileName>_LFP_<kind>[_<area>][_<align>][_<tags>]` with `.mat` and `.png`
differing only in extension. `saveLFPresult` attaches `provenance` via
`buildLFPProvenance`, including the `normalize` and `statistics` info structs,
and saves `-v7.3` so a TFR with trials is not truncated.

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

### CSD (Phase 2)
- `opt.lfp.session.csd` — gate for NGL07 (h)
- `opt.lfp.csd.conductivity` (0.3 S/m), `.smoothPasses` (1), `.vaknin` (true)
- `opt.lfp.csd.trials` — condition field restricting the trials averaged; `''` = all
- `opt.lfp.csd.plot` (true)

### Contrasts (Phase 1)
- `opt.lfp.session.contrast` — gate for NGL07 (g)
- `opt.lfp.contrast.pairs` — e.g. `{'correct vs error', 'stim2 vs ~stim2'}`
- `opt.lfp.contrast.areas` — `{}` = every area in `chanArea`
- `opt.lfp.contrast.plot` — write the three-panel figure (default true)

### Normalisation (Phase 0)
- `opt.lfp.norm.method` (`'db'` | `'relchange'` | `'percent'` | `'z'` | `'absolute'` | `'none'`; default `'db'`)
- `opt.lfp.norm.baseline` (default `[-0.5 0]`)
- `opt.lfp.norm.singleTrial` (default true)

### Statistics (Phase 0)
- `opt.lfp.stats.design` (`'trials'` | `'paired'` | `'baseline'`; default `'trials'`)
- `opt.lfp.stats.numrand` (default 1000), `.alpha` (0.05), `.clusteralpha` (0.05)
- `opt.lfp.stats.clusterstatistic` (`'maxsum'`), `.tail` (0 = two-sided)
- `opt.lfp.stats.avgoverchan` (default true), `.minnbchan` (0)
- `opt.lfp.stats.latency` / `.frequency` (`'all'` or `[lo hi]`)

### Figures (Phase 0)
- `opt.lfp.plot.divergingColormap` (`''` = built-in blue-white-red)
- `opt.lfp.plot.figSize` (`[1100 700]`), `.fontSize` (10)

### Behaviour regression (stub)
- `opt.lfp.behReg.covariates` (default `{'x','y','head_dir'}`)
- `opt.lfp.behReg.numrand` (default 1000)
- `opt.lfp.behReg.alpha` (default 0.05)

## Neural frequency tagging — NGL08_NFT (project stage, Sep 2026)

Project-specific (Q1 convention): analyses in
`functions/analysis/projects/NFT/`, behind `opt.nft.do`, driven by
`NGL08_NFT.m`. Input is the **continuous** FT file, not the trial-parsed one:
tagging is read from long stretches of periodic stimulation, not from
event-locked epochs.

### The idea

A tagged response sits at a frequency the experiment chose, so it is judged
against its own neighbourhood — the neighbouring bins carry the same noise at
almost the same frequency, measured at the same time. Two numbers per
frequency: **SNR** (amplitude ÷ mean neighbour, 1 when nothing is there) and
**z** ((amplitude − mean) ÷ SD of neighbours).

Amplitude alone cannot separate driving from ongoing rhythm, so **ITPC**
across epochs is computed from the same transform: a driven response keeps its
phase every epoch and the unit vectors add; ongoing activity at the same
frequency drifts and cancels.

### The chain

| Step | Function | What matters |
|---|---|---|
| Epoch | `nftEpochs` | Epochs hold a **whole number of stimulation cycles**, so the tag lands exactly on an FFT bin. A fractional epoch smears the peak *and* raises the baseline it is measured against — both errors push the response down. |
| Spectrum | `computeTaggingSpectrum` | `mtmfft`, amplitude averaged over epochs, complex spectrum kept. Amplitude is averaged (not the complex values), so a response drifting in phase survives here and is judged by ITPC instead. |
| Response | `taggingResponse` | SNR/z per bin; harmonics of the base, minus any colliding with mains or with an excluded frequency; summed baseline-corrected amplitude over significant harmonics. |
| Phase | `computeITPC` | ITPC, Rayleigh z and p. Biased upward at small n — the function warns under 10 epochs. |
| Figure | `plotTaggingSpectrum` | z spectrum with harmonics marked, response per harmonic, ITPC per harmonic. |

### Event convention — to confirm against real data (Sep 2026)

Stated by the owner, **not yet verified against a recording**; treat as the
working assumption until a real session is analysed:

- Every stimulus carries a **`stimOn`** event followed by its own code. Codes
  are arbitrary but continuous per stimulus (A–E ≈ 8001–8006), defined in
  `analysisCode/eventDefinitions.m` (copied from `configfiles/`; reserved
  codes 0–15 are hardware-locked, project codes start at 16).
- Stimuli are paired into fixed **"syllables"** (AB, CD, FE): the second
  member always follows the first, while what precedes the first varies.
- That asymmetry is the experiment: the **transition probability** creates
  structure at half the component rate, and finding neural energy at that
  frequency is the point of the analysis.

What this means for the analysis, once confirmed:

- Block windows can be derived from the first and last `stimOn` of each block
  rather than supplied by hand, which is what `opt.nft.blocks` currently
  wants. A helper turning event codes into blocks is the obvious next piece.
- Two frequencies per block, not one: the **component rate** (2.6 Hz) and the
  **syllable rate** (1.3 Hz = component/2, the pairing). Run the block twice,
  with `.base` set to each. The component peak says the stimuli drove a
  response at all; the syllable peak says the transition structure was picked
  up. They are different claims.
- The **random-sequence blocks are the control**: with no fixed pairing there
  should be a component peak and *no* syllable peak. A syllable-rate peak in a
  random block means something other than learning produced it — a periodicity
  in the stimulus set, or a block-window boundary artifact.
- Beware the harmonic collision: the syllable rate's 2nd harmonic *is* the
  component rate. `taggingResponse` already drops harmonics that land on an
  excluded frequency — pass the component rate in `'exclude'` when testing the
  syllable rate, or the component response will be read as evidence of
  structure.

### Block windows are yours to supply

`opt.nft.blocks` is a struct array of `.name`, `.base` (Hz) and `.window`
`[t0 t1]` in seconds. Mapping blocks to windows depends on the stimulation
log and the study's event codes, so the stage does not guess: a session whose
blocks are unknown is skipped with a warning rather than analysed against
invented boundaries.

For the three-block design in use (1.3 Hz stream, 2.6 Hz stream, 2.6 Hz
components in pairs with a gap), note that the pair structure of block 3
repeats at base/3 — run that block a second time with `.base` set to the pair
rate to test it, as the demo does.

### Caveats worth keeping in mind

- `zThreshold` (default 1.64, one-sided p = 0.05) is a **screen, not a
  corrected test**: across 6 harmonics a pure-noise channel has roughly a 1-in-4
  chance of one harmonic crossing it. Judge a channel by its base frequency and
  by ITPC, not by "at least one significant harmonic".
- `sumAmp` sums only significant harmonics, so it cannot go below zero. Compare
  channels by it only when both have a significant harmonic; otherwise compare z.
- Overlapping epochs (`opt.nft.overlap` > 0) are not independent: both the
  Rayleigh test and the ITPC bias correction become optimistic.

### `opt.nft.*`
- `nft.do`, `nft.blocks`, `nft.areas` (`{}` = all), `nft.plot`
- `nft.epochSeconds` (20), `nft.overlap` (0), `nft.taper` (`'hanning'`), `nft.fmax` (40)
- `nft.neighbours` (12), `nft.gap` (1), `nft.zThreshold` (1.64)
- `nft.maxHarmonic` (8), `nft.lineFreq` (50; 0 = ignore)

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
