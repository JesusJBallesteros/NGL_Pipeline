# NGL06_VideoProcess — design plan

Status: **draft / not implemented**. This document captures the scope, layout,
inputs/outputs and open decisions for the new video-processing block.
Implementation is queued as a separate task.

---

## 1  Purpose

NGL06 is an **independent stream** that runs from raw video data, in
parallel with NGL01 (raw-ephys preprocessing). Its outputs feed back
into the same per-session analysis tree the spike / LFP pipeline uses,
so that downstream condition tables and analysis scripts can read pose,
gaze and tracking products alongside the neural data.

It bundles two existing branches that currently live as scattered
options:

- **Pose / gaze estimation** (B.16 in `NGL_SetAndRunMe`).  
  Cleans a DeepLabCut `.csv` pose table and renders a gaze-cone video.
  Wrapper: `functions/video/process_gaze.m` →
  `configfiles/master_gaze.py` → `toolboxes/GazEstim/{pose_clean,pose_render}.py`.
- **Blob / social tracking** (B.12.a `opt.offlineTrack`, `opt.useTrack`).  
  Detects bird / object blobs in raw video and indexes them against
  trial events for spike-vs-social-event analysis. Lives in
  `functions/video/{TreatVideo,processAndTrack_video}.m` and friends.

NGL06 orchestrates both: discover the video files, dispatch the
configured branch(es), and write outputs into a predictable per-session
folder.

---

## 2  Inputs

### 2.1  Filesystem inputs

NGL06 looks primarily into `<datadrive>:\<studyName>\data\behaviour\`:

```
data/behaviour/
    <subject>/
        <YYYYMMDD>/
            session.mp4                                <- raw video (any name)
            session_DLC.csv                            <- DeepLabCut output (any name; we glob *.csv)
            optional: session.json / metadata.txt      <- per-session annotations
```

Discovery contract:

- For every `(subject, date)` in `input.subjects` × `input.dates`,
  walk `data/behaviour/<subject>/<date>/` and return:
  - all `*.mp4` / `*.avi` / `*.mkv` (raw video sources)
  - all `*.csv` whose name contains `DLC` or `DeepCut` (pose tables)
- Allow per-paradigm overrides via `opt.video.discoveryGlob` etc.

### 2.2  Configuration inputs (from `NGL_SetAndRunMe`)

| group              | switches |
|--------------------|----------|
| Master gates       | `opt.gaze.do`, `opt.offlineTrack`, `opt.useTrack` |
| Gaze knobs (B.16)  | `opt.gaze.*` (fps, pCut, monoFOV, ... see optSchema) |
| Blob / social      | existing `opt.offlineTrack`, `opt.useTrack` + any future `opt.blob.*` |
| Project paradigm   | `opt.proj_socialLearning` etc. (already present) |
| Python executable  | `input.GAZEpythonExe` (set by `NGL_machineConfig`) |
| Script overrides   | `opt.gaze.masterScript` (per-project copy in `analysisCode/`) |

NGL06 itself adds **no new opt fields**; it consumes the existing ones.

---

## 3  Per-paradigm branches

Two dispatch axes:

| axis          | values                                          |
|---------------|-------------------------------------------------|
| **Detection** | `pose` (DeepLabCut CSV) ∣ `blob` (in-video)     |
| **Use**       | `gaze` ∣ `social` ∣ both                        |

Routing matrix (rough first draft — refine on first run-through):

| Project flag                       | Detection | Use     | Script(s) called                              |
|------------------------------------|-----------|---------|-----------------------------------------------|
| `opt.gaze.do = true`               | pose      | gaze    | `process_gaze(csv, gazeParams)`               |
| `opt.offlineTrack = true`          | blob      | social  | `TreatVideo` / `processAndTrack_video`        |
| `opt.useTrack = true`              | n/a       | indexing| `sort2trials_blob` (already exists)           |
| `opt.proj_socialLearning = true`   | blob      | social  | + ASL-specific helpers                        |

A single (subject, date) can trigger multiple branches; NGL06 runs them
in a fixed order (detection → cleaning → indexing) so per-trial outputs
are consistent.

---

## 4  Outputs

Mirroring the spike-side layout (`data/analysis/<subject>/<sess>/`):

```
data/analysis/<subject>/<sess>/video/
    pose/
        <video>_clean.csv                  <- cleaned DLC table
        <video>_gaze.mp4                   <- rendered gaze video
        <video>_gaze.json                  <- pose_clean info struct (frames, body_px, flagged, canon, ...)
    blob/
        <video>_blob.mat                   <- TreatVideo output
        <video>_blob_track.mat             <- processAndTrack_video output
    indexing/
        <video>_socialEvents.mat           <- sort2trials_blob output
    NGL06_manifest.txt                     <- human-readable summary of what ran
```

Per-run summary appended to `data/analysis/<subject>/<sess>/preprocInfo.mat`
under field `videoMeta` so the existing preprocInfo loaders see it.

---

## 5  Pipeline sequence inside NGL06

```
00.  NGL00_Prep + set_default + findSessions  (standard scaffolding)
01.  For each (subject, sess):
       01.a  Discover raw videos + DLC csv(s)        -> videoMeta struct
       01.b  Mark a stale-output check               -> skip if up-to-date
       01.c  Dispatch:
                if opt.gaze.do:        process_gaze(csv, gazeParams)
                if opt.offlineTrack:   TreatVideo / processAndTrack_video
                if opt.useTrack:       sort2trials_blob (after blob track)
       01.d  Save NGL06_manifest.txt + update preprocInfo.mat
02.  Optional cross-session aggregation (LATER — same pattern as NGL03)
```

The stale-output check (01.b) lets re-runs skip already-processed
sessions, mirroring the firepools cache idea.

---

## 6  Integration touchpoints

### 6.1  Path & python wiring (DONE)
- `toolboxes/GazEstim/` already added to MATLAB path in `set_default` Section 7.
- `input.GAZEpythonExe` already populated in `set_default` Section 5 with a `'python'` fallback.
- `process_gaze.m` already accepts `params.masterScript / .toolboxes / .pythonExe` overrides; NGL06 will pass:
  ```matlab
  params.pythonExe    = input.GAZEpythonExe;
  params.toolboxes    = fullfile(input.toolbox, 'toolboxes', 'GazEstim');
  params.masterScript = fullfile(input.analysisCode, 'master_gaze.py');
  params.background   = resolveBg(opt.gaze.background, input);   % helper TBD
  ```

### 6.2  Existing blob / social functions (TO AUDIT)
- `functions/video/TreatVideo.m`
- `functions/video/processAndTrack_video.m`
- `functions/events/sort2trials_blob.m`

NGL06 will not rewrite these; it will only call them through a thin
adapter so the file naming / output folder convention from §4 holds.

### 6.3  Downstream wiring
- `conditions_script(_)` (project-specific) can read
  `data/analysis/<subj>/<sess>/video/indexing/*.mat` to attach
  per-trial social-event tags onto the `condition` struct, the same
  way `attachMergeSession` does for the INTAN merger output.

---

## 7  Open questions (block implementation start)

1. **Video discovery convention.** Are video filenames standardised across
   projects (e.g. always `<date>_<cam>.mp4`)? If yes, we can make
   discovery deterministic; if no, NGL06 needs `opt.video.discoveryGlob`
   per project.
2. **DLC pre-step.** Does NGL06 *also* run DeepLabCut, or is the CSV
   assumed to be already on disk? (My read: CSV is pre-computed; NGL06
   only consumes it. Confirm.)
3. **Multi-video per session.** Some sessions have multiple cameras /
   files. Process each independently with the same params, or merge?
4. **Output cache key.** Per-CSV (so re-running for a new gaze param
   regenerates) or per-(session, paradigm)?
5. **Cross-session aggregation.** Do we want an NGL03-style
   `aggregated_video.mat` collecting pose / blob outputs across the
   study, or is per-session enough for now?
6. **Schema additions.** Beyond the current `opt.gaze.*`, do we need
   `opt.blob.*` (currently the blob/social branch is configured via
   loose flags). If yes, this is the right moment to define them.
7. **DeepLabCut model registry.** When per-project DLC models live
   somewhere persistent (e.g. `analysisCode/DLC_models/`), the
   discovery needs to know which model maps to which paradigm.

---

## 8  Tasks to file once implementation starts

- [ ] Build `functions/video/discoverVideoFiles.m` — returns a struct per (subject, session).
- [ ] Build `functions/video/resolveGazeBackground.m` — picks analysisCode > configfiles fallback.
- [ ] Write `NGL06_VideoProcess.m` top-level orchestrator.
- [ ] Audit `TreatVideo.m` / `processAndTrack_video.m` interfaces, write adapters.
- [ ] Add `data/analysis/<subj>/<sess>/video/` to the standard output tree.
- [ ] Update the diagrams in `docs/diagrams/` to add NGL06 between NGL01 and NGL02.
- [ ] Schema fields for the blob/social half (`opt.blob.*`) once we have a working `TreatVideo` interface.

---

## 9  Cross-references

- `functions/video/process_gaze.m` — the MATLAB wrapper.
- `configfiles/master_gaze.py` — the Python orchestrator copied to user `analysisCode/`.
- `toolboxes/GazEstim/{pose_clean,pose_render,make_plate}.py` — the Python toolbox.
- B.12.a `opt.offlineTrack` / `opt.useTrack` — pre-existing blob/social switches.
- B.16 `opt.gaze.*` — current gaze knobs (added in this integration pass).
