# NGL06_videoAnalysis / GazEstim gaze pipeline

Per-session gaze-cone rendering from DeepLabCut pose tables. Wraps the GazEstim Python pipeline (`pose_clean` + `pose_render`, driven by `configfiles/master_gaze.py`) into the standard NGL subject × session fan-out.

## Files (inside the toolbox)

```
functions/video/process_gaze.m        MATLAB wrapper (JSON config -> Python subprocess)
configfiles/master_gaze.py            Python master (template; copy to analysisCode/)
configfiles/HexArena.png              default arena background (template; copy to analysisCode/)
toolboxes/GazEstim/pose_clean.py      DLC read + physically-constrained cleaning
toolboxes/GazEstim/pose_render.py     skeleton + gaze cones over the arena background
toolboxes/GazEstim/make_plate.py      one-off utility: bird-out photo -> clean plate
NGL06_videoAnalysis.m                 top-level stage (per subject x session fan-out)
docs/examples/gaze/                   bundled Example_data.csv + reference _rendered.mp4
```

## One-off project setup

Two files must be copied into your project's `analysisCode\` once per study. NGL06 hard-errors with the copy-from location if either is missing (rationale: per-project tweaks to the arena background or a patched master must live with the project, not the toolbox).

From MATLAB (or your OS shell):

```matlab
copyfile( fullfile(input.toolbox, 'configfiles', 'master_gaze.py'), ...
          fullfile(input.analysisCode, 'master_gaze.py') );
copyfile( fullfile(input.toolbox, 'configfiles', 'HexArena.png'), ...
          fullfile(input.analysisCode, 'HexArena.png') );
```

If you have a custom arena, regenerate the plate from a photo:

```bash
python <toolbox>\toolboxes\GazEstim\make_plate.py photo.png clean_plate.png x0 y0 x1 y1
```

`x0 y0 x1 y1` = bounding box of the bird in the photo (pixel coordinates). Save the resulting PNG as `<analysisCode>\HexArena.png` (or point `opt.gaze.background` at it).

## Machine configuration

`NGL_machineConfig.m` may set `GAZEpythonExe` to a specific Python interpreter (e.g. an env with numpy / matplotlib / pillow installed and `ffmpeg` on PATH). If unset, NGL06 falls back to `python` (system PATH). `opt.gaze.pythonExe` overrides per project.

At the start of every run NGL06 does a fail-fast smoke test that imports `numpy` / `matplotlib` / `PIL` and checks for `ffmpeg` on PATH; if anything is missing you get one clear error naming the missing dependency, not a cryptic subprocess failure per session.

## Input layout

```
<datadrive>:\<studyname>\data\behaviour\<subject>\<session>\*.csv
```

Exactly one DLC .csv per session. Zero → NGL06 writes `gazeSkipped_noCsv.txt` next to where the CSV should have been and moves on. Multiple → NGL06 writes `gazeSkipped_multiCsv.txt`, halts THIS session, and continues the batch (the safest behaviour: force disambiguation without killing the whole run).

## Output

Alongside each CSV, one of:

```
<csvName>_gaze.mp4        default: rendered video
<csvName>_gaze.png        when opt.gaze.previewFrame is set (single-frame QC)
```

Idempotent by default: `opt.gaze.overwrite = false` (the schema default) skips a session when the expected output file already exists. Set `opt.gaze.overwrite = true` to force re-render.

## `opt.gaze.*` reference

Master gate

- `opt.gaze.do` — turn NGL06 on.

Path / environment overrides (empty → sensible default)

- `opt.gaze.masterScript` — path to `master_gaze.py`. Empty → `<analysisCode>/master_gaze.py`.
- `opt.gaze.background` — path to background image. Empty → `<analysisCode>/HexArena.png`.
- `opt.gaze.pythonExe` — Python interpreter. Empty → `input.GAZEpythonExe` → `python`.

Windowing / decimation

- `opt.gaze.fps` — input video fps (default 59.94).
- `opt.gaze.downsampleStep` — keep 1 frame per N (default 2 → out_fps = fps/2).
- `opt.gaze.targetFps` — alternative to `downsampleStep`; mutually exclusive.
- `opt.gaze.startTime` / `opt.gaze.endTime` — seconds; both empty → whole clip.
- `opt.gaze.maxFrames` / `opt.gaze.maxSeconds` — cap on output length; mutually exclusive.
- `opt.gaze.previewFrame` — render one PNG at this frame index instead of the mp4.

Video geometry

- `opt.gaze.videoWidth` / `opt.gaze.videoHeight` — source DLC video pixel dimensions (default 1250 × 1160).

Cleaning thresholds (physically-constrained pose repair)

- `opt.gaze.pCut` — DLC likelihood gate (default 0.5).
- `opt.gaze.devFac` — jump threshold = `devFac * bodyLength` (default 0.6).
- `opt.gaze.smooth` — temporal smoothing window (frames, default 5).
- `opt.gaze.wMed` — rolling-median window (frames, default 9).
- `opt.gaze.boneTolFrac` — bone-length tolerance fraction (default 0.4).
- `opt.gaze.boneTolMad` — bone-length tolerance (× MAD, default 5.0).
- `opt.gaze.orderMargin` — head-behind-wing clamp margin (× body, default 0.10).

Gaze cones

- `opt.gaze.drawCones` — render cones (default true).
- `opt.gaze.monoFOV` — monocular field per eye (degrees, default 170).
- `opt.gaze.binoHalf` — binocular half-angle (degrees, default 15).
- `opt.gaze.coneMult` — cone length = `coneMult * birdLength` (default 2.5).
- `opt.gaze.eyeFwdFrac` — eye base fraction of head→beak from head (default 1/3).
- `opt.gaze.eyeLatFrac` — eye lateral offset fraction of back→wing (default 1/5).

Encoding

- `opt.gaze.dpi` — render DPI (default 120).
- `opt.gaze.crf` — ffmpeg CRF (default 24; lower = better, larger).
- `opt.gaze.preset` — ffmpeg preset (default `'veryfast'`).

Behaviour

- `opt.gaze.overwrite` — skip if output exists (default false → skip); true to re-render.

## Cross-field guards (in `optPostChecks.m`)

- `downsampleStep` and `targetFps` are mutually exclusive.
- `maxFrames` and `maxSeconds` are mutually exclusive.
- `startTime < endTime` when both are set.
- `previewFrame` set → info warning that `crf` / `preset` are ignored.

## Smoke test

`docs/examples/gaze/` bundles `Example_data.csv` and `Example_data_rendered.mp4`. To smoke-test the wrapper (outside the NGL fan-out):

```matlab
p                 = struct();
p.masterScript    = fullfile(input.analysisCode, 'master_gaze.py');
p.background      = fullfile(input.analysisCode, 'HexArena.png');
p.pythonExe       = input.GAZEpythonExe;
r = process_gaze( fullfile(input.toolbox, 'docs', 'examples', 'gaze', 'Example_data.csv'), p );
```

`r.output` should be a new mp4 alongside the CSV; compare visually to `Example_data_rendered.mp4`.

## Failure modes

- `NGL06:noMaster` / `NGL06:noBackground` — copy the template as shown above.
- `NGL06:noToolbox` — `input.toolbox` doesn't point at the ephys-data-pipeline root; check `NGL_machineConfig.m`.
- `NGL06:pythonEnv` — one of numpy/matplotlib/pillow is missing, or ffmpeg isn't on PATH. Install via `pip install numpy matplotlib pillow` in the Python env pointed at by `GAZEpythonExe`.
- `NGL06:multipleCsvs` — the session has more than one .csv; keep exactly one and re-run (or curate manually).
- Per-session failure — check the `<csvName>_gazeFailed.txt` sidebar next to the CSV; it holds the MATLAB error message and full stack.

## See also

- `functions/video/process_gaze.m` — the wrapper, with the full parameter surface documented in its header.
- `configfiles/master_gaze.py` — the Python master; also runnable standalone with `python master_gaze.py myconfig.json`.
- `docs/examples/gaze/GAZE_PIPELINE_README.md` — original standalone README (kept for reference).
