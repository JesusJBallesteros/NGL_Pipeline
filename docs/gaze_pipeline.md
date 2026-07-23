# Jackdaw gaze pipeline — wrapper

Run the full pose-cleaning + gaze-animation pipeline on any DeepLabCut table, from
MATLAB or from Python.

## Layout
```
<root>/
  functions/video/process_gaze.m   MATLAB entry point (prepares config, runs Python)
  configfiles/master_gaze.py       Python master: reads JSON config, runs the tools
  toolboxes/                       Python tools (the actual logic)
    pose_clean.py                  robust DLC read + physically-constrained cleaning
    pose_render.py                 skeleton + gaze cones over the arena background
    make_plate.py                  utility: remove the bird from a photo (cv2 inpaint)
  Jackdaw_clean_plate.png          default background (bird removed)
```
The three live in separate folders by request; `process_gaze.m` and `master_gaze.py`
derive each other's location from a shared `<root>`. Override with `params.root`,
`params.toolboxes`, `params.masterScript` if you relocate them.

## Run from MATLAB
```matlab
% simplest — all defaults
process_gaze('C:\data\bird01.csv');

% with options
p = struct();
p.startTime = 10;  p.endTime = 40;     % analyse 10–40 s of video
p.pCut = 0.6;                          % stricter likelihood gate
p.monoFOV = 160;  p.coneMult = 3;      % gaze-cone tuning
p.output = 'C:\out\bird01_gaze.mp4';
r = process_gaze('C:\data\bird01.csv', p);   % r.output, r.frames, r.flagged, ...
```
`help process_gaze` lists every parameter and its default.

## Run from Python (no MATLAB)
```bash
python configfiles/master_gaze.py myconfig.json
```
`myconfig.json` mirrors the params (snake_case): `csv`, `output`, `background`,
`video_w/h`, `fps_in`, `downsample_step`/`target_fps`, `start_time`/`end_time`,
`max_frames`, `clean:{p_cut,dev_fac,smooth,w_med,bone_tol_frac,bone_tol_mad,order_margin}`,
`gaze:{enabled,mono_fov,bino_half,cone_mult,eye_fwd_frac,eye_lat_frac}`, `dpi/crf/preset`,
and `preview_frame` (render one still PNG instead of a video).

## Body parts & roles
Labels are read from the csv, never hardcoded. See what a table contains:
```matlab
process_gaze('C:\data\bird01.csv', struct('listParts', true));
% -> available body parts: beak, tail, left_wing, right_wing, head, back
```
Use any subset:
```matlab
p.parts = {'beak','head','back'};      % wings/tail not tracked, or too unreliable
```
If your DLC project uses different names, map them onto the pipeline's anatomical
**roles** (`beak head back left_wing right_wing tail`):
```matlab
p.parts = {'bill','nape','back','tail'};
p.roles.beak = 'bill';   p.roles.head = 'nape';
```

What each role buys you:

| role(s) | used for |
|---|---|
| beak + head | gaze direction (**required** for cones) and eye placement |
| head + back | body-length scale, i.e. the jump threshold |
| back + wings | wing rigid bones, eye lateral offset, head-behind-wings rule |
| back + tail  | tail rigid bone, bird length that sets cone size |

Everything degrades gracefully — a missing role just switches off the checks that need
it (no wings → ordering rule skipped; no tail → bird length ≈ 4× body length; no
beak/head → cones skipped). Selected labels with **no** role are still cleaned and
drawn, with a generic colour and no bone links. Ask for a label that isn't in the file
and the error lists exactly what is available.

## Key facts baked in
- DLC coords are in the **video** pixel space (default 1250×1160); the background PNG
  is stretched to that box so the skeleton overlays 1:1. Set `videoWidth/videoHeight`
  for other cameras.
- 60→30 fps downsample keeps the higher-likelihood frame of each pair.
- Gaze: eyes at ⅓ head→beak, ±⅕ back-wing laterally; monocular 170°, binocular ±15°,
  cones 2.5× body length.

## Derived features (`.mat` + figure)
```matlab
p.features = true;         % -> <output>_estimated_features.mat  (+ a .png figure)
p.video    = false;        % optional: features only, no mp4 (much faster)
r = process_gaze(csv, p);  % r.features, r.figure
```
Loads in MATLAB as struct `features`; every measure is a column vector:

| field | meaning |
|---|---|
| `angle_deg` | head direction, **0° = video vertical (up), clockwise**, head→beak |
| `angvel_deg_s` | angular velocity of that vector (deg/s) |
| `centroid_x`, `centroid_y` | centre of gravity = unweighted mean of all selected labels |
| `speed_px_s` | speed of the centroid (`velocity_x/y_px_s` for the components) |
| `*_likelihood` | DLC likelihood: `min` for pairs (weakest link), `mean` for the centroid |
| `*_estimated` | 1 where the value depends on a joint the cleaner re-estimated |
| `*_sd`, `*_ci95` | approximate propagated uncertainty (caveats below) |
| `timestamp_ms` | **ms from recording start**; `frame_index` = original raw frame |
| `lever_px` | beak–head distance; short levers ⇒ noisier angles |

plus metadata (`convention`, `parts`, `beak_label`, `head_label`, `fps_in/out`,
`centroid_sigma_px`). `p.headDirection` is accepted as the old name for `p.features`.

The **figure** (`p.featuresFigure`, default on) has four panels: polar rose of head
direction, angular velocity vs direction, log-scaled occupancy map of the centre of
gravity, and speed over time.

Four honest caveats:
- **Sample spacing is not uniform** (≈16.7 / 33.4 / 50.0 ms at step 2) because the
  downsampler keeps the *higher-likelihood* frame of each pair. Derivatives already use
  the real timestamps; resample if you need a fixed grid.
- The angle is **circular**: it is unwrapped before differentiating. Head turns faster
  than ~half a turn between samples would still alias.
- `*_sd` is **geometric error propagation**, not a calibrated statistical CI, and it
  **understates** error where `*_estimated` is 1. Gate on `*_likelihood`/`*_estimated`.
- Speed is in **px/s**, and at 30 Hz its floor is set by residual tracking jitter —
  smooth or threshold before interpreting small values, and apply your own px→cm scale.

Needs `scipy` (for writing the .mat).

## Performance
Measured on the example table (15,930 frames after downsampling, 912×846 output):

| stage | cost |
|---|---|
| clean | **3.4 s total** (0.2 ms/frame) — negligible |
| render | **~27 ms/frame (~37 fps)** — dominates the runtime |
| memory | **~85 MB peak**, single CPU core (Agg is single-threaded) |

So a 900-frame/30 s clip takes **~28 s**, and the full table **~7 min**.

Rendering caches the static background once and redraws only the ~14 moving artists
(blitting), piping raw RGBA straight to ffmpeg — **9.4× faster** than redrawing the whole
figure per frame, and verified **pixel-identical** (max diff 0). Put `"legacy_render": true`
in the config to fall back to the slow matplotlib writer for comparison.

**GPU encoding is not worth enabling here.** ffmpeg encodes these frames at ~550 fps —
under 1% of the runtime — because the cost is matplotlib rasterisation, not encoding.
Note that `ffmpeg -encoders` lists `nvenc`/`qsv`/`amf` even when the hardware is absent;
only a real test encode tells you (`ffmpeg -f lavfi -i nullsrc=s=256x256:d=0.2 -c:v
h264_qsv -f null -`). The remaining per-frame cost is mostly the legend redraw, which
must sit above the cones.

## Dependencies / first run on a new machine
Python 3 with `numpy`, `matplotlib`, `pillow` (plus `opencv-python` only to regenerate
a plate), and **ffmpeg** on PATH for video output.

`process_gaze.m` does **not** trust bare `python`. It probes candidates in order —
`params.pythonExe` → the cached choice → MATLAB's `pyenv` interpreter → `python` /
`py -3` / `python3` → common install dirs (`…\Programs\Python\Python3*`, `anaconda3`,
`miniconda3`, `C:\Python3*`) — and picks the first that can **actually import** the
packages, caching it in MATLAB prefs so later runs are instant.

**`pipeline failed: No module named 'numpy'`** means python ran but was the *wrong*
interpreter (often the Microsoft Store stub, or a system python while your packages
live in Anaconda). Fixes, easiest first:
```matlab
p.autoInstall = true;   process_gaze(csv, p)   % pip install into the python it found
p.pythonExe = 'C:\Users\me\anaconda3\python.exe';   process_gaze(csv, p)
p.forcePythonSearch = true;                     % ignore a stale cached interpreter
```
…or once from a terminal: `python -m pip install numpy matplotlib pillow`.

**Canvas size.** The rendered canvas is `outWidth` px wide (default 912) with the height
derived from `videoWidth:videoHeight`. Both sides are forced **even**, because H.264 with
`yuv420p` rejects odd dimensions — e.g. a 640×512 video used to land on 912×**729** and
ffmpeg refused to start. The adjustment is at most 1 px. Set `params.outWidth` for a
bigger/smaller render.

**Render fails part-way with an ffmpeg message?** The renderer pipes frames to ffmpeg, so
if ffmpeg exits the write fails (on Windows as `Errno 22`, not a broken pipe). The error
now quotes ffmpeg's own message; the usual causes are the **output .mp4 being open in a
video player**, a non-writable folder, a full disk, or an output path without a video
extension. Close the player and re-run.

**Progress in MATLAB** streams live to the command window: `process_gaze.m` runs
`python -u` and does not capture stdout, reading the outcome from a result file instead.

Missing **ffmpeg** now fails fast with instructions instead of dying inside matplotlib:
`winget install Gyan.FFmpeg` (or `conda install -c conda-forge ffmpeg`), then restart
MATLAB. `params.previewFrame` renders a still PNG and needs no ffmpeg.

## Regenerate the background plate (new photo)
```bash
python toolboxes/make_plate.py photo.png clean_plate.png x0 y0 x1 y1
```
`x0 y0 x1 y1` = the bird's bounding box in photo pixels.
