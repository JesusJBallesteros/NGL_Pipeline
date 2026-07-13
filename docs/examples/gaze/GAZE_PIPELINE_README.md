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

## Key facts baked in
- DLC coords are in the **video** pixel space (default 1250×1160); the background PNG
  is stretched to that box so the skeleton overlays 1:1. Set `videoWidth/videoHeight`
  for other cameras.
- 60→30 fps downsample keeps the higher-likelihood frame of each pair.
- Gaze: eyes at ⅓ head→beak, ±⅕ back-wing laterally; monocular 170°, binocular ±15°,
  cones 2.5× body length.

## Dependencies
Python: `numpy`, `matplotlib`, `pillow`, `opencv-python` (only for make_plate), and
**ffmpeg** on PATH. Set `params.pythonExe` if `python` is not your interpreter.

## Regenerate the background plate (new photo)
```bash
python toolboxes/make_plate.py photo.png clean_plate.png x0 y0 x1 y1
```
`x0 y0 x1 y1` = the bird's bounding box in photo pixels.
