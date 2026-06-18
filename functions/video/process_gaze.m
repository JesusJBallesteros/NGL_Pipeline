function result = process_gaze(csvPath, params)
%PROCESS_GAZE  Clean a DeepLabCut pose table and render a gaze animation.
%
%   result = PROCESS_GAZE(csvPath)
%   result = PROCESS_GAZE(csvPath, params)
%
% PURPOSE:
%   Thin MATLAB wrapper for the GazEstim Python pipeline. Writes a JSON
%   config from MATLAB inputs, then shells out to
%       <pythonExe> <masterScript> <config.json>
%   which loads pose_clean -> pose_render from toolboxes/GazEstim.
%   The cleaning/rendering logic lives entirely in Python.
%
% USAGE FROM THE NGL PIPELINE:
%   The expected caller is NGL06_VideoProcess (or a project script in
%   analysisCode). The orchestrator passes resolved paths via `params`:
%       params.pythonExe    = input.GAZEpythonExe;                       % from NGL_machineConfig
%       params.toolboxes    = fullfile(input.toolbox, 'toolboxes', 'GazEstim');
%       params.masterScript = fullfile(input.analysisCode, 'master_gaze.py');
%       params.background   = fullfile(input.analysisCode, 'HexArena.png');
%
%   If those fields are missing, this function falls back to locations
%   discovered from its own file path:
%       <toolbox>/functions/video/process_gaze.m       <- this file
%       <toolbox>/toolboxes/GazEstim/                  <- default toolboxes
%       <toolbox>/configfiles/master_gaze.py           <- example master
%       <toolbox>/configfiles/HexArena.png             <- example background
%   so it can also be run standalone from the toolbox root for testing.
%
% INPUTS:
%   csvPath (required)  Path to the DLC .csv table (raw or filtered; the
%                       reader auto-detects header rows / index column /
%                       part columns).
%
%   params (optional struct).  All fields optional; sensible defaults
%   shown below in parentheses.  See pose_clean.py / pose_render.py for
%   the underlying knobs.
%
%   Acquisition / windowing:
%     .fps            input video fps                          (59.94)
%     .downsampleStep keep 1 frame per N                       (2)  -> out fps = fps/N
%     .targetFps      alternative to downsampleStep            (unset)
%     .startTime      analyse from this time [s]               (unset = start)
%     .endTime        analyse up to this time [s]              (unset = end)
%     .maxFrames      cap rendered output frames               (unset = all)
%     .maxSeconds     cap rendered output by duration [s]      (unset = all)
%
%   Pose cleaning thresholds:
%     .pCut           DLC likelihood gate                      (0.5)
%     .devFac         jump threshold = devFac*bodyLength       (0.6)
%     .smooth         temporal smoothing window [frames]       (5)
%     .wMed           rolling-median window [frames]           (9)
%     .boneTolFrac    bone-length tolerance fraction           (0.4)
%     .boneTolMad     bone-length tolerance (xMAD)             (5.0)
%     .orderMargin    head-behind-wing clamp margin (xbody)    (0.10)
%
%   Render geometry:
%     .videoWidth     DLC video width  [px]                    (1250)
%     .videoHeight    DLC video height [px]                    (1160)
%     .background     background image (bird removed)         (configfiles/HexArena.png)
%     .output         output .mp4 (or .png if previewFrame)    (<csv>_gaze.mp4)
%
%   Gaze cones:
%     .gaze           draw gaze cones                          (true)
%     .monoFOV        monocular field per eye [deg]            (170)
%     .binoHalf       binocular half-angle [deg]               (15)
%     .coneMult       cone length = coneMult*birdLength        (2.5)
%     .eyeFwdFrac     eye base: frac of head->beak from head   (1/3)
%     .eyeLatFrac     eye lateral offset: frac of back-wing    (1/5)
%
%   Encoding:
%     .dpi / .crf / .preset  ffmpeg encoding options           (120 / 24 / veryfast)
%     .previewFrame   render a single still (frame index) -> PNG, no video (unset)
%
%   Component / executable overrides (used by NGL06_VideoProcess):
%     .pythonExe      python executable                       (input.GAZEpythonExe; falls back to 'python')
%     .root           ephys-data-pipeline root folder         (auto-detected)
%     .toolboxes      path to GazEstim python modules         (<root>/toolboxes/GazEstim)
%     .masterScript   path to master_gaze.py                  (<root>/configfiles/master_gaze.py)
%
% OUTPUT:
%   result struct: .status .output .frames .out_fps .body_px .flagged .canon
%
% Last modified 18.06.2026 (Jesus) - integrated into ephys-data-pipeline;
%                                     paths resolved via pipeline lookup
%                                     or explicit params overrides.

    if nargin < 2 || isempty(params), params = struct(); end
    assert(isfile(csvPath), 'process_gaze:csv', 'CSV not found: %s', csvPath);

    % Locate pipeline components relative to this file, then accept
    % overrides via params. The auto-detect lets process_gaze run
    % standalone from a freshly-cloned toolbox; NGL06_VideoProcess
    % always passes the resolved analysisCode-side paths.
    here     = fileparts(mfilename('fullpath'));        % <root>/functions/video
    rootAuto = fileparts(fileparts(here));              % <root>
    root      = getp(params, 'root',         rootAuto);
    toolboxes = getp(params, 'toolboxes',    fullfile(root, 'toolboxes', 'GazEstim'));
    master    = getp(params, 'masterScript', fullfile(root, 'configfiles', 'master_gaze.py'));
    pythonExe = getp(params, 'pythonExe',    'python');

    assert(isfolder(toolboxes), 'process_gaze:toolboxes', ...
        ['GazEstim Python toolboxes not found: %s\n', ...
         'Pass params.toolboxes pointing to ephys-data-pipeline/toolboxes/GazEstim.'], toolboxes);
    assert(isfile(master), 'process_gaze:master', ...
        ['master_gaze.py not found: %s\n', ...
         'Pass params.masterScript pointing to your analysisCode/master_gaze.py.'], master);

    fpsIn = getp(params, 'fps', 59.94);
    step  = getp(params, 'downsampleStep', 2);

    % defaults for output / background
    [csvDir, csvName] = fileparts(csvPath);
    output     = getp(params, 'output',     fullfile(csvDir, [csvName '_gaze.mp4']));
    background = getp(params, 'background', fullfile(root, 'configfiles', 'HexArena.png'));

    % assemble config (consumed by master_gaze.py / pose_clean / pose_render)
    cfg = struct();
    cfg.csv        = absify(csvPath);
    cfg.background = absify(background);
    cfg.toolboxes  = toolboxes;
    cfg.video_w    = getp(params, 'videoWidth',  1250);
    cfg.video_h    = getp(params, 'videoHeight', 1160);
    cfg.fps_in     = fpsIn;
    cfg.downsample_step = step;
    cfg.dpi    = getp(params, 'dpi', 120);
    cfg.crf    = getp(params, 'crf', 24);
    cfg.preset = getp(params, 'preset', 'veryfast');

    % preview still (PNG) vs video (mp4)
    if isfield(params, 'previewFrame')
        cfg.preview_frame = params.previewFrame;
        [od, on] = fileparts(output); output = fullfile(od, [on '.png']);
    end
    cfg.output = absify(output);

    % optional windowing
    if isfield(params, 'targetFps'), cfg.target_fps = params.targetFps; end
    if isfield(params, 'startTime'), cfg.start_time = params.startTime; end
    if isfield(params, 'endTime'),   cfg.end_time   = params.endTime;   end
    if isfield(params, 'maxFrames')
        cfg.max_frames = params.maxFrames;
    elseif isfield(params, 'maxSeconds')
        cfg.max_frames = round(params.maxSeconds * fpsIn / step);
    end

    % cleaning thresholds
    cl = struct();
    cl = addif(cl, 'p_cut',         params, 'pCut');
    cl = addif(cl, 'dev_fac',       params, 'devFac');
    cl = addif(cl, 'smooth',        params, 'smooth');
    cl = addif(cl, 'w_med',         params, 'wMed');
    cl = addif(cl, 'bone_tol_frac', params, 'boneTolFrac');
    cl = addif(cl, 'bone_tol_mad',  params, 'boneTolMad');
    cl = addif(cl, 'order_margin',  params, 'orderMargin');
    cfg.clean = cl;

    % gaze
    gz = struct('enabled', logical(getp(params, 'gaze', true)));
    gz = addif(gz, 'mono_fov',     params, 'monoFOV');
    gz = addif(gz, 'bino_half',    params, 'binoHalf');
    gz = addif(gz, 'cone_mult',    params, 'coneMult');
    gz = addif(gz, 'eye_fwd_frac', params, 'eyeFwdFrac');
    gz = addif(gz, 'eye_lat_frac', params, 'eyeLatFrac');
    cfg.gaze = gz;

    % write config + run python master
    tmp = [tempname '.json'];
    fid = fopen(tmp, 'w');  assert(fid > 0, 'cannot open temp config');
    fwrite(fid, jsonencode(cfg));  fclose(fid);

    cmd = sprintf('"%s" "%s" "%s"', pythonExe, master, tmp);
    fprintf('process_gaze: running\n  %s\n', cmd);
    [st, out] = system(cmd);
    fprintf('%s\n', out);
    if isfile(tmp), delete(tmp); end

    result = parseResult(out);
    if st ~= 0 || ~strcmp(getp2(result, 'status', 'error'), 'ok')
        error('process_gaze:failed', 'pipeline failed: %s', ...
            getp2(result, 'message', 'see output above'));
    end
    fprintf('process_gaze: wrote %s  (%d frames @ %.2f fps)\n', ...
            result.output, result.frames, result.out_fps);
end

% --- helpers ----------------------------------------------------------
function v = getp(p, f, d)
    if isstruct(p) && isfield(p, f) && ~isempty(p.(f)), v = p.(f); else, v = d; end
end

function s = addif(s, jsonField, p, matField)
    if isfield(p, matField) && ~isempty(p.(matField)), s.(jsonField) = p.(matField); end
end

function p = absify(p)
    p = char(p);
    if isempty(regexp(p, '^([A-Za-z]:|\\\\|/)', 'once')), p = fullfile(pwd, p); end
end

function v = getp2(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end

function r = parseResult(out)
    r = struct('status', 'error', 'message', 'no RESULT line in output');
    lines = regexp(out, '\r?\n', 'split');
    for i = 1:numel(lines)
        ln = strtrim(lines{i});
        if startsWith(ln, 'RESULT ')
            try, r = jsondecode(strtrim(ln(8:end))); catch, end
        end
    end
end
