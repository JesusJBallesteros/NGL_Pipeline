function result = process_gaze(csvPath, params)
%PROCESS_GAZE  Clean a DeepLabCut pose table and render the gaze animation.
%   result = PROCESS_GAZE(csvPath)
%   result = PROCESS_GAZE(csvPath, params)
%
%   Thin MATLAB wrapper around the GazEstim Python pipeline (pose_clean +
%   pose_render). Writes a JSON config from `params`, calls the Python
%   master script, and parses the RESULT line. All cleaning / rendering
%   logic lives in the Python side; this file exists only to bridge
%   MATLAB workspace params -> JSON -> subprocess.
%
%   NGL PIPELINE INTEGRATION:
%     Called by NGL06_videoAnalysis, which resolves masterScript and
%     background from <input.analysisCode>/ (project-local copies of the
%     toolbox's configfiles/master_gaze.py + configfiles/HexArena.png)
%     and toolboxes from <input.toolbox>/toolboxes/GazEstim/.
%
%   csvPath (required)  Path to the DLC .csv table (raw or filtered; the
%                       Python reader auto-detects header rows / index
%                       column / part columns).
%
%   params (struct). All fields optional (defaults shown) unless flagged
%   [REQUIRED FROM NGL06]:
%     .toolboxes      folder holding pose_clean.py / pose_render.py.
%                     Default: <toolbox>/toolboxes/GazEstim/ derived from
%                     this file's location (functions/video/process_gaze.m
%                     lives two levels below the toolbox root).
%     .masterScript   Python master script (configfiles/master_gaze.py).
%                     [REQUIRED FROM NGL06 — no toolbox-side fallback.
%                     The runtime copy MUST live under the user's
%                     analysisCode/ so per-project tweaks are visible.]
%     .background     background image (arena, bird removed).
%                     [REQUIRED FROM NGL06 — same rationale as masterScript.]
%     .pythonExe      python executable                        ('python')
%     .fps            input video fps                          (59.94)
%     .downsampleStep keep 1 frame per N (higher likelihood)   (2)
%     .targetFps      alternative to downsampleStep            (unset)
%     .startTime      analyse from this time [s]               (unset)
%     .endTime        analyse up to this time [s]              (unset)
%     .maxFrames      cap rendered output frames               (unset)
%     .maxSeconds     cap rendered output by duration [s]      (unset)
%     .pCut           likelihood limit                         (0.5)
%     .devFac         jump threshold = devFac*bodyLength       (0.6)
%     .smooth         temporal smoothing window [frames]       (5)
%     .wMed           rolling-median window [frames]           (9)
%     .boneTolFrac    bone-length tolerance fraction           (0.4)
%     .boneTolMad     bone-length tolerance (xMAD)             (5.0)
%     .orderMargin    head-behind-wing clamp margin (xbody)    (0.10)
%     .videoWidth     true DLC video width  [px]               (1250)
%     .videoHeight    true DLC video height [px]               (1160)
%     .output         output .mp4 (or .png if previewFrame)    (<csv>_gaze.mp4)
%     .gaze           draw gaze cones                          (true)
%     .monoFOV        monocular field per eye [deg]            (170)
%     .binoHalf       binocular half-angle [deg]               (15)
%     .coneMult       cone length = coneMult*birdLength        (2.5)
%     .eyeFwdFrac     eye base: frac of head->beak from head   (1/3)
%     .eyeLatFrac     eye lateral offset: frac of back-wing    (1/5)
%     .dpi/.crf/.preset  ffmpeg encoding options               (120/24/veryfast)
%     .previewFrame   render a single .png, no video           (unset)
%
%   result is a struct: .status .output .frames .out_fps .body_px .flagged .canon
%
%   Example (standalone, outside the pipeline):
%     p.masterScript = 'C:\Code\ephys-data-pipeline\configfiles\master_gaze.py';
%     p.background   = 'C:\Code\ephys-data-pipeline\configfiles\HexArena.png';
%     p.startTime = 10;  p.endTime = 40;
%     p.pCut = 0.6;      p.monoFOV = 160;
%     r = process_gaze('C:\data\bird01.csv', p);
%
% Last modified 26.06.2026 (Jesus) - integrated into ephys-data-pipeline:
%                                     drop hardcoded root, derive default
%                                     toolboxes path from mfilename;
%                                     masterScript + background become
%                                     required-via-params (no toolbox
%                                     configfiles fallback).

    if nargin < 2 || isempty(params), params = struct(); end
    assert(isfile(csvPath), 'process_gaze:csv', 'CSV not found: %s', csvPath);

    % Default toolboxes path: two levels up from this file
    % (functions/video/process_gaze.m -> <toolbox>/toolboxes/GazEstim/).
    thisFile        = mfilename('fullpath');
    thisDir         = fileparts(thisFile);                        % .../functions/video
    toolboxRoot     = fileparts(fileparts(thisDir));              % .../<toolbox>
    defaultTools    = fullfile(toolboxRoot, 'toolboxes', 'GazEstim');
    toolboxes       = getp(params, 'toolboxes', defaultTools);
    pythonExe       = getp(params, 'pythonExe', 'python');

    % masterScript + background are required-via-params. NGL06 sets them
    % from <input.analysisCode>/. If a standalone caller forgets them the
    % error is loud and points at the copy-from location.
    master     = getp(params, 'masterScript', '');
    background = getp(params, 'background',   '');
    assert(~isempty(master), 'process_gaze:noMaster', ...
        ['params.masterScript is required. Point it at your project''s ', ...
         '<analysisCode>\\master_gaze.py (copy the template once from ', ...
         '%s\\configfiles\\master_gaze.py).'], toolboxRoot);
    assert(~isempty(background), 'process_gaze:noBackground', ...
        ['params.background is required. Point it at your project''s ', ...
         '<analysisCode>\\HexArena.png (copy the template once from ', ...
         '%s\\configfiles\\HexArena.png).'], toolboxRoot);
    assert(isfile(master),     'process_gaze:master',     'master script not found: %s', master);
    assert(isfile(background), 'process_gaze:background', 'background not found: %s',   background);

    fpsIn = getp(params, 'fps', 59.94);
    step  = getp(params, 'downsampleStep', 2);

    % Default output naming: <csvName>_gaze.mp4 next to the CSV.
    [csvDir, csvName] = fileparts(csvPath);
    output = getp(params, 'output', fullfile(csvDir, [csvName '_gaze.mp4']));

    % Assemble config (schema matches master_gaze.py + pose_clean +
    % pose_render). MATLAB-side camelCase -> Python-side snake_case here.
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

    % Preview still (PNG) vs video (mp4). When previewFrame is set the
    % output extension is coerced to .png regardless of what the caller
    % passed.
    if isfield(params, 'previewFrame') && ~isempty(params.previewFrame)
        cfg.preview_frame = params.previewFrame;
        [od, on] = fileparts(output); output = fullfile(od, [on '.png']);
    end
    cfg.output = absify(output);

    % Optional windowing
    if isfield(params, 'targetFps') && ~isempty(params.targetFps), cfg.target_fps = params.targetFps; end
    if isfield(params, 'startTime') && ~isempty(params.startTime), cfg.start_time = params.startTime; end
    if isfield(params, 'endTime')   && ~isempty(params.endTime),   cfg.end_time   = params.endTime;   end
    if isfield(params, 'maxFrames') && ~isempty(params.maxFrames)
        cfg.max_frames = params.maxFrames;
    elseif isfield(params, 'maxSeconds') && ~isempty(params.maxSeconds)
        cfg.max_frames = round(params.maxSeconds * fpsIn / step);
    end

    % Cleaning thresholds (only forwarded when the caller set them).
    cl = struct();
    cl = addif(cl, 'p_cut',         params, 'pCut');
    cl = addif(cl, 'dev_fac',       params, 'devFac');
    cl = addif(cl, 'smooth',        params, 'smooth');
    cl = addif(cl, 'w_med',         params, 'wMed');
    cl = addif(cl, 'bone_tol_frac', params, 'boneTolFrac');
    cl = addif(cl, 'bone_tol_mad',  params, 'boneTolMad');
    cl = addif(cl, 'order_margin',  params, 'orderMargin');
    cfg.clean = cl;

    % Gaze cone tuning.
    gz = struct('enabled', logical(getp(params, 'gaze', true)));
    gz = addif(gz, 'mono_fov',     params, 'monoFOV');
    gz = addif(gz, 'bino_half',    params, 'binoHalf');
    gz = addif(gz, 'cone_mult',    params, 'coneMult');
    gz = addif(gz, 'eye_fwd_frac', params, 'eyeFwdFrac');
    gz = addif(gz, 'eye_lat_frac', params, 'eyeLatFrac');
    cfg.gaze = gz;

    % Write config to temp JSON + run python master.
    tmp = [tempname '.json'];
    fid = fopen(tmp, 'w');  assert(fid > 0, 'process_gaze:tempJson', 'cannot open temp config');
    fwrite(fid, jsonencode(cfg));  fclose(fid);

    cmd = sprintf('"%s" "%s" "%s"', pythonExe, master, tmp);
    fprintf('process_gaze: running\n  %s\n', cmd);
    [st, out] = system(cmd);
    fprintf('%s\n', out);
    if isfile(tmp), delete(tmp); end

    result = parseResult(out);
    if st ~= 0 || ~strcmp(getp2(result, 'status', 'error'), 'ok')
        error('process_gaze:failed', 'pipeline failed: %s', getp2(result, 'message', 'see output above'));
    end
    fprintf('process_gaze: wrote %s  (%d frames @ %.2f fps)\n', ...
            result.output, result.frames, result.out_fps);
end

% -----------------------------------------------------------------------
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
