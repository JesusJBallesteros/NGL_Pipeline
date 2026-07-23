function result = process_gaze(csvPath, params)
%PROCESS_GAZE  Clean a DeepLabCut pose table and render the gaze animation.
%   result = PROCESS_GAZE(csvPath)
%   result = PROCESS_GAZE(csvPath, params)
%
%   MATLAB wrapper: it prepares a JSON config from your inputs and runs the
%   Python master pipeline (../../configfiles/master_gaze.py), which uses the tools
%   in ../../toolboxes. The cleaning/rendering lives entirely in Python.
%
%   csvPath (required)  Path to the DLC .csv table (raw or filtered; the reader
%                       auto-detects header rows / index column / part columns).
%
%   params (optional struct).  All fields optional; sensible defaults shown:
%     .parts          cellstr of labels to use    (default: every part in the csv)
%     .roles          struct role->label for non-standard names (identity by name)
%                       roles: beak head back left_wing right_wing tail
%                       e.g. p.roles.beak = 'bill';  p.roles.head = 'nape';
%     .bodyPx         body scale [px] if no rigid bone is available     (30)
%     .listParts      print the labels available in the csv, then return
%     .features       save derived features as <output>_estimated_features.mat
%                       (true | path).  Loads as struct `features` with: angle_deg,
%                       angvel_deg_s, centroid_x/y, speed_px_s, each with
%                       *_likelihood / *_estimated / *_sd / *_ci95, plus
%                       timestamp_ms (from recording start) and frame_index.
%                       Angle: 0 deg = video vertical (up), clockwise, head->beak.
%                       (.headDirection is accepted as the old name for this.)
%     .featuresFigure summary figure next to the .mat            (true | path)
%     .video          render the mp4                            (true)
%                       set false to compute headDirection only (much faster)
%     .fps            input video fps                          (59.94)
%     .downsampleStep keep 1 frame per N (higher likelihood)   (2)  -> final output fps/N
%     .targetFps      alternative to downsampleStep            (unset)
%     .startTime      analyse from this time [s]               (unset = start)
%     .endTime        analyse up to this time [s]              (unset = end)
%     .maxFrames      cap rendered output frames               (unset = all)
%     .maxSeconds     cap rendered output by duration [s]      (unset = all)
%     .pCut           likelihood limit                         (0.5)
%     .devFac         jump threshold = devFac*bodyLength       (0.6)
%     .smooth         temporal smoothing window [frames]       (5)
%     .wMed           rolling-median window [frames]           (9)
%     .boneTolFrac    bone-length tolerance fraction           (0.4)
%     .boneTolMad     bone-length tolerance (xMAD)             (5.0)
%     .orderMargin    head-behind-wing clamp margin (xbody)    (0.10)
%     .videoWidth     true DLC video width  [px]               (1250)
%     .videoHeight    true DLC video height [px]               (1160)
%     .background     background image                         (<root>/<filename>.png)
%     .output         output .mp4 (or .png if previewFrame)    (<csv>_gaze.mp4)
%     .gaze           draw gaze cones                          (true)
%     .monoFOV        monocular field per eye [deg]            (170)
%     .binoHalf       binocular half-angle [deg]               (15)
%     .coneMult       cone length = coneMult*birdLength        (2.5)
%     .eyeFwdFrac     eye base: frac of head->beak from head   (1/3)
%     .eyeLatFrac     eye lateral offset: frac of back-wing    (1/5)
%     .dpi/.crf/.preset  DLC's encoding options                (120/24/veryfast)
%     .previewFrame   render a single .png, no video           (unset)
%     .pythonExe      python executable      (auto-detected & cached if unset)
%     .autoInstall    pip-install missing python packages      (false)
%     .forcePythonSearch  re-search, ignoring the cached interpreter  (false)
%     .root/.toolboxes/.masterScript  override component locations
%
%   Python is NOT taken from PATH blindly: every candidate is verified by actually
%   importing numpy/matplotlib/PIL, and the first working one is cached in MATLAB
%   prefs (setpref 'process_gaze'). Candidates: params.pythonExe, the cached choice,
%   MATLAB's pyenv interpreter, python / py -3 / python3, then common install dirs
%   (…\Programs\Python\Python3*, anaconda3, miniconda3, C:\Python3*).
%
%   result is a struct: .status .output .frames .out_fps .body_px .flagged .canon
%
%   Example:
%     process_gaze('C:\data\bird01.csv');
%       p.startTime = 10;
%       p.endTime = 40;
%       p.pCut = 0.6;
%       p.monoFOV = 160;
%     r = process_gaze('C:\data\bird01.csv', p);
%
% Jesus. rev. 13.07.2026

    if nargin < 2 || isempty(params), params = struct(); end
    assert(isfile(csvPath), 'process_gaze:csv', 'CSV not found: %s', csvPath);

    % locate components relative to this file: <root>/functions/video/process_gaze.m
    % thisDir = pwd;
    root = 'C:\Code\Gaze estimation';
    toolboxes = getp(params, 'toolboxes', fullfile(root, 'toolboxes'));
    master    = getp(params, 'masterScript', fullfile(root, 'configfiles', 'master_gaze.py'));
    assert(isfile(master), 'process_gaze:master', 'master script not found: %s', master);
    [pythonExe, hasFFmpeg] = resolve_python(params);

    % just report which labels the table contains, then stop
    if isfield(params, 'listParts') && params.listParts
        result = run_master(pythonExe, master, struct('csv', absify(csvPath), 'list_parts', true));
        fprintf('process_gaze: available body parts: %s\n', strjoin(cellstr(result.parts), ', '));
        return;
    end

    fpsIn = getp(params, 'fps', 59.94);
    step  = getp(params, 'downsampleStep', 2);

    % defaults for output / background ---
    [csvDir, csvName] = fileparts(csvPath);
    output     = getp(params, 'output',     fullfile(csvDir, [csvName '_gaze.mp4']));
    background = getp(params, 'background', fullfile(root, 'configfiles', 'HexArena.png'));

    % assemble config (matches master_gaze.py / pose_clean / pose_render) ---
    cfg = struct();
    cfg.csv        = absify(csvPath);
    cfg.background = absify(background);
    cfg.toolboxes  = toolboxes;
    cfg.video_w    = getp(params, 'videoWidth',  1250);
    cfg.video_h    = getp(params, 'videoHeight', 1160);
    if isfield(params, 'outWidth') && ~isempty(params.outWidth)
        cfg.out_width = params.outWidth;      % rendered pixel width (default 912)
    end
    cfg.fps_in     = fpsIn;
    cfg.downsample_step = step;
    cfg.dpi    = getp(params, 'dpi', 120);
    cfg.crf    = getp(params, 'crf', 24);
    cfg.preset = getp(params, 'preset', 'veryfast');

    % must be set before the ffmpeg check below uses it
    wantVideo = logical(getp(params, 'video', true));
    cfg.video = wantVideo;

    % preview still (PNG) vs video (mp4)
    if isfield(params, 'previewFrame')
        cfg.preview_frame = params.previewFrame;
        [od, on] = fileparts(output); output = fullfile(od, [on '.png']);
    elseif wantVideo && ~hasFFmpeg
        error('process_gaze:ffmpeg', ['ffmpeg was not found by matplotlib, so no video can be written.\n' ...
              'Install it, then restart MATLAB:\n' ...
              '    winget install Gyan.FFmpeg        (or)  conda install -c conda-forge ffmpeg\n' ...
              'Tip: params.previewFrame = 0 renders a single PNG and needs no ffmpeg.']);
    end
    cfg.output = absify(output);

    % body parts: which labels to use, and how they map to anatomical roles
    if isfield(params, 'parts') && ~isempty(params.parts)
        pp = params.parts;
        if ischar(pp), pp = {pp}; end          % a bare 'beak' must not split into letters
        cfg.parts = cellstr(pp(:))';
    end
    if isfield(params, 'roles') && ~isempty(params.roles)
        cfg.roles = params.roles;          % struct: .beak='bill', .head='nape', ...
    end
    if isfield(params, 'bodyPx') && ~isempty(params.bodyPx), cfg.body_px = params.bodyPx; end

    % extra output: derived features .mat + summary figure (true, or an explicit path)
    featOpt = getp(params, 'features', getp(params, 'headDirection', []));   % old name still works
    if ~isempty(featOpt)
        if islogical(featOpt) || isnumeric(featOpt)
            cfg.features = logical(featOpt);        % accepts true or 1
        else
            cfg.features = absify(featOpt);         % an explicit path
        end
        figOpt = getp(params, 'featuresFigure', true);
        if islogical(figOpt) || isnumeric(figOpt)
            cfg.features_figure = logical(figOpt);
        else
            cfg.features_figure = absify(figOpt);
        end
    end
    if ~wantVideo && ~isfield(cfg, 'features') && ~isfield(cfg, 'preview_frame')
        error('process_gaze:nothingToDo', ['params.video = false but nothing else was ' ...
              'requested, so there would be no output.\nSet params.features = true ' ...
              '(or params.previewFrame) as well.']);
    end

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

    % write config + run python master ---
    result = run_master(pythonExe, master, cfg);
    if ~isempty(getp2(result, 'output', ''))
        fprintf('process_gaze: wrote %s  (%d frames @ %.2f fps)\n', ...
                result.output, result.frames, result.out_fps);
    end
    if ~isempty(getp2(result, 'features', ''))
        fprintf('process_gaze: features -> %s\n', result.features);
    end
    if ~isempty(getp2(result, 'figure', ''))
        fprintf('process_gaze: figure   -> %s\n', result.figure);
    end
end

function result = run_master(pythonExe, master, cfg)
    tmp = [tempname '.json'];
    res = [tempname '_result.json'];
    cfg.result_json = res;                 % python writes the outcome here
    fid = fopen(tmp, 'w');  assert(fid > 0, 'cannot open temp config');
    fwrite(fid, jsonencode(cfg));  fclose(fid);

    % -u = unbuffered python, and system() WITHOUT an output argument, so the
    % [clean]/[render] progress streams to the command window as it happens
    % instead of arriving in one lump at the end.
    cmd = sprintf('%s -u "%s" "%s"', qexe(pythonExe), master, tmp);
    if ispc && startsWith(strtrim(cmd), '"')
        cmd = ['"' cmd '"'];   % only when the exe itself needed quoting (path with spaces)
    end
    fprintf('process_gaze: running\n  %s\n', cmd);
    st = system(cmd);
    if isfile(tmp), delete(tmp); end

    if isfile(res)
        result = jsondecode(fileread(res));
        delete(res);
    else
        result = struct('status', 'error', ...
                        'message', sprintf('no result file written (exit %d) - see output above', st));
    end
    if st ~= 0 || ~strcmp(getp2(result, 'status', 'error'), 'ok')
        error('process_gaze:failed', 'pipeline failed: %s', getp2(result, 'message', 'see output above'));
    end
end

%  helpers
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

% ---------------- python discovery ----------------
% "No module named numpy" means python ran but was the WRONG interpreter.
% So we never trust bare `python`: every candidate is verified by actually
% importing the packages, and the winner is cached in MATLAB prefs.
function [exe, hasFFmpeg] = resolve_python(params)
    req = 'numpy, matplotlib, pillow';

    % 1) explicit override always wins
    if isfield(params, 'pythonExe') && ~isempty(params.pythonExe)
        exe = params.pythonExe;
        [ok, hasFFmpeg] = probe_python(exe);
        if ok, return; end
        error('process_gaze:python', ['pythonExe "%s" cannot import %s.\nInstall them:\n    %s -m pip install %s'], ...
              exe, req, qexe(exe), req);
    end

    % 2) interpreter cached from a previous successful run
    if ~(isfield(params, 'forcePythonSearch') && params.forcePythonSearch)
        cached = getpref('process_gaze', 'pythonExe', '');
        if ~isempty(cached)
            [ok, hasFFmpeg] = probe_python(cached);
            if ok, exe = cached; return; end
        end
    end

    % 3) search, verifying the imports for each candidate
    cands = python_candidates();
    runnable = {};
    for i = 1:numel(cands)
        [ok, ff, runs] = probe_python(cands{i});
        if runs, runnable{end+1} = cands{i}; end %#ok<AGROW>
        if ok
            exe = cands{i}; hasFFmpeg = ff;
            setpref('process_gaze', 'pythonExe', exe);
            fprintf('process_gaze: using python -> %s\n', exe);
            return;
        end
    end

    % 4) optional auto-install into the first python that at least runs
    if ~isempty(runnable) && isfield(params, 'autoInstall') && params.autoInstall
        fprintf('process_gaze: installing %s into %s ...\n', req, runnable{1});
        system(sprintf('%s -m pip install %s', qexe(runnable{1}), req));
        [ok, ff] = probe_python(runnable{1});
        if ok
            exe = runnable{1}; hasFFmpeg = ff;
            setpref('process_gaze', 'pythonExe', exe);
            return;
        end
    end

    % 5) nothing worked -> actionable message
    if isempty(runnable)
        msg = sprintf(['No working Python found (tried %d candidates).\n' ...
            'Install Python 3 from https://www.python.org/downloads/ (tick "Add python.exe to PATH"),\n' ...
            'then:  p.pythonExe = ''C:\\path\\to\\python.exe''; process_gaze(csv, p)'], numel(cands));
    else
        msg = sprintf(['Python found, but it is missing %s.\nEasiest fix - run this once:\n' ...
            '    %s -m pip install %s\n' ...
            'Or let MATLAB do it:  p.autoInstall = true; process_gaze(csv, p)\n' ...
            'Or point at another interpreter:  p.pythonExe = ''...\\python.exe'''], ...
            req, qexe(runnable{1}), req);
    end
    error('process_gaze:python', '%s', msg);
end

function [ok, hasFFmpeg, runs] = probe_python(exe)
    code = ['import numpy,matplotlib,PIL;import matplotlib.animation as A;' ...
            'print(''GAZE_OK'',A.FFMpegWriter.isAvailable())'];
    [st, out] = system(sprintf('%s -c "%s"', qexe(exe), code));
    ok        = (st == 0) && contains(out, 'GAZE_OK');
    runs      = ok || contains(out, 'Traceback') || contains(out, 'ModuleNotFoundError');
    hasFFmpeg = ok && contains(out, 'GAZE_OK True');
end

function c = python_candidates()
    c = {};
    try                                     % MATLAB's own configured interpreter
        pe = pyenv;
        if strlength(pe.Executable) > 0, c{end+1} = char(pe.Executable); end
    catch
    end
    if ispc
        c = [c, {'python', 'py -3', 'python3'}];
        pats = { fullfile(getenv('LOCALAPPDATA'), 'Programs', 'Python', 'Python3*', 'python.exe'), ...
                 fullfile(getenv('USERPROFILE'), 'anaconda3', 'python.exe'), ...
                 fullfile(getenv('USERPROFILE'), 'miniconda3', 'python.exe'), ...
                 fullfile(getenv('LOCALAPPDATA'), 'anaconda3', 'python.exe'), ...
                 'C:\ProgramData\Anaconda3\python.exe', ...
                 'C:\Python3*\python.exe' };
    else
        c = [c, {'python3', 'python'}];
        pats = {'/usr/local/bin/python3', '/opt/homebrew/bin/python3', '/usr/bin/python3'};
    end
    for i = 1:numel(pats)
        d = dir(pats{i});
        for j = 1:numel(d)
            if ~d(j).isdir, c{end+1} = fullfile(d(j).folder, d(j).name); end %#ok<AGROW>
        end
    end
    c = unique(c, 'stable');
end

function s = qexe(exe)
    exe = strtrim(char(exe));
    isPath = contains(exe, '\') || contains(exe, '/') || contains(exe, ':');
    if isPath && any(exe == ' ') && ~startsWith(exe, '"')
        s = ['"' exe '"'];      % quote real paths with spaces; leave `python` / `py -3` bare
    else
        s = exe;
    end
end
