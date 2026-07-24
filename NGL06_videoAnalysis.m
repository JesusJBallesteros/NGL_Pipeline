%% NGL06_videoAnalysis. Video-based gaze analysis over per-session DLC CSVs.
%
% PURPOSE:
%   Stage 6: fans out the GazEstim Python pipeline (pose_clean -> features ->
%   pose_render, driven by configfiles/master_gaze.py) across every discovered
%   (subject, session) DLC csv. Per csv it can produce, next to the csv:
%     - a gaze-cone-overlaid mp4 (or a single PNG when opt.gaze.previewFrame
%       is set),
%     - derived per-frame kinematic features as a .mat plus a summary .png
%       figure (opt.gaze.features), and can run features-only with the video
%       off (opt.gaze.video = false) for a much faster pass.
%   Body-part selection (opt.gaze.parts) and role remap (opt.gaze.roles.*) are
%   forwarded too. There is no area loop and no aggregated data - the stage is
%   a pure fan-out over the on-disk csv list.
%
% USAGE:
%   Do NOT run or edit this script directly. Configure via your project's
%   NGL_SetAndRunMe.m (section 6) and invoke it from there.
%
% REQUIRED WORKSPACE VARIABLES (set by NGL_SetAndRunMe -> NGL00_Prep):
%   datadrive, studyname, subjects, dates, opt
%
% PROJECT-LOCAL COPIES REQUIRED at <input.analysisCode>/:
%   master_gaze.py     copy once from <toolbox>/configfiles/master_gaze.py
%   HexArena.png       copy once from <toolbox>/configfiles/HexArena.png
%   Rationale: per-project tweaks (custom arena background, patched
%   master) must be visible in analysisCode. NGL06 hard-errors with the
%   copy-from location if either file is missing.
%
% INPUT LAYOUT:
%   <input.bhvfolder>/<subject>/<session>/*.csv    DLC pose table(s)
%   Per-session convention: exactly ONE csv. Zero csv -> warn + skip.
%   Multiple csvs -> NGL06:multipleCsvs error for THAT session (batch
%   continues; user disambiguates by keeping one csv per session).
%
% PIPELINE:
%   00.  NGL00_Prep + set_default + findSessions
%   00b. Gate on opt.gaze.do
%   01.  Resolve toolbox / analysisCode / python paths + smoke-test
%   02.  Build the params struct once from opt.gaze.* (incl. video / features /
%        parts / roles); fail fast if the config would produce no output.
%   03.  Loop subjects x sessions:
%          - discover csv, decide skip / process / halt-per-session
%          - process_gaze() -> mp4/png and/or features .mat + figure next to csv
%          - failures are caught: warning + per-csv <name>_gazeFailed.txt
%   04.  Print summary (nProcessed / nFeatures / nSkippedExists / nSkippedNoCsv /
%        nSkippedMultiCsv / nFailed).
%
% OUTPUT (per successful csv; any subset, depending on opt.gaze.*):
%   <csv_folder>/<csvName>_gaze.mp4                      (video; ..._gaze.png if previewFrame)
%   <csv_folder>/<csvName>_gaze_estimated_features.mat   (if opt.gaze.features)
%   <csv_folder>/<csvName>_gaze_estimated_features.png   (features summary figure)
%
% DEPENDENCIES:
%   functions/video/process_gaze.m (this wrapper)
%   toolboxes/GazEstim/{pose_clean, features, pose_render, make_plate}.py
%   configfiles/master_gaze.py     (template; runtime copy at analysisCode)
%   configfiles/HexArena.png       (template; runtime copy at analysisCode)
%
% SEE ALSO:
%   docs/gaze_pipeline.md, docs/examples/gaze/ (bundled smoke-test csv +
%   reference mp4).
%
% Last modified 24.07.2026 (Jesus) - forward opt.gaze.video / features / parts /
%                                     roles; features-only mode; per-mode
%                                     skip-if-exists (mp4 / png /
%                                     _estimated_features.mat); capture + report
%                                     produced files. Tracks GazEstim features.py.

%% 00. Standard scaffolding.
NGL00_Prep
[input, opt] = set_default(input, opt);
input.sessions = findSessions(input);

%% 00b. Early exit if the master gate is off.
if ~isfield(opt,'gaze') || ~isfield(opt.gaze,'do') || ~opt.gaze.do
    warning('NGL06:nothingToDo', ...
        'opt.gaze.do is false or missing; NGL06_videoAnalysis has nothing to do.');
    return
end

%% 01. Resolve paths and smoke-test Python.
toolboxDir   = fullfile(input.toolbox, 'toolboxes', 'GazEstim');
templateDir  = fullfile(input.toolbox, 'configfiles');
masterScript = localResolvePath(opt.gaze.masterScript, ...
                    fullfile(input.analysisCode, 'master_gaze.py'));
background   = localResolvePath(opt.gaze.background, ...
                    fullfile(input.analysisCode, 'HexArena.png'));

assert(isfolder(toolboxDir), 'NGL06:noToolbox', ...
    ['GazEstim toolbox folder missing at %s. The toolbox ships it under ', ...
     'toolboxes/GazEstim/; check your input.toolbox (currently %s).'], ...
    toolboxDir, input.toolbox);

if ~isfile(masterScript)
    error('NGL06:noMaster', ...
        ['master_gaze.py not found at\n   %s\n', ...
         'This must live under your project''s analysisCode/. Copy the ', ...
         'template once:\n   copy "%s\\master_gaze.py"  "%s\\master_gaze.py"'], ...
        masterScript, templateDir, input.analysisCode);
end
if ~isfile(background)
    error('NGL06:noBackground', ...
        ['background image not found at\n   %s\n', ...
         'This must live under your project''s analysisCode/. Copy the ', ...
         'template once:\n   copy "%s\\HexArena.png"  "%s\\HexArena.png"\n', ...
         'For a custom arena, generate a new plate via ', ...
         '%s\\make_plate.py photo.png clean_plate.png x0 y0 x1 y1'], ...
        background, templateDir, input.analysisCode, toolboxDir);
end

% Python executable: opt override -> input.GAZEpythonExe -> 'python'.
if isfield(opt.gaze,'pythonExe') && ~isempty(opt.gaze.pythonExe)
    pythonExe = opt.gaze.pythonExe;
elseif isfield(input,'GAZEpythonExe') && ~isempty(input.GAZEpythonExe)
    pythonExe = input.GAZEpythonExe;
else
    pythonExe = 'python';
end

% Fail-fast smoke test: python starts, imports the deps master_gaze.py
% needs, and ffmpeg is on PATH. Cheaper to fail here than to fail per
% session with a cryptic subprocess error.
localSmokeTestPython(pythonExe, toolboxDir);

%% 02. Build the params struct once (mirrored to process_gaze's signature).
params = struct();
params.pythonExe    = pythonExe;
params.toolboxes    = toolboxDir;
params.masterScript = masterScript;
params.background   = background;

% Map the schema-defined opt.gaze fields into MATLAB-side param names.
% Fields with empty defaults stay unset so process_gaze's own defaults
% apply. This keeps NGL06 as a thin fan-out rather than a re-writer of
% the wrapper's contract.
fwdList = { ...
    'fps',            'fps'; ...
    'downsampleStep', 'downsampleStep'; ...
    'targetFps',      'targetFps'; ...
    'startTime',      'startTime'; ...
    'endTime',        'endTime'; ...
    'maxFrames',      'maxFrames'; ...
    'maxSeconds',     'maxSeconds'; ...
    'pCut',           'pCut'; ...
    'devFac',         'devFac'; ...
    'smooth',         'smooth'; ...
    'wMed',           'wMed'; ...
    'boneTolFrac',    'boneTolFrac'; ...
    'boneTolMad',     'boneTolMad'; ...
    'orderMargin',    'orderMargin'; ...
    'videoWidth',     'videoWidth'; ...
    'videoHeight',    'videoHeight'; ...
    'drawCones',      'gaze'; ...
    'monoFOV',        'monoFOV'; ...
    'binoHalf',       'binoHalf'; ...
    'coneMult',       'coneMult'; ...
    'eyeFwdFrac',     'eyeFwdFrac'; ...
    'eyeLatFrac',     'eyeLatFrac'; ...
    'dpi',            'dpi'; ...
    'crf',            'crf'; ...
    'preset',         'preset'; ...
    'previewFrame',   'previewFrame' };
for k = 1:size(fwdList, 1)
    optField   = fwdList{k, 1};
    paramField = fwdList{k, 2};
    if isfield(opt.gaze, optField) && ~isempty(opt.gaze.(optField))
        params.(paramField) = opt.gaze.(optField);
    end
end

% --- New GazEstim outputs / toggles (features.py, video off, parts, roles) ---
% Handled explicitly rather than via the generic forwarder above: roles is a
% nested struct (drop the empty identity leaves before forwarding), and
% features/video are booleans the generic '~isempty' rule would forward even
% when false.
if isfield(opt.gaze, 'video')
    params.video = logical(opt.gaze.video);          % false -> features-only (skip mp4)
end
featuresOn = false;
if isfield(opt.gaze, 'features')
    ff = opt.gaze.features;
    if (islogical(ff) && isscalar(ff) && ff) || (ischar(ff) && ~isempty(ff))
        params.features = ff;                        % true, or an explicit .mat path
        featuresOn = true;                           % process_gaze adds the summary figure by default
    end
end
if isfield(opt.gaze, 'parts') && ~isempty(opt.gaze.parts)
    params.parts = opt.gaze.parts;                   % cellstr subset of DLC labels
end
if isfield(opt.gaze, 'roles') && isstruct(opt.gaze.roles)
    rr = localCompactStruct(opt.gaze.roles);         % keep only non-empty role->label maps
    if ~isempty(fieldnames(rr)), params.roles = rr; end
end

% Which primary file marks a session as already done (drives skip-if-exists)?
%   previewFrame -> _gaze.png ; video on -> _gaze.mp4 ;
%   video off (features-only) -> _gaze_estimated_features.mat
wantPreview = isfield(opt.gaze,'previewFrame') && ~isempty(opt.gaze.previewFrame);
videoOff    = isfield(opt.gaze,'video') && islogical(opt.gaze.video) && ~opt.gaze.video;
if wantPreview
    primarySuffix = '_gaze.png';
elseif ~videoOff
    primarySuffix = '_gaze.mp4';
else
    primarySuffix = '_gaze_estimated_features.mat';
end

% Fail fast: a config asking for no video, no features and no preview would make
% process_gaze error once per session. Catch it up front (mirrors the other
% NGL06 pre-flight checks).
if videoOff && ~featuresOn && ~wantPreview
    error('NGL06:nothingToDo', ...
        ['opt.gaze.video = false but neither opt.gaze.features nor ', ...
         'opt.gaze.previewFrame is set: nothing would be produced. Enable ', ...
         'opt.gaze.features (features-only pass) or set opt.gaze.previewFrame.']);
end

skipIfExists = ~(isfield(opt.gaze,'overwrite') && opt.gaze.overwrite);

fprintf('\nNGL06_videoAnalysis: python=%s\n', pythonExe);
fprintf('  toolboxes = %s\n', toolboxDir);
fprintf('  master    = %s\n', masterScript);
fprintf('  background= %s\n', background);
fprintf('  skipIfExists = %s\n', localBool2Str(skipIfExists));
fprintf('  video = %s ; features = %s ; primary = *%s\n', ...
        localBool2Str(~videoOff), localBool2Str(featuresOn), primarySuffix);

%% 03. Fan out over (subject, session).
nProcessed        = 0;
nFeatures         = 0;
nSkippedExists    = 0;
nSkippedNoCsv     = 0;
nSkippedMultiCsv  = 0;
nFailed           = 0;
for x = 1:input.nsubjects
    subject = input.subjects(x).name;
    for y = 1:input.sessions(x).nsessions
        session = input.sessions(x).list{y};
        sessDir = fullfile(input.bhvfolder, subject, session);
        if ~isfolder(sessDir)
            fprintf('NGL06: [%s/%s] behaviour folder missing (%s); skipping.\n', ...
                    subject, session, sessDir);
            nSkippedNoCsv = nSkippedNoCsv + 1;
            continue
        end

        hits = dir(fullfile(sessDir, '*.csv'));
        if isempty(hits)
            fprintf('NGL06: [%s/%s] no *.csv under %s; skipping.\n', ...
                    subject, session, sessDir);
            localWriteNote(sessDir, 'gazeSkipped_noCsv.txt', ...
                'No DLC .csv found in this behaviour folder.', subject, session);
            nSkippedNoCsv = nSkippedNoCsv + 1;
            continue
        end
        if numel(hits) > 1
            names = strjoin({hits.name}, ', ');
            warning('NGL06:multipleCsvs', ...
                ['[%s/%s] found %d CSVs (%s). One CSV per session ', ...
                 'is required; NGL06 halts this session but continues ', ...
                 'the batch. Keep exactly one CSV, then re-run.'], ...
                subject, session, numel(hits), names);
            localWriteNote(sessDir, 'gazeSkipped_multiCsv.txt', ...
                sprintf(['Multiple DLC .csv files found (%d): %s\n', ...
                         'Keep exactly one CSV per session and re-run NGL06.'], ...
                        numel(hits), names), subject, session);
            nSkippedMultiCsv = nSkippedMultiCsv + 1;
            continue
        end

        csvPath   = fullfile(hits(1).folder, hits(1).name);
        [~, cnm] = fileparts(csvPath);
        % Primary output for skip-if-exists depends on the mode (set in 02):
        % _gaze.png (preview) / _gaze.mp4 (video) / _gaze_estimated_features.mat.
        expectedOut = fullfile(sessDir, [cnm primarySuffix]);
        if skipIfExists && isfile(expectedOut)
            fprintf('NGL06: [%s/%s] output already exists (%s); skipping.\n', ...
                    subject, session, expectedOut);
            nSkippedExists = nSkippedExists + 1;
            continue
        end

        fprintf('\nNGL06: ===== [%s/%s] %s =====\n', subject, session, hits(1).name);
        try
            result = process_gaze(csvPath, params);
            nProcessed = nProcessed + 1;
            % process_gaze / master_gaze already wrote each artifact next to
            % the csv; report exactly what came back so the log names them.
            if isfield(result,'output') && ~isempty(result.output)
                fprintf('  video    -> %s\n', result.output);
            end
            if isfield(result,'features') && ~isempty(result.features)
                fprintf('  features -> %s\n', result.features);
                nFeatures = nFeatures + 1;
            end
            if isfield(result,'figure') && ~isempty(result.figure)
                fprintf('  figure   -> %s\n', result.figure);
            end
        catch ME
            warning('NGL06:gazeFailed', ...
                '[%s/%s] gaze pipeline failed: %s', subject, session, ME.message);
            failFile = fullfile(sessDir, [cnm '_gazeFailed.txt']);
            localWriteNote(sessDir, [cnm '_gazeFailed.txt'], ...
                sprintf('process_gaze failed on %s\n\n%s\n\nStack:\n%s', ...
                        csvPath, ME.message, getReport(ME, 'extended', 'hyperlinks', 'off')), ...
                subject, session);
            fprintf('  wrote %s\n', failFile);
            nFailed = nFailed + 1;
            continue
        end
    end
end

%% 04. Summary.
fprintf(['\n=== NGL06_videoAnalysis summary ===\n', ...
         '  processed        : %d\n', ...
         '  features saved   : %d\n', ...
         '  skipped (exists) : %d\n', ...
         '  skipped (no csv) : %d\n', ...
         '  skipped (multi)  : %d\n', ...
         '  failed           : %d\n'], ...
        nProcessed, nFeatures, nSkippedExists, nSkippedNoCsv, nSkippedMultiCsv, nFailed);


%% ===================================================================
%% Local helpers.
%% ===================================================================
function p = localResolvePath(userVal, defaultVal)
% Empty user override -> use default; else use user override verbatim.
    if isempty(userVal), p = defaultVal; else, p = char(userVal); end
end

function s = localBool2Str(b)
    if b, s = 'true'; else, s = 'false'; end
end

function s = localCompactStruct(in)
% Keep only the fields of `in` whose value is a non-empty char. Used to strip
% the empty identity leaves from opt.gaze.roles before forwarding them to
% process_gaze (which expects .roles = struct of role->label maps).
    s = struct();
    if ~isstruct(in), return; end
    fn = fieldnames(in);
    for i = 1:numel(fn)
        v = in.(fn{i});
        if ischar(v) && ~isempty(v)
            s.(fn{i}) = v;
        end
    end
end

function localWriteNote(dstDir, fname, body, subject, session)
% Write a small diagnostic text file next to the session data.
    if ~isfolder(dstDir), return; end
    fpath = fullfile(dstDir, fname);
    fid = fopen(fpath, 'w');
    if fid < 0, return; end
    fprintf(fid, 'NGL06_videoAnalysis note\n');
    fprintf(fid, 'subject : %s\n', subject);
    fprintf(fid, 'session : %s\n', session);
    fprintf(fid, 'time    : %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
    fprintf(fid, '%s\n', body);
    fclose(fid);
end

function localSmokeTestPython(pythonExe, toolboxDir)
% Fail-fast pre-flight: python starts, master_gaze deps import, ffmpeg
% is on PATH. If any step trips, raise a single NGL06:pythonEnv error
% that names the missing piece and the fix.
    scriptLines = { ...
        'import sys, shutil, importlib', ...
        'missing = []', ...
        'for m in ("numpy","matplotlib","PIL"):', ...
        '    try:', ...
        '        importlib.import_module(m)', ...
        '    except Exception as e:', ...
        '        missing.append(m + ": " + str(e))', ...
        'if shutil.which("ffmpeg") is None:', ...
        '    missing.append("ffmpeg: not found on PATH")', ...
        'if missing:', ...
        '    print("MISSING " + " | ".join(missing)); sys.exit(2)', ...
        'print("OK")' };
    scriptFile = [tempname '.py'];
    fid = fopen(scriptFile, 'w');
    assert(fid > 0, 'NGL06:tempPy', 'cannot write temp python probe');
    fwrite(fid, strjoin(scriptLines, sprintf('\n')));
    fclose(fid);
    cmd = sprintf('"%s" "%s"', pythonExe, scriptFile);
    [st, out] = system(cmd);
    if isfile(scriptFile), delete(scriptFile); end
    ln = strtrim(regexprep(out, '.*(OK|MISSING [^\r\n]*).*', '$1'));
    if st == 0 && startsWith(ln, 'OK')
        fprintf('NGL06: python smoke test OK (%s).\n', pythonExe);
        return
    end
    if startsWith(ln, 'MISSING ')
        error('NGL06:pythonEnv', ...
            ['Python environment is missing required dependencies.\n   %s\n', ...
             'Fix: install numpy / matplotlib / pillow (pip install ', ...
             'numpy matplotlib pillow) and put ffmpeg on PATH. Python: %s'], ...
            ln(9:end), pythonExe);
    end
    error('NGL06:pythonEnv', ...
        ['Python probe failed (exit %d). Output:\n%s\n', ...
         'Fix: verify pythonExe is a valid interpreter (currently %s). ', ...
         'Set opt.gaze.pythonExe or NGL_machineConfig.GAZEpythonExe to ', ...
         'override.'], st, out, pythonExe);
    %#ok<*NASGU>
    %#ok<*UNRCH>
    %#ok<*ASGLU>
    %#ok<*UNRCH>
end
