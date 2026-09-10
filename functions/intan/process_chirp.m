function result = process_chirp(input, opt)
%PROCESS_CHIRP  Quick pre-sorter survey of a session's raw INTAN channels.
%   result = PROCESS_CHIRP(input, opt)
%
%   MATLAB wrapper: it prepares a JSON config from the session's paths and
%   opt.chirp.* settings and runs the Python master pipeline
%   (../../configfiles/master_chirp.py), which drives the CHIRP engine in
%   ../../toolboxes/CHIRP. Detection, clustering and the statistics live
%   entirely in Python.
%
%   CHIRP reads the raw one-file-per-channel .dat files BEFORE any pipeline
%   conversion or filtering, band-passes them itself, detects spikes at
%   -opt.chirp.negK sigma, splits them by trough amplitude and reports one row
%   per channel x window x cluster. The point is to describe what Kilosort is
%   about to be handed, in the units Kilosort's own settings are written in:
%   trough-to-peak and half-width size the template window, sigma and the
%   cross-channel sharing measure flag channels worth dropping from the map,
%   and the SNR distribution says whether the detection thresholds will be
%   comfortable. It is not a spike sorter and does not pretend to be one.
%
%   INPUTS:
%     input  - NGL input struct after prepforsession. Uses:
%                .sessions(run(1)).info.amplifier_sample_rate  acquisition rate
%                .sessions(run(1)).info.fileformat             must be 'fileperch'
%                .areaMap (optional)                           multi-area groups
%                .run                                          [subject session]
%     opt    - NGL opt struct after prepforsession. Uses:
%                .PathRaw            session raw folder (holds amp*.dat)
%                .FolderProcDataMat  preprocessing output folder
%                .chirp.*            settings, see optSchema / NGL_SetAndRunMe
%
%   OUTPUT (per session, in <FolderProcDataMat>\chirp\):
%     chirp_cluster_stats.csv  one row per area x channel x window x cluster
%     chirp_report.txt         the same run summarised per area
%     *.mp4                    only when opt.chirp.video is true
%
%   result is a struct: .status .csv .report .videos .n_rows .n_channels
%                       .groups .isolated .noise .chirp_version
%   plus .table, the CSV read back as a MATLAB table.
%
%   MULTI-AREA: when input.areaMap is present each area is surveyed as its own
%   channel group. CHIRP's noise test calls a waveform common when it appears
%   on more than max(8, 0.20 x N) of the N channels it saw, so the group it is
%   measured over decides the verdict; per-area grouping judges an area against
%   its own channel count rather than against the whole headstage.
%
%   For an ad-hoc check outside the pipeline, CHIRP's own CLI is simpler:
%       python toolboxes\CHIRP\dat_to_stats.py --folder <raw> --all
%
% Jesus. rev. 10.09.2026

    result = struct('status', 'skipped');

    % --- gate: CHIRP reads one-file-per-channel INTAN only ------------------
    fmt = input.sessions(input.run(1)).info.fileformat;
    if ~strcmp(fmt, 'fileperch')
        warning('process_chirp:format', ...
            ['CHIRP needs INTAN ''fileperch'' (amp-*.dat, one file per ' ...
             'channel); this session is ''%s''. Skipping.'], fmt);
        return
    end
    assert(isfolder(opt.PathRaw), 'process_chirp:raw', ...
        'raw folder not found: %s', opt.PathRaw);

    files = dir(fullfile(opt.PathRaw, 'amp*.dat'));
    if isempty(files)
        warning('process_chirp:noFiles', ...
            'no amp*.dat in %s. Skipping.', opt.PathRaw);
        return
    end
    [~, order] = sort({files.name});          % same order the .bin is built in
    files = files(order);
    fullPaths = arrayfun(@(f) fullfile(f.folder, f.name), files, ...
                         'UniformOutput', false);

    % locate components relative to THIS file: <root>/functions/intan/process_chirp.m
    here      = fileparts(mfilename('fullpath'));
    root      = fileparts(fileparts(here));
    toolboxes = getp(opt.chirp, 'toolboxes', fullfile(root, 'toolboxes', 'CHIRP'));
    master    = getp(opt.chirp, 'masterScript', ...
                     fullfile(root, 'configfiles', 'master_chirp.py'));
    assert(isfile(master), 'process_chirp:master', ...
        'master script not found: %s', master);

    wantVideo = logical(getp(opt.chirp, 'video', false));
    [pythonExe, hasFFmpeg] = resolve_python(opt.chirp);
    if wantVideo && ~hasFFmpeg
        warning('process_chirp:ffmpeg', ...
            ['ffmpeg was not found, so no video can be written. Continuing ' ...
             'with statistics only.\nInstall it with:  winget install Gyan.FFmpeg']);
        wantVideo = false;
    end

    % --- channel groups: one per area, or one for the lot -------------------
    groups = buildGroups(input, opt, fullPaths);

    % --- sample rate: from the session header, not a default ---------------
    % opt.sampleRate is only set later, by the INTAN wrapper; at this point in
    % NGL01 the header read by findSetting is the authoritative source.
    info = input.sessions(input.run(1)).info;
    assert(isfield(info, 'amplifier_sample_rate') && ~isempty(info.amplifier_sample_rate), ...
        'process_chirp:fs', ...
        'amplifier_sample_rate missing from session info; cannot survey.');
    fs = double(info.amplifier_sample_rate);

    % --- assemble config ----------------------------------------------------
    outDir = fullfile(opt.FolderProcDataMat, 'chirp');
    cfg = struct();
    cfg.toolboxes  = toolboxes;
    cfg.groups     = groups;
    cfg.out_dir    = outDir;
    cfg.fs         = fs;
    cfg.duration   = getp(opt.chirp, 'duration',   10);
    cfg.band       = getp(opt.chirp, 'band',       [450 8000]);
    cfg.neg_k      = getp(opt.chirp, 'negK',       5);
    cfg.pos_k      = getp(opt.chirp, 'posK',       8);
    cfg.artifact_k = getp(opt.chirp, 'artifactK',  18);
    cfg.max_k      = getp(opt.chirp, 'maxClusters', 3);
    cfg.segments   = getp(opt.chirp, 'segments',   3);
    cfg.scan_step  = getp(opt.chirp, 'scanStep',   60);
    cfg.video      = wantVideo;
    cfg.video_top  = getp(opt.chirp, 'videoTop',   0);
    cfg.label      = sprintf('%s / %s', input.subjects(input.run(1)).name, ...
                             opt.SavFileName);
    scanRange = getp(opt.chirp, 'scanRange', []);
    if ~isempty(scanRange), cfg.scan_range = scanRange; end
    startAt = getp(opt.chirp, 'start', []);
    if ~isempty(startAt), cfg.start = startAt; end

    fprintf('\nCHIRP: %d channel(s), %d group(s), %g s windows at %g Hz -> %s\n', ...
            numel(fullPaths), numel(groups), cfg.duration, fs, outDir);

    % --- run ----------------------------------------------------------------
    result = run_master(pythonExe, master, cfg);

    % Read the table back so callers can filter without touching the CSV.
    result.table = table();
    if isfield(result, 'csv') && isfile(result.csv)
        try
            result.table = readtable(result.csv);
        catch ME
            warning('process_chirp:readtable', ...
                'stats CSV written but could not be read back: %s', ME.message);
        end
    end
    fprintf(['CHIRP: %d cluster row(s) over %d channel(s). ' ...
             '%d channel(s) look isolated, %d look like noise.\n'], ...
            getp2(result, 'n_rows', 0), getp2(result, 'n_channels', 0), ...
            getp2(result, 'isolated', 0), getp2(result, 'noise', 0));
    fprintf('CHIRP: %s\n', getp2(result, 'report', '(no report)'));
end

% ---------------- channel grouping ----------------
function groups = buildGroups(input, opt, fullPaths) %#ok<INUSD>
% One group per area when a multi-area map exists, otherwise a single group.
% areaMap.chanIdx_per_area holds 1-based indices into the channel ordering,
% which is the alphabetical amp*.dat order the .bin is written from.
%
% Returns a CELL array of structs, not a struct array: jsonencode turns a
% 1x1 struct into a JSON object but a cell into a JSON array, and the master
% script iterates groups either way only if it is an array. The single group
% is named 'all' rather than left blank so the CSV's area column always holds
% a real label; readtable turns an empty one into <missing>.
    if ~isfield(input, 'areaMap') || isempty(input.areaMap)
        groups = {oneGroup('all', fullPaths)};
        return
    end

    am     = input.areaMap;
    maxIdx = max(cellfun(@max, am.chanIdx_per_area));
    if maxIdx > numel(fullPaths)
        % Happens when reduceChanMap trimmed the run but the map still
        % describes the full probe. Grouping would silently mis-assign
        % channels, so survey everything as one group instead.
        warning('process_chirp:areaMismatch', ...
            ['areaMap covers %d channels but only %d amp*.dat files are ' ...
             'present; surveying all channels as one group.'], ...
            maxIdx, numel(fullPaths));
        groups = {oneGroup('all', fullPaths)};
        return
    end

    groups = cell(1, numel(am.uniqueAreas));
    for a = 1:numel(am.uniqueAreas)
        idx = am.chanIdx_per_area{a};
        groups{a} = oneGroup(am.uniqueAreas{a}, fullPaths(idx(:)));
    end
end

function g = oneGroup(name, paths)
% Built by assignment, so the cell of paths stays one field value instead of
% being spread into a struct array the way struct('files', {c}) would.
    g = struct();
    g.name  = name;
    g.files = paths(:)';
end

% ---------------- python plumbing ----------------
function result = run_master(pythonExe, master, cfg)
    tmp = [tempname '.json'];
    res = [tempname '_result.json'];
    cfg.result_json = res;                 % python writes the outcome here
    fid = fopen(tmp, 'w');  assert(fid > 0, 'cannot open temp config');
    fwrite(fid, jsonencode(cfg));  fclose(fid);

    % -u = unbuffered python, and system() WITHOUT an output argument, so the
    % scan progress streams to the command window as it happens instead of
    % arriving in one lump at the end.
    cmd = sprintf('%s -u "%s" "%s"', qexe(pythonExe), master, tmp);
    if ispc && startsWith(strtrim(cmd), '"')
        cmd = ['"' cmd '"'];   % only when the exe itself needed quoting
    end
    fprintf('process_chirp: running\n  %s\n', cmd);
    st = system(cmd);
    if isfile(tmp), delete(tmp); end

    if isfile(res)
        result = jsondecode(fileread(res));
        delete(res);
    else
        result = struct('status', 'error', 'message', sprintf( ...
            'no result file written (exit %d) - see output above', st));
    end
    if st ~= 0 || ~strcmp(getp2(result, 'status', 'error'), 'ok')
        error('process_chirp:failed', 'CHIRP survey failed: %s', ...
              getp2(result, 'message', 'see output above'));
    end
end

function v = getp(p, f, d)
    if isstruct(p) && isfield(p, f) && ~isempty(p.(f)), v = p.(f); else, v = d; end
end

function v = getp2(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end

% "No module named numpy" means python ran but was the WRONG interpreter.
% So we never trust bare `python`: every candidate is verified by actually
% importing the packages, and the winner is cached in MATLAB prefs.
function [exe, hasFFmpeg] = resolve_python(params)
    req = 'numpy, scipy, matplotlib';

    % 1) explicit override always wins
    if isstruct(params) && isfield(params, 'pythonExe') && ~isempty(params.pythonExe)
        exe = params.pythonExe;
        [ok, hasFFmpeg] = probe_python(exe);
        if ok, return; end
        error('process_chirp:python', ...
              'pythonExe "%s" cannot import %s.\nInstall them:\n    %s -m pip install %s', ...
              exe, req, qexe(exe), req);
    end

    % 2) interpreter cached from a previous successful run
    if ~(isstruct(params) && isfield(params, 'forcePythonSearch') && params.forcePythonSearch)
        cached = getpref('process_chirp', 'pythonExe', '');
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
            setpref('process_chirp', 'pythonExe', exe);
            fprintf('process_chirp: using python -> %s\n', exe);
            return;
        end
    end

    % 4) optional auto-install into the first python that at least runs
    if ~isempty(runnable) && isstruct(params) && isfield(params, 'autoInstall') ...
            && params.autoInstall
        fprintf('process_chirp: installing %s into %s ...\n', req, runnable{1});
        system(sprintf('%s -m pip install %s', qexe(runnable{1}), req));
        [ok, ff] = probe_python(runnable{1});
        if ok
            exe = runnable{1}; hasFFmpeg = ff;
            setpref('process_chirp', 'pythonExe', exe);
            return;
        end
    end

    % 5) nothing worked -> actionable message
    if isempty(runnable)
        msg = sprintf(['No working Python found (tried %d candidates).\n' ...
            'Install Python 3 from https://www.python.org/downloads/ (tick "Add python.exe to PATH"),\n' ...
            'then:  opt.chirp.pythonExe = ''C:\\path\\to\\python.exe'''], numel(cands));
    else
        msg = sprintf(['Python found, but it is missing %s.\nEasiest fix - run this once:\n' ...
            '    %s -m pip install %s\n' ...
            'Or let MATLAB do it:  opt.chirp.autoInstall = true;\n' ...
            'Or point at another interpreter:  opt.chirp.pythonExe = ''...\\python.exe'''], ...
            req, qexe(runnable{1}), req);
    end
    error('process_chirp:python', '%s', msg);
end

function [ok, hasFFmpeg, runs] = probe_python(exe)
    code = ['import numpy,scipy,matplotlib;import matplotlib.animation as A;' ...
            'print(''CHIRP_OK'',A.FFMpegWriter.isAvailable())'];
    [st, out] = system(sprintf('%s -c "%s"', qexe(exe), code));
    ok        = (st == 0) && contains(out, 'CHIRP_OK');
    runs      = ok || contains(out, 'Traceback') || contains(out, 'ModuleNotFoundError');
    hasFFmpeg = ok && contains(out, 'CHIRP_OK True');
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
        s = ['"' exe '"'];      % quote real paths with spaces; leave `python` bare
    else
        s = exe;
    end
end
