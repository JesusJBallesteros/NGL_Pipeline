function savePreprocInfo(input, opt, mode)
% savePreprocInfo  Persist the resolved opt struct used by NGL01 preprocessing.
%
% PURPOSE:
%   Captures the validated, post-set_default options that governed an
%   NGL01_Main preprocessing run. NGL02_postPhy reloads this snapshot so
%   downstream analyses use exactly the options that were used during
%   preprocessing, even after MATLAB has been restarted between stages.
%
% USAGE:
%   savePreprocInfo(input, opt, 'master')   % once, before the session loop
%   savePreprocInfo(input, opt, 'session')  % once per session iteration
%
% INPUTS:
%   input  - struct from set_default. Required fields:
%              .analysisCode    (both modes)
%              .toolbox         (both modes, used for git version probe)
%              .subjects        (session mode)
%              .sessions        (session mode)
%              .run             (session mode) [subj_idx sess_idx]
%   opt    - resolved options struct. Post-set_default for 'master'.
%            Post-prepforsession (carrying session paths) for 'session'.
%   mode   - 'master'  -> writes <analysisCode>\preprocInfo_lastRun.mat,
%                         stripping session-specific path fields from opt.
%            'session' -> writes <opt.preprocessing>\preprocInfo.mat for the
%                         current session (input.run-derived).
%
% OUTPUT (on disk):
%   preprocInfo struct with fields:
%     .opt              the resolved opt that was actually used
%     .Areas            input.Areas (cell array) or []
%     .subject          (session mode only) subject name
%     .session          (session mode only) session name
%     .runDate          datetime when this snapshot was written
%     .toolboxVersion   struct(.hash, .branch) of the toolbox git checkout,
%                       or 'unknown' if git probe fails (e.g. not a git
%                       checkout, git missing, or shell error)
%     .MATLABversion    output of version()
%     .hardwareFilters  (session mode, INTAN) struct from chckV with the
%                       hardware notch / HPF / LPF / DSP cutoff actually
%                       used by the headstage. NaN where unset.
%     .recordingTime    (session mode, INTAN) wall-clock recording start
%                       (settings.xml Date + Time), char.
%     .recordingNotes   (session mode, INTAN) cellstr of non-empty user
%                       notes (Note1/Note2/Note3) from RHX.
%     .rhxVersion       (session mode, INTAN) RHX recorder version string.
%     .controllerType   (session mode, INTAN) RHX controller hardware
%                       (e.g. 'ControllerRecordUSB3').
%     .info             (session mode) full snapshot of
%                       input.sessions(x).info as chckV / recoverInfoForRegen
%                       returned it. Forward-looking: lets future regen runs
%                       read whatever chckV captured without re-scanning raw.
%     .regenFromPreproc (session mode, only when regen path ran) logical true
%                       so a study-tree grep tells you which sessions were
%                       rebuilt from preprocessed-only data.
%     .regenFromSystem  (session mode, only when regen path ran) 'INTAN'
%                       or 'Deuteron'.
%
% NOTES:
%   - Master mode strips session-specific path fields (PathRaw, SavFileName,
%     FolderProcDataMat, behavFiles, spikeSorted, trialSorted, analysis,
%     KSfolder, KSfolders) because they only carry meaning inside a session.
%   - This is a write-only persistence helper. Loading is the responsibility
%     of loadPreprocInfo (consumed by NGL02_postPhy).
%   - Errors here should never abort NGL01; the caller is expected to wrap
%     this in a try/catch if it wants belt-and-braces protection. Internal
%     git probing is already best-effort and never throws.
%
% Last modified 18.06.2026 (Jesus) - session mode now surfaces
%                                     hardwareFilters / recordingTime /
%                                     recordingNotes / rhxVersion from
%                                     input.sessions(x).info (populated
%                                     by chckV from INTAN settings.xml).

%% Validate mode
assert(ischar(mode) && ismember(mode, {'master','session'}), ...
    'NGL:savePreprocInfo', ...
    'mode must be ''master'' or ''session''.');

%% Common metadata
preprocInfo               = struct();
preprocInfo.runDate       = datetime('now');
preprocInfo.MATLABversion = version();

if isfield(input, 'Areas')
    preprocInfo.Areas = input.Areas;
else
    preprocInfo.Areas = [];
end

% Toolbox version via git (best-effort; non-fatal if anything fails).
preprocInfo.toolboxVersion = struct('hash','unknown','branch','unknown');
if isfield(input, 'toolbox') && exist(input.toolbox, 'dir')
    try
        cwd     = pwd;
        cleanup = onCleanup(@() cd(cwd));
        cd(input.toolbox)
        [s1, h] = system('git rev-parse --short HEAD');
        [s2, b] = system('git rev-parse --abbrev-ref HEAD');
        if s1 == 0, preprocInfo.toolboxVersion.hash   = strtrim(h); end
        if s2 == 0, preprocInfo.toolboxVersion.branch = strtrim(b); end
    catch
        % swallow — leave 'unknown'
    end
end

%% Mode-specific payload and destination
switch mode
    case 'master'
        % Strip per-session path fields so the master snapshot represents
        % only study-wide, post-set_default options.
        sessionFields = {'PathRaw','SavFileName','FolderProcDataMat', ...
                         'behavFiles','spikeSorted','trialSorted', ...
                         'analysis','KSfolder','KSfolders'};
        optMaster = opt;
        for k = 1:numel(sessionFields)
            if isfield(optMaster, sessionFields{k})
                optMaster = rmfield(optMaster, sessionFields{k});
            end
        end
        preprocInfo.opt = optMaster;

        assert(isfield(input,'analysisCode') && ~isempty(input.analysisCode), ...
            'NGL:savePreprocInfo', ...
            'Master save requires input.analysisCode to be set.');
        target = fullfile(input.analysisCode, 'preprocInfo_lastRun.mat');

    case 'session'
        preprocInfo.opt = opt;

        assert(isfield(input,'run') && numel(input.run) == 2, ...
            'NGL:savePreprocInfo', ...
            'Session save requires input.run = [subj_idx sess_idx].');
        x = input.run(1); y = input.run(2);

        assert(isfield(input,'subjects') && x <= numel(input.subjects), ...
            'NGL:savePreprocInfo', 'input.subjects index out of range.');
        assert(isfield(input,'sessions') && x <= numel(input.sessions) && ...
               y <= numel(input.sessions(x).list), ...
            'NGL:savePreprocInfo', 'input.sessions index out of range.');

        preprocInfo.subject = input.subjects(x).name;
        preprocInfo.session = input.sessions(x).list{y};

        % Pass-through provenance fields from the session's `info`
        % (populated by chckV in case 3 / INTAN). They're forensic only:
        % future me trying to figure out why a re-analysis differs from
        % the original run can read recording wall-clock, the notes the
        % experimenter typed into RHX, the hardware filter chain, and
        % the RHX recorder version straight from preprocInfo.mat.
        if isfield(input,'sessions') && numel(input.sessions) >= x ...
                && isstruct(input.sessions(x)) && isfield(input.sessions(x),'info')
            sInfo = input.sessions(x).info;
            for f = {'hardwareFilters','recordingTime','recordingNotes','rhxVersion','controllerType'}
                if isfield(sInfo, f{1})
                    preprocInfo.(f{1}) = sInfo.(f{1});
                end
            end
            % Full info-struct snapshot so future regens can read
            % everything chckV produced without rescanning raw. Costs a
            % few hundred bytes and saves a lot of guessing later.
            preprocInfo.info = sInfo;
        end

        % Audit trail for regen runs: stamp whether this snapshot came
        % from a regenFrom.preproc=true execution so a `dir | grep`
        % across the study tells you which sessions were rebuilt.
        if isfield(opt,'regenFrom') && isfield(opt.regenFrom,'preproc') ...
                && opt.regenFrom.preproc
            preprocInfo.regenFromPreproc = true;
            if isfield(opt.regenFrom,'system')
                preprocInfo.regenFromSystem = opt.regenFrom.system;
            end
        end

        assert(isfield(opt,'FolderProcDataMat') && ~isempty(opt.FolderProcDataMat), ...
            'NGL:savePreprocInfo', ...
            'Session save requires opt.preprocessing (set by prepforsession).');
        if ~exist(opt.FolderProcDataMat, 'dir'), mkdir(opt.FolderProcDataMat); end
        target = fullfile(opt.FolderProcDataMat, 'preprocInfo.mat');
end

%% Write
save(target, 'preprocInfo', '-mat');
fprintf('preprocInfo written: %s\n', target);
end
