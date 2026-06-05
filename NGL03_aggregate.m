%% NGL03_aggregate. Cross-session and cross-subject aggregator.
%
% PURPOSE:
%   Stage 3 aggregator: pulls together the per-session outputs that
%   NGL02_postPhy and NGL02_LFP produced and packs them into cell
%   arrays indexed by (subject, session). Two gates, with a dependency:
%       opt.aggregateSessions  -> build per-subject aggregated .mat under
%                                 data\analysis\<subject>\<subject>_aggregated.mat
%       opt.aggregateSubjects  -> build study-level aggregated .mat at
%                                 data\analysis\aggregated.mat
%                                 (requires aggregateSessions=true)
%   Both default to false; the aggregation is opt-in.
%
% USAGE:
%   Do NOT run or edit this script directly. Configure via your project's
%   NGL_SetAndRunMe.m (template section 3.5) and invoke it from there.
%
% REQUIRED WORKSPACE VARIABLES (set by NGL_SetAndRunMe -> NGL00_Prep):
%   datadrive, studyname, subjects, dates, opt
%
% PIPELINE:
%   00.  NGL00_Prep        - parse inputs into input/opt structs
%   00b. Areas recovery    - pull input.Areas from master preprocInfo
%   01.  set_default       - validate opt, build paths
%   02.  findSessions      - discover session folders on disk
%   03.  Build aggregation target list from opt flags
%   04.  Per-subject session aggregation (if opt.aggregateSessions)
%   05.  Cross-subject aggregation       (if opt.aggregateSubjects)
%
% AGGREGATION SHAPE:
%   Cell arrays of size (Nsubj, Nsess) where each cell contains the
%   loaded variable for that (subject, session). Cells corresponding to
%   missing files stay empty.
%   - Per-subject file: (1 x nSess_for_this_subject) cells per variable.
%   - Study-level file: (Nsubj x maxSess_across_subjects) cells, padded
%     with empties where a subject ran fewer sessions than maxSess.
%
% VARIABLES AGGREGATED (this version — SPIKE SIDE ONLY):
%   Always (when opt.doSpikething): spike, neurons, fireRate,
%                                    condition, events, trialdef.
%   Conditional:
%     - neuralDynamics  if opt.popDyn.do
%     - blob            if opt.offlineTrack or opt.useTrack
%
%   LFP-side aggregation is planned but not yet here. The LFP outputs
%   (continuous, per-alignment trial-parsed, artifact-rejected copies)
%   are multi-file per session and need a separate model. The current
%   spike-side template is the foundation; LFP entries will join the
%   varList in a follow-up.
%
% OUTPUTS:
%   - data\analysis\<subject>\<subject>_aggregated.mat
%       (one per subject, when aggregateSessions=true)
%   - data\analysis\aggregated.mat
%       (study-level, when aggregateSubjects=true)
%
% Last modified 29.05.2026 (Jesus) - new script (task #23)

%% 00. Check current inputs.
NGL00_Prep

%% 00b. Pre-flight: recover input.Areas from the NGL01 master snapshot.
% set_default needs input.Areas BEFORE it can build the area map. If the
% user did not re-set it in NGL_SetAndRunMe (common case when NGL03 runs
% in a fresh MATLAB session), pull it from analysisCode\preprocInfo_lastRun.mat.
if ~isfield(input,'Areas') || isempty(input.Areas)
    analysisCodePath = fullfile(input.datadrive, input.studyName, 'analysisCode');
    if ~contains(analysisCodePath, ':\') && ~isempty(input.datadrive)
        analysisCodePath = fullfile([input.datadrive(1) ':\'], input.studyName, 'analysisCode');
    end
    try
        masterInfo = loadPreprocInfo(analysisCodePath, 'master');
        if isfield(masterInfo,'Areas') && ~isempty(masterInfo.Areas)
            input.Areas = masterInfo.Areas;
            fprintf('NGL03: recovered input.Areas = {%s} from preprocInfo_lastRun.mat\n', ...
                    strjoin(input.Areas, ', '));
        end
    catch ME
        if strcmp(ME.identifier, 'NGL:loadPreprocInfo:notFound')
            % No master snapshot. Continue in single-area mode.
        else
            warning('NGL03:preflight', 'Could not read master preprocInfo: %s', ME.message);
        end
    end
end

[input, opt] = set_default(input, opt);

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

%% 02. Early-exit if neither aggregation gate is enabled.
if ~opt.aggregateSessions && ~opt.aggregateSubjects
    warning('NGL03:nothingToDo', ...
        ['Neither opt.aggregateSessions nor opt.aggregateSubjects is true. ', ...
         'NGL03_acrossSession has nothing to do.']);
    return
end

%% 03. Build the list of variables to aggregate from opt flags.
% Each row: {varName, sourceKey, fileName}
%   sourceKey is one of 'analysis' | 'spikeSorted' | 'trialSorted',
%   resolved per-session via input.<sourceKey>/<subject>/<session>.
varList = {};

if isfield(opt,'doSpikething') && opt.doSpikething
    varList = [varList; {
        'spike',     'spikeSorted', 'spike.mat'    ;
        'neurons',   'analysis',    'neurons.mat'  ;
        'fireRate',  'analysis',    'fireRate.mat' ;
        'condition', 'trialSorted', 'condition.mat';
        'events',    'trialSorted', 'events.mat'   ;
        'trialdef',  'trialSorted', 'trialdef.mat' }];
    if isfield(opt,'popDyn') && isfield(opt.popDyn,'do') && opt.popDyn.do
        varList = [varList; {'neuralDynamics', 'analysis', 'neuralDynamics.mat'}];
    end
    if (isfield(opt,'offlineTrack') && opt.offlineTrack) || ...
       (isfield(opt,'useTrack')     && opt.useTrack)
        varList = [varList; {'blob', 'analysis', 'blob.mat'}];
    end
end

% TODO LFP-side aggregation (task: extend NGL03 with LFP outputs once the
% LFP path stabilises). The LFP outputs are multi-file per session
% (per-alignment trial-parsed, continuous, artifact-rejected copies);
% they need a richer aggregation model than the single-file rows above.

if isempty(varList)
    warning('NGL03:emptyVarList', ...
        ['No aggregation targets enabled. Set opt.doSpikething=true to ', ...
         'aggregate spike-side files. (LFP-side aggregation is planned.)']);
    return
end

fprintf('NGL03_acrossSession: will aggregate %d variables: %s\n', ...
        size(varList,1), strjoin(varList(:,1)', ', '));

%% 04. Session aggregation per subject.
if opt.aggregateSessions
    fprintf('NGL03_acrossSession: aggregating sessions per subject...\n');

    for x = 1:input.nsubjects
        subject = input.subjects(x).name;
        nSess   = input.sessions(x).nsessions;

        % Pre-allocate one (1 x nSess) cell per target.
        perSubject = struct();
        for v = 1:size(varList,1)
            perSubject.(['all' varList{v,1}]) = cell(1, nSess);
        end

        % Walk sessions, load each existing file into its cell.
        for y = 1:nSess
            session = input.sessions(x).list{y};
            for v = 1:size(varList,1)
                varName = varList{v,1};
                fname   = varList{v,3};
                folder  = localSessionFolder(input, varList{v,2}, subject, session);
                fpath   = fullfile(folder, fname);
                if ~isfile(fpath), continue, end

                S = load(fpath);
                if isfield(S, varName)
                    payload = S.(varName);
                elseif numel(fieldnames(S)) == 1
                    % Fallback: file exists but variable inside has a
                    % different name (e.g. legacy 'conditions' vs 'condition').
                    fn = fieldnames(S);
                    payload = S.(fn{1});
                else
                    warning('NGL03:varNotInFile', ...
                        '%s does not contain variable ''%s''; skipping.', ...
                        fpath, varName);
                    continue
                end
                perSubject.(['all' varName]){1, y} = payload;
            end
        end

        % Save the per-subject aggregated file.
        outDir = fullfile(input.analysis, subject);
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        outFile = fullfile(outDir, [subject '_aggregated.mat']);
        save(outFile, '-struct', 'perSubject', '-v7.3');
        fprintf('  %s\n', outFile);

        clear perSubject
    end
end

%% 05. Cross-subject aggregation.
% Stitches the per-subject files (just produced in §04) into a single
% study-level (Nsubj x maxSess) cell array per variable. We re-read from
% disk rather than carry the per-subject structs in memory, so this also
% works if §04 ran in a previous MATLAB session and the workspace is fresh.
if opt.aggregateSubjects
    fprintf('NGL03_acrossSession: aggregating subjects...\n');

    maxSess = max(arrayfun(@(s) s.nsessions, input.sessions));

    studyLevel = struct();
    for v = 1:size(varList,1)
        studyLevel.(['all' varList{v,1}]) = cell(input.nsubjects, maxSess);
    end

    for x = 1:input.nsubjects
        subject        = input.subjects(x).name;
        perSubjectFile = fullfile(input.analysis, subject, [subject '_aggregated.mat']);
        if ~isfile(perSubjectFile)
            warning('NGL03:missingPerSubject', ...
                ['Per-subject aggregated file missing for %s: %s. Run with ', ...
                 'opt.aggregateSessions=true first, or process the missing ', ...
                 'subject. Leaving that row empty.'], ...
                subject, perSubjectFile);
            continue
        end
        S = load(perSubjectFile);
        for v = 1:size(varList,1)
            fld = ['all' varList{v,1}];
            if isfield(S, fld)
                rowLen = size(S.(fld), 2);
                rowLen = min(rowLen, maxSess);
                studyLevel.(fld)(x, 1:rowLen) = S.(fld)(1, 1:rowLen);
            end
        end
    end

    outFile = fullfile(input.analysis, 'aggregated.mat');
    save(outFile, '-struct', 'studyLevel', '-v7.3');
    fprintf('  %s\n', outFile);
end

% ----------------------------------------------------------------------
function folder = localSessionFolder(input, sourceKey, subject, session)
% Resolve <input.<sourceKey>>/<subject>/<session>.
    switch sourceKey
        case 'analysis',    folder = fullfile(input.analysis,    subject, session);
        case 'spikeSorted', folder = fullfile(input.spikeSorted, subject, session);
        case 'trialSorted', folder = fullfile(input.trialSorted, subject, session);
        otherwise
            error('NGL03:badSourceKey', ...
                'Unknown source key ''%s'' in varList. Use analysis|spikeSorted|trialSorted.', ...
                sourceKey);
    end
end
