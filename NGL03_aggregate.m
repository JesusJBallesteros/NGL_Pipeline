%% NGL03_aggregate. Cross-session and cross-subject aggregator (per-area).
%
% PURPOSE:
%   Stage 3 aggregator: pulls together the per-session outputs that
%   NGL02_postPhy and NGL02_LFP produced and packs them into per-area
%   cell-array files indexed by (subject, session). Two gates, with a
%   dependency:
%       opt.aggregateSessions  -> per-subject per-area files at
%                                 data\analysis\<subject>\
%                                 <subject>_aggregated_<area>.mat
%       opt.aggregateSubjects  -> study-level per-area files at
%                                 data\analysis\aggregated_<area>.mat
%                                 (requires aggregateSessions=true)
%   Both default to false; aggregation is opt-in.
%
% ARCHITECTURE NOTE (refactor 26.06.2026):
%   Pre-26.06.2026 NGL03 wrote ONE nested file per scope:
%       aggregated.mat                 with .allspike{x,y}.NCL / .STR
%       <subj>_aggregated.mat          per-subject equivalent
%   That shape forced every downstream consumer to know about area-nesting
%   and silently broke buildRequestCatSets (which iterated cells as if
%   they were flat spike structs). The new layout writes ONE file PER
%   AREA with a uniformly FLAT shape:
%       aggregated_NCL.mat             allspike{x,y} is a flat spike struct
%       aggregated_STR.mat             ditto, for STR
%       <subj>_aggregated_NCL.mat      per-subject equivalent
%   Single-area runs (input.Areas empty or one entry with no name) write
%   aggregated_main.mat. Consumers (loadAggregatedSpikes / NGL04) load
%   one file per area and never see nesting.
%
% LEGACY:
%   If you have legacy aggregated.mat / <subj>_aggregated.mat from before
%   26.06.2026, run migrate_aggregated_to_perArea.m once to split them
%   into the new per-area form. Re-running NGL03 from the per-session
%   files always produces the new shape directly.
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
%   03b. Resolve areas-to-write list (default 'main' if none declared)
%   04.  Per-subject session aggregation (if opt.aggregateSessions)
%        -> one file per area
%   05.  Cross-subject aggregation       (if opt.aggregateSubjects)
%        -> one file per area
%
% AGGREGATION SHAPE (per area):
%   Cell arrays of size (Nsubj, Nsess) where each cell contains the
%   loaded variable for that (subject, session), with any area-level
%   nesting in the source file DRILLED INTO at load time so the saved
%   payload is the flat single-area shape.
%   - Per-subject file: (1 x nSess_for_this_subject) cells per variable.
%   - Study-level file: (Nsubj x maxSess_across_subjects) cells, padded
%     with empties where a subject ran fewer sessions than maxSess.
%
% VARIABLES AGGREGATED (this version - SPIKE SIDE ONLY):
%   Always (when opt.doSpikething): spike, neurons, fireRate,
%                                    condition, events, trialdef.
%   Conditional:
%     - neuralDynamics  if opt.popDyn.do
%     - blob            if opt.offlineTrack or opt.useTrack
%
%   "Global" payloads (condition / events / trialdef) are not area-
%   nested in the per-session files. They are written verbatim into
%   EVERY per-area aggregated file. Disk overhead is minor and keeps
%   consumers symmetric (they read one file per area and find everything
%   they need there).
%
%   LFP-side aggregation is planned but not yet here. The LFP outputs
%   (continuous, per-alignment trial-parsed, artifact-rejected copies)
%   are multi-file per session and need a separate model. The current
%   spike-side template is the foundation; LFP entries will join the
%   varList in a follow-up.
%
% OUTPUTS:
%   - data\analysis\<subject>\<subject>_aggregated_<area>.mat
%       (one per (subject, area), when aggregateSessions=true)
%   - data\analysis\aggregated_<area>.mat
%       (one per area, when aggregateSubjects=true)
%
% Last modified 26.06.2026 (Jesus)

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

%% 03b. Resolve the list of area tags to write one file per.
% Multi-area: input.Areas is e.g. {'NCL','NCL','STR'}; the unique-stable
% list is {'NCL','STR'} and we write one file per entry.
% Single-area: input.Areas absent or empty -> use 'main' as the sentinel
% so file naming stays uniform (aggregated_main.mat).
if isfield(input,'Areas') && ~isempty(input.Areas)
    areasToWrite = unique(input.Areas, 'stable');
else
    areasToWrite = {'main'};
end
fprintf('NGL03_acrossSession: writing per-area files for: %s\n', ...
        strjoin(areasToWrite, ', '));

%% 04. Session aggregation per subject.
if opt.aggregateSessions
    fprintf('NGL03_acrossSession: aggregating sessions per subject...\n');

    for x = 1:input.nsubjects
        subject = input.subjects(x).name;
        nSess   = input.sessions(x).nsessions;

        % Pre-allocate one (1 x nSess) cell per (area, target).
        perArea = struct();
        for a = 1:numel(areasToWrite)
            ak = areasToWrite{a};
            for v = 1:size(varList,1)
                perArea.(ak).(['all' varList{v,1}]) = cell(1, nSess);
            end
        end

        % Walk sessions, load each existing file into its (area-drilled)
        % per-area accumulator.
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

                % Split per area. Variables that aren't area-nested in
                % the source (e.g. condition / events / trialdef) get
                % written verbatim into every per-area file.
                for a = 1:numel(areasToWrite)
                    ak       = areasToWrite{a};
                    extracted = localExtractAreaPayload(payload, ak, areasToWrite);
                    perArea.(ak).(['all' varName]){1, y} = extracted;
                end
            end
        end

        % Save the per-subject per-area aggregated files.
        outDir = fullfile(input.analysis, subject);
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        for a = 1:numel(areasToWrite)
            ak         = areasToWrite{a};
            outFile    = fullfile(outDir, [subject '_aggregated_' ak '.mat']);
            perSubject = perArea.(ak);
            save(outFile, '-struct', 'perSubject', '-v7.3');
            fprintf('  %s\n', outFile);
        end

        clear perArea perSubject
    end
end

%% 05. Cross-subject aggregation.
% Stitches the per-subject per-area files (just produced in §04) into a
% single per-area (Nsubj x maxSess) cell array per variable. We re-read
% from disk rather than carry the per-subject structs in memory, so this
% also works if §04 ran in a previous MATLAB session and the workspace
% is fresh.
if opt.aggregateSubjects
    fprintf('NGL03_acrossSession: aggregating subjects...\n');

    maxSess = max(arrayfun(@(s) s.nsessions, input.sessions));

    for a = 1:numel(areasToWrite)
        ak = areasToWrite{a};

        studyLevel = struct();
        for v = 1:size(varList,1)
            studyLevel.(['all' varList{v,1}]) = cell(input.nsubjects, maxSess);
        end

        for x = 1:input.nsubjects
            subject        = input.subjects(x).name;
            perSubjectFile = fullfile(input.analysis, subject, ...
                                      [subject '_aggregated_' ak '.mat']);
            if ~isfile(perSubjectFile)
                warning('NGL03:missingPerSubject', ...
                    ['Per-subject aggregated file missing for %s/%s: %s. ', ...
                     'Run with opt.aggregateSessions=true first, or process ', ...
                     'the missing subject. Leaving that row empty.'], ...
                    subject, ak, perSubjectFile);
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

        outFile = fullfile(input.analysis, ['aggregated_' ak '.mat']);
        save(outFile, '-struct', 'studyLevel', '-v7.3');
        fprintf('  %s\n', outFile);
    end
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

function out = localExtractAreaPayload(payload, areaKey, areasToWrite)
% Per-area extractor used at NGL03 write time. Behaviour:
%   * payload empty or non-struct -> pass through verbatim (caller's
%     job to handle nil).
%   * payload is a struct with `areaKey` as a struct-valued field
%     (the multi-area shape, e.g. spike.NCL.HumanLabel) -> drill in
%     and return payload.(areaKey).
%   * payload is a struct WITHOUT `areaKey` but with one of the other
%     areas in areasToWrite -> it IS area-nested but not for this area;
%     return [] (no data for this (subj, sess, area)).
%   * otherwise (flat single-area struct, no area-nesting) -> pass
%     through verbatim. Single-area variables (condition / events /
%     trialdef) take this path and get copied into every area file.
    if isempty(payload), out = []; return; end
    if ~isstruct(payload), out = payload; return; end

    if isfield(payload, areaKey) && isstruct(payload.(areaKey))
        out = payload.(areaKey);
        return
    end

    fns        = fieldnames(payload);
    otherAreas = setdiff(areasToWrite, {areaKey});
    if any(ismember(otherAreas, fns))
        % Nested by area but not THIS area -> empty for this slot.
        out = [];
        return
    end

    % Flat / single-area / non-area-nested payload (e.g. condition).
    out = payload;
end
