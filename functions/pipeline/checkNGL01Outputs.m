function checkNGL01Outputs(input, opt)
% checkNGL01Outputs  Pre-flight check that the NGL01-stage outputs needed
%                    by NGL02 actually exist for the current session.
%
% PURPOSE:
%   When a user runs NGL02_postPhy on a session that NGL01 never finished
%   (or that has not been curated in Phy yet), the failure typically only
%   shows up several steps into the pipeline as a confusing "file not
%   found" deep inside loadSpikes / sort2trials. This helper looks for
%   every NGL01 output that NGL02 will demand for THIS session, given
%   the option flags in opt, and errors up front with a single clear
%   message listing everything that is missing.
%
% USAGE:
%   checkNGL01Outputs(input, opt)
%   Called from NGL02_postPhy right after prepforsession, before the
%   doSpikething / doLFPthing branches do any work.
%
% INPUTS:
%   input  - struct from set_default. Required: input.run = [x y],
%            input.subjects, input.sessions. Optional: input.areaMap
%            (multi-area mode).
%   opt    - resolved options struct, POST-prepforsession (so the
%            session-specific paths opt.FolderProcDataMat,
%            opt.trialSorted, opt.KSfolder / opt.KSfolders are populated).
%            Read flags: opt.doSpikething, opt.doLFPthing,
%            opt.trialparsed, opt.SavFileName.
%
% CHECKS:
%   Always:
%     - opt.FolderProcDataMat  exists (the preprocessing folder).
%     - <opt.trialSorted>/trialdef.mat
%     - <opt.trialSorted>/events.mat
%   If opt.doSpikething:
%     - For each KS output folder (single or per-area):
%       * params.py
%       * cluster_info.tsv  (Phy curation produced this)
%   If opt.doLFPthing:
%     - <opt.FolderProcDataMat>/<SavFileName>_FTcont.mat     (when ~opt.trialparsed)
%       (continuous FieldTrip data is written by MAT2FieldTrip to the
%        preprocessing folder, not to trialSorted.)
%     - <opt.trialSorted>/<SavFileName>_<align>.mat          (when opt.trialparsed)
%       (one file per alignment in opt.alignto; e.g. with
%        opt.alignto = {'itiOn','stimOn1','rwd'} we expect
%        <SavFileName>_itiOn.mat, <SavFileName>_stimOn1.mat,
%        <SavFileName>_rwd.mat. Every missing file is reported.)
%
% ERROR:
%   NGL02:missingNGL01Output, listing every missing path. The error message
%   names the subject and session so users can spot which one's broken at
%   a glance.
%
% NOTES:
%   - condition.mat is not checked because the spike branch is forgiving
%     about it (it will fall back to loading on demand); a missing
%     condition is signalled later, more gracefully, by the conditions
%     handling inside calculate_fireRate_general.
%   - blob.mat is not checked because it is only required when
%     opt.offlineTrack / opt.useTrack are on, and those project-specific
%     paths have their own guards.
%
% SEE ALSO:
%   prepforsession, loadPreprocInfo, NGL02_postPhy.
%
% Last modified 29.05.2026 (Jesus) - audit item N

x       = input.run(1);
y       = input.run(2);
subject = input.subjects(x).name;
session = input.sessions(x).list{y};

missing = {};

%% Always required.
if ~exist(opt.FolderProcDataMat, 'dir')
    missing{end+1} = sprintf('preprocessing folder: %s', opt.FolderProcDataMat);
end
trialdefMat = fullfile(opt.trialSorted, 'trialdef.mat');
if ~isfile(trialdefMat),  missing{end+1} = trialdefMat; end
eventsMat   = fullfile(opt.trialSorted, 'events.mat');
if ~isfile(eventsMat),    missing{end+1} = eventsMat;   end

%% Required if running the spike branch.
if isfield(opt,'doSpikething') && opt.doSpikething
    if isfield(input,'areaMap') && ~isempty(input.areaMap)
        for a = 1:numel(input.areaMap.uniqueAreas)
            areaName = input.areaMap.uniqueAreas{a};
            ksFolder = opt.KSfolders.(areaName);
            missing  = localCheckKSFolder(missing, ksFolder, areaName);
        end
    else
        missing = localCheckKSFolder(missing, opt.KSfolder, '');
    end
end

%% Required if running the LFP branch.
if isfield(opt,'doLFPthing') && opt.doLFPthing
    if isfield(opt,'trialparsed') && opt.trialparsed
        % Trial-parsed: one FT file per alignment in opt.alignto, sitting
        % in opt.trialSorted as <SavFileName>_<alignment>.mat. Check each.
        assert(isfield(opt,'alignto') && iscell(opt.alignto) && ~isempty(opt.alignto), ...
            'NGL02:badAlignto', ...
            ['checkNGL01Outputs needs opt.alignto (a non-empty cell of ', ...
             'alignment names) when opt.trialparsed=true.']);
        nFound = 0;
        for k = 1:numel(opt.alignto)
            ftFile = fullfile(opt.trialSorted, ...
                              [opt.SavFileName '_' opt.alignto{k} '.mat']);
            if isfile(ftFile)
                nFound = nFound + 1;
            else
                missing{end+1} = ftFile; %#ok<AGROW>
            end
        end
        % Sanity report: explicitly note count mismatch in case the
        % alignment names match but extras exist (e.g. user added an
        % alignment in NGL02 that NGL01 never produced).
        if nFound ~= numel(opt.alignto)
            missing{end+1} = sprintf(['expected %d trial-parsed FT files ', ...
                                       'in %s, found %d'], ...
                                       numel(opt.alignto), opt.trialSorted, nFound);
        end
    else
        % Continuous: single _FTcont.mat in the preprocessing folder
        % (MAT2FieldTrip writes it to opt.FolderProcDataMat, NOT trialSorted).
        ftFile = fullfile(opt.FolderProcDataMat, [opt.SavFileName '_FTcont.mat']);
        if ~isfile(ftFile),  missing{end+1} = ftFile; end
    end
end

%% Report.
if ~isempty(missing)
    msg = sprintf(['NGL02 cannot proceed for %s / %s. The following NGL01 ', ...
                   'outputs are missing:\n  - %s\n', ...
                   'Run NGL01_Main for this session and make sure Phy curation ', ...
                   'has been saved (cluster_info.tsv) before re-running NGL02.'], ...
                   subject, session, strjoin(missing, sprintf('\n  - ')));
    error('NGL02:missingNGL01Output', '%s', msg);
end
end

% ----------------------------------------------------------------------
function missing = localCheckKSFolder(missing, ksFolder, areaTag)
% Verify params.py and cluster_info.tsv exist in a Kilosort/Phy folder.
    if isempty(areaTag), tag = ''; else, tag = sprintf(' (area %s)', areaTag); end

    if ~exist(ksFolder, 'dir')
        missing{end+1} = sprintf('Kilosort folder: %s%s', ksFolder, tag);
        return  % skip file-level checks if the folder itself is missing
    end
    paramsPy = fullfile(ksFolder, 'params.py');
    if ~isfile(paramsPy),    missing{end+1} = sprintf('%s%s', paramsPy, tag); end
    clusterInfo = fullfile(ksFolder, 'cluster_info.tsv');
    if ~isfile(clusterInfo), missing{end+1} = sprintf('%s%s', clusterInfo, tag); end
end
