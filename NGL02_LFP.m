%% NGL02_LFP. NGL Electrophysiology LFP Analysis (Stage 2 — LFP path)
%
% PURPOSE:
%   Per-session LFP processing. Loads the FieldTrip-formatted LFP file
%   produced by NGL01, optionally rejects artifacts, and runs time-
%   frequency analysis (continuous multitaper spectrogram or trial-
%   parsed). Independent of Phy curation — can be run any time after
%   NGL01 finishes.
%
% USAGE:
%   Do NOT run or edit this script directly. Configure and run from your
%   project's copy of NGL_SetAndRunMe.m. NGL_SetAndRunMe defines
%   'datadrive', 'studyname', 'subjects', 'dates', and 'opt' before
%   calling this script.
%
% REQUIRED WORKSPACE VARIABLES (set in NGL_SetAndRunMe):
%   datadrive, studyname, subjects, dates, opt — same as NGL02_postPhy.
%
% PIPELINE:
%   00.  NGL00_Prep            - parse inputs into input/opt structs
%   00b. Areas recovery        - pull input.Areas from master preprocInfo
%   01.  set_default           - validate opt, build paths
%   02.  findSessions          - discover session folders on disk
%     per session:
%   03.  prepSession           - prepforsession + pre-flight + applyPreprocInfo
%   04.  Load FT_data          - continuous or trial-parsed
%   05.  Artifact rejection    - if opt.artifdet
%   06.  Time-frequency        - if opt.spectrogram
%
% OUTPUTS (per session, paths set in prepforsession):
%   <SavFileName>_FT_data_NoArtif.mat  - artifact-rejected FT_data
%                                        (if opt.artifdet)
%   <SavFileName>_TFR_*.mat            - per-session TFR (if opt.spectrogram)
%
%   Cross-session TFR aggregation (allTFR_continuous, allTFR_trialparsed)
%   is NOT done in this script. That responsibility now lives in
%   NGL03_acrossSession (planned; task #7).
%
% DEPENDENCIES:
%   set_default, findSessions, prepSession, loadPreprocInfo,
%   applyPreprocInfo, artifact_detRej_lfp, continous_MTspectrogram,
%   trialparsed_MTspectrogram.
%   All toolbox paths are added automatically by set_default.
%
% Last modified 29.05.2026 (Jesus) - split from NGL02_postPhy (#21)

%% 00. Check current inputs.
NGL00_Prep

%% 00b. Pre-flight: recover input.Areas from the NGL01 master snapshot.
% set_default needs input.Areas BEFORE it can build the area map. If the
% user did not re-set it in NGL_SetAndRunMe (common case when NGL02_LFP
% runs in a fresh MATLAB session), pull it from analysisCode\preprocInfo_
% lastRun.mat. A user-supplied input.Areas always overrides.
if ~isfield(input,'Areas') || isempty(input.Areas)
    analysisCodePath = fullfile(input.datadrive, input.studyName, 'analysisCode');
    if ~contains(analysisCodePath, ':\') && ~isempty(input.datadrive)
        analysisCodePath = fullfile([input.datadrive(1) ':\'], input.studyName, 'analysisCode');
    end
    try
        masterInfo = loadPreprocInfo(analysisCodePath, 'master');
        if isfield(masterInfo,'Areas') && ~isempty(masterInfo.Areas)
            input.Areas = masterInfo.Areas;
            fprintf('NGL02_LFP: recovered input.Areas = {%s} from preprocInfo_lastRun.mat\n', ...
                    strjoin(input.Areas, ', '));
        end
    catch ME
        if strcmp(ME.identifier, 'NGL:loadPreprocInfo:notFound')
            % No master snapshot. Continue in single-area mode.
        else
            warning('NGL02_LFP:preflight', 'Could not read master preprocInfo: %s', ME.message);
        end
    end
end

% Validate opt and build the IKN path structure.
[input, opt] = set_default(input, opt);

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

%% 02. Proceed with data per session
for x = 1:input.nsubjects
    for y = 1:input.sessions(x).nsessions
        %% 03. Per-session scaffolding (prepforsession + pre-flight + preprocInfo overlay).
        input.run = [x y];
        [input, opt] = prepSession(input, opt);

        %% 04. LFP DATA. UNDER DEVELOPMENT
        if opt.doLFPthing
            % --- Common per-session loads (used by both modes) ---
            % trialdef.mat must have been produced by NGL01; the pre-flight
            % in checkNGL01Outputs already enforces this, re-check defensively.
            trialdefPath = fullfile(opt.trialSorted, 'trialdef.mat');
            if ~exist('trialdef','var')
                if isfile(trialdefPath)
                    load(trialdefPath);
                else
                    error('NGL02_LFP:missingTrialdef', ...
                        ['trialdef.mat not found at %s. This file is produced ', ...
                         'by NGL01_Main; rerun NGL01 for this session before ', ...
                         'NGL02_LFP.'], trialdefPath);
                end
            end

            % param is left undefined by user setups that have moved away
            % from postPhy_param.m; guarantee an empty struct so downstream
            % functions can fall back to their inline defaults.
            if ~exist('param','var'), param = struct(); end

            % condition is loaded on demand (some functions consume it).
            if ~exist('condition','var') && isfile(fullfile(opt.trialSorted, "condition.mat"))
                load(fullfile(opt.trialSorted, "condition.mat"));
                if ~exist('condition','var') && exist('conditions','var')
                    condition = conditions; clear conditions
                end
            end
            if ~exist('condition','var'), condition = struct(); end

            %  Dispatch on continuous vs trial-parsed
            if ~opt.trialparsed
                %% 04a. CONTINUOUS mode: single FT file in FolderProcDataMat.
                ftFile = fullfile(opt.FolderProcDataMat, ...
                                  [opt.SavFileName '_FTcont.mat']);
                fprintf('NGL02_LFP: loading continuous FT file %s\n', ftFile);
                load(ftFile, "-mat", 'FT_data');
                if isfield(FT_data,"FT_data"), FT_data = FT_data.FT_data; end
                FT_data.cfg.continuous = 'yes';

                %% 05a. Artifact rejection (continuous).
                if opt.artifdet
                    FT_data = artifact_detRej_lfp(FT_data, opt);
                    save(fullfile(opt.FolderProcDataMat, ...
                                  [opt.SavFileName '_FT_data_NoArtif.mat']), ...
                         'FT_data', '-mat');
                end

                %% 06a. Time-frequency (continuous).
                if opt.spectrogram
                    TFR = continous_MTspectrogram(FT_data, condition, param, opt);
                    save(fullfile(opt.analysis, ...
                                  [opt.SavFileName '_TFR_continuous.mat']), ...
                         'TFR', 'param', 'opt', '-mat');
                end

            else
                %% 04b. TRIAL-PARSED mode: one FT file per opt.alignto entry.
                % Each alignment is loaded and processed independently;
                % outputs are tagged by alignment so they don't collide.
                assert(isfield(opt,'alignto') && iscell(opt.alignto) && ~isempty(opt.alignto), ...
                    'NGL02_LFP:badAlignto', ...
                    'opt.alignto must be a non-empty cell array of alignment names.');

                for k = 1:numel(opt.alignto)
                    alignName = opt.alignto{k};
                    ftFile = fullfile(opt.trialSorted, ...
                                      [opt.SavFileName '_' alignName '.mat']);
                    if ~isfile(ftFile)
                        warning('NGL02_LFP:missingTrialparsedFT', ...
                            ['Trial-parsed FT file not found for alignment ', ...
                             '''%s'': %s. Skipping this alignment.'], ...
                            alignName, ftFile);
                        continue
                    end
                    fprintf('NGL02_LFP: loading trial-parsed FT file for alignment %s\n', ...
                            alignName);
                    load(ftFile, "-mat", 'FT_data');
                    if isfield(FT_data,"FT_data"), FT_data = FT_data.FT_data; end

                    %% 05b. Artifact rejection (per alignment).
                    if opt.artifdet
                        FT_data = artifact_detRej_lfp(FT_data, opt);
                        save(fullfile(opt.trialSorted, ...
                                      [opt.SavFileName '_' alignName '_FT_data_NoArtif.mat']), ...
                             'FT_data', '-mat');
                    end

                    %% 06b. Time-frequency (per alignment).
                    if opt.spectrogram
                        % TODO proj_ASL gate (task #8): the testname default
                        % and the input.analysisCode forwarding are
                        % SocialLearning-specific workarounds that will be
                        % moved behind opt.proj_socialLearning.
                        if ~isfield(param,'testname'), param.testname = 'trial_TFR_'; end
                        opt.analysisCode = input.analysisCode;  % TODO clean up

                        [TFR, TFRcfg] = trialparsed_MTspectrogram( ...
                                           FT_data, condition, param, opt);
                        save(fullfile(opt.analysis, ...
                                      [opt.SavFileName '_' alignName '_TFR.mat']), ...
                             'TFR', 'TFRcfg', 'param', 'opt', '-mat');
                    end

                    % Drop this alignment's FT_data before the next iteration
                    % so we don't accidentally reuse it.
                    clear FT_data TFR TFRcfg
                end
            end
        end

        %% Clean up to move on to next session.
        clear FT_data trialdef condition param
    end
end
