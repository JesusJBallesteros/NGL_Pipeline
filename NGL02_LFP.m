%% NGL02_LFP. LFP quick-look.
% PURPOSE:
%   Per-session LFP quick-look processing. Loads the FieldTrip-formatted
%   LFP file produced by NGL01, optionally rejects artifacts, and runs
%   the minimal time-frequency analyses needed for immediate QC:
%       * continuous mode (opt.spectrogram + ~opt.trialparsed):
%           computeContinuousTFR -> saves TFR .mat + provenance
%           plotContinuousTFR    -> writes the summary PNG
%       * trial-parsed mode (opt.spectrogram + opt.trialparsed):
%           trialparsed_MTspectrogram, called per opt.alignto with the
%           alignment tag threaded into the save filename.
%   Independent of Phy curation - can be run any time after NGL01
%   finishes.
%
% SCOPE (26.06.2026):
%   This stage is intentionally kept minimal. It's the "does the LFP
%   look sane on this session" pass. Research-grade LFP analyses (trial-
%   parsed TFR per condition, oscillation / burst detection, phase
%   extraction, spike-field coupling / PPC, LFP-behavior regression)
%   are the domain of the PLANNED NGL07_LFPanalysis, which runs AFTER
%   NGL02_postPhy AND NGL06_videoAnalysis so it can consume spikes and
%   behavioral covariates.
%
% USAGE:
%   Do NOT run or edit this script directly. Configure and run from your
%   project's copy of NGL_SetAndRunMe.m. NGL_SetAndRunMe defines
%   'datadrive', 'studyname', 'subjects', 'dates', and 'opt' before
%   calling this script.
%
% REQUIRED WORKSPACE VARIABLES (set in NGL_SetAndRunMe):
%   datadrive, studyname, subjects, dates, opt - same as NGL02_postPhy.
%
% PIPELINE:
%   00.  NGL00_Prep            - parse inputs into input/opt structs
%   00b. Areas recovery        - pull input.Areas from master preprocInfo
%   01.  set_default           - validate opt, build paths
%   02.  findSessions          - discover session folders on disk
%     per session:
%   03.  prepSession           - prepforsession + pre-flight + applyPreprocInfo
%   04.  Load FT_data          - continuous or trial-parsed
%   05.  Artifact rejection    - if opt.artifdet (once, upstream)
%   06.  Time-frequency        - if opt.spectrogram (compute + plot)
%
% OUTPUTS (per session, paths set in prepforsession):
%   <SavFileName>_FT_data_NoArtif.mat            - artifact-rejected FT_data
%                                                  (if opt.artifdet)
%   <SavFileName>_TFR_continuous.mat             - continuous-mode TFR
%                                                  (if opt.spectrogram && ~opt.trialparsed)
%   <SavFileName>_<align>_TFR.mat                - trial-parsed TFR, one per
%                                                  alignment in opt.alignto
%                                                  (if opt.spectrogram && opt.trialparsed)
%   <plots>/TFR/Cont_allCh.png                   - continuous TFR heatmap
%   plus trial-parsed plots when opt.chbych / opt.trialbytrial gates fire.
%
%   Cross-session TFR aggregation (allTFR_continuous, allTFR_trialparsed)
%   is NOT done in this script. That responsibility is deferred to
%   NGL03_aggregate (see the LFP-aggregation TODO in that file's header).
%
% DEPENDENCIES:
%   set_default, findSessions, prepSession, loadPreprocInfo,
%   applyPreprocInfo, artifact_detRej_lfp,
%   computeContinuousTFR, plotContinuousTFR,
%   trialparsed_MTspectrogram.
%   All toolbox paths are added automatically by set_default.
%
% SEE ALSO:
%   NGL07_LFPanalysis  (planned; research-grade session LFP stage).
%   functions/_deprecated/LFP_Fieldtrip.m         (dead scaffolding).
%   functions/_deprecated/continous_MTspectrogram.m (superseded by
%       computeContinuousTFR + plotContinuousTFR pair).
%
% Last modified 26.06.2026 (Jesus)

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

        if opt.doLFPthing
            % Common per-session loads (used by both modes)
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

                % old FT files have no chanArea tag; ensureChanArea attaches it
                %  from input.areaMap or defaults to 'main'.
                areaMapForBackfill = [];
                if isfield(input, 'areaMap'), areaMapForBackfill = input.areaMap; end
                FT_data = ensureChanArea(FT_data, areaMapForBackfill);

                % 05a. Artifact rejection (continuous).
                if opt.artifdet
                    FT_data = artifact_detRej_lfp(FT_data, opt);
                    save(fullfile(opt.FolderProcDataMat, ...
                                  [opt.SavFileName '_FT_data_NoArtif.mat']), ...
                         'FT_data', '-mat');
                end

                % 06a. Time-frequency (continuous).
                if opt.spectrogram
                    TFR = computeContinuousTFR(FT_data, opt);
                    plotContinuousTFR(TFR, param, opt);
                    save(fullfile(opt.analysis, ...
                                  [opt.SavFileName '_TFR_continuous.mat']), ...
                         'TFR', 'param', 'opt', '-mat');
                end

            else
                %% 04b. TRIAL-PARSED mode: one FT file per opt.alignto entry.
                % Each alignment is loaded and processed independently;
                assert(isfield(opt,'alignto') && iscell(opt.alignto) && ~isempty(opt.alignto), ...
                    'NGL02_LFP:badAlignto', ...
                    'opt.alignto must be a non-empty cell array of alignment names.');

                % opt.lfp.alignSubset (Pass 3 follow-up): let LFP runs
                % restrict to a subset of the study-wide opt.alignto
                % without changing the spike side. Empty -> use all.
                alignsToRun = opt.alignto;
                if isfield(opt,'lfp') && isfield(opt.lfp,'alignSubset') && ~isempty(opt.lfp.alignSubset)
                    keep = ismember(opt.alignto, opt.lfp.alignSubset);
                    if ~any(keep)
                        warning('NGL02_LFP:emptyAlignSubset', ...
                            ['opt.lfp.alignSubset = {%s} did not intersect opt.alignto = {%s}; ', ...
                             'nothing to run for this session.'], ...
                            strjoin(opt.lfp.alignSubset, ', '), strjoin(opt.alignto, ', '));
                        continue
                    end
                    alignsToRun = opt.alignto(keep);
                    fprintf('NGL02_LFP: opt.lfp.alignSubset restricts LFP to {%s}.\n', ...
                            strjoin(alignsToRun, ', '));
                end

                for k = 1:numel(alignsToRun)
                    alignName = alignsToRun{k};
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

                    % attach chanArea if missing
                    areaMapForBackfill = [];
                    if isfield(input, 'areaMap'), areaMapForBackfill = input.areaMap; end
                    FT_data = ensureChanArea(FT_data, areaMapForBackfill);

                    %% 05b. Artifact rejection (per alignment).
                    if opt.artifdet
                        FT_data = artifact_detRej_lfp(FT_data, opt);
                        save(fullfile(opt.trialSorted, ...
                                      [opt.SavFileName '_' alignName '_FT_data_NoArtif.mat']), ...
                             'FT_data', '-mat');
                    end

                    %% 06b. Time-frequency (per alignment).
                    if opt.spectrogram
                        % SocialLearning-specific 
                        if opt.proj_socialLearning
                            if ~isfield(param,'testname')
                                param.testname = 'ASL_Clean_Final_Correct';
                            end
                            opt.analysisCode = input.analysisCode; % required by ASL TFR path
                        else
                            if ~isfield(param,'testname'), param.testname = 'trial_TFR_'; end
                        end

                        % the GENERIC path, loop per area so the
                        % quick-look plotter has area-scoped input.
                        if opt.proj_socialLearning || ~isfield(input,'areaMap') || isempty(input.areaMap)
                            areasToRun = {''};   % '' -> no filter
                        else
                            areasToRun = input.areaMap.uniqueAreas(:).';
                        end

                        % Precompute which areas are actually present in
                        % THIS session's FT_data. input.areaMap is study-
                        % wide, but individual sessions can carry only a
                        % subset of areas (e.g. a single-headstage day in
                        % a multi-area study). Trying to filter to an
                        % absent area would hit NGL:computeTrialparsedTFR:noChans
                        % and kill the run; skip with a warning instead.
                        if opt.proj_socialLearning || isempty(areasToRun) || isempty(areasToRun{1})
                            presentInSession = containers.Map();
                        else
                            presentInSession = containers.Map( ...
                                unique(FT_data.chanArea), ...
                                num2cell(true(1, numel(unique(FT_data.chanArea)))));
                        end

                        for aIdx = 1:numel(areasToRun)
                            areaTag = areasToRun{aIdx};
                            if ~opt.proj_socialLearning
                                if isempty(areaTag)
                                    opt.lfp.tfrAreaFilter = '';
                                else
                                    if ~isKey(presentInSession, areaTag)
                                        warning('NGL02_LFP:areaAbsent', ...
                                            ['Area ''%s'' from input.areaMap has 0 channels in this session ', ...
                                             '(FT_data.chanArea unique = {%s}); skipping this area for [%s/%s/%s].'], ...
                                            areaTag, ...
                                            strjoin(unique(FT_data.chanArea), ', '), ...
                                            input.subjects(input.run(1)).name, ...
                                            input.sessions(input.run(1)).list{input.run(2)}, ...
                                            alignName);
                                        continue
                                    end
                                    opt.lfp.tfrAreaFilter = areaTag;
                                end
                            end
                            [TFR, TFRcfg] = trialparsed_MTspectrogram( ...
                                FT_data, condition, param, opt, alignName);

                            % Filename: append area on the generic path
                            % so per-area outputs don't clobber each other.
                            if opt.proj_socialLearning || isempty(areaTag)
                                outStem = [opt.SavFileName '_' alignName];
                            else
                                outStem = [opt.SavFileName '_' alignName '_' areaTag];
                            end
                            save(fullfile(opt.analysis, [outStem '_TFR.mat']), ...
                                    'TFR', 'TFRcfg', 'param', 'opt', '-mat');

                            % Quick-look condition-mean spectrogram
                            if ~opt.proj_socialLearning
                                try
                                    % Plot spec pulls trialFilter + baseline
                                    % from opt.lfp.plot.* (schema-defaulted).
                                    plotTrialFilter = 'correct';
                                    plotBaseline    = [];
                                    if isfield(opt,'lfp') && isfield(opt.lfp,'plot')
                                        if isfield(opt.lfp.plot,'trialFilter')
                                            plotTrialFilter = opt.lfp.plot.trialFilter;
                                        end
                                        if isfield(opt.lfp.plot,'baseline')
                                            plotBaseline = opt.lfp.plot.baseline;
                                        end
                                    end
                                    plotSpec = struct( ...
                                        'alignName',   alignName, ...
                                        'areaTag',     areaTag, ...
                                        'trialFilter', plotTrialFilter, ...
                                        'baseline',    plotBaseline);
                                    plotTrialparsedTFR_example(TFR, condition, opt, plotSpec);
                                catch ME
                                    warning('NGL02_LFP:quickPlot', ...
                                        'plotTrialparsedTFR_example failed for align=%s, area=%s: %s', ...
                                        alignName, areaTag, ME.message);
                                end
                            end

                            clear TFR TFRcfg
                        end
                    end

                    % Drop this alignment's FT_data before the next iteration
                    % so we don't accidentally reuse it.
                    clear FT_data
                end
            end
        end

        %% Clean up to move on to next session.
        clear FT_data trialdef condition param
    end
end
