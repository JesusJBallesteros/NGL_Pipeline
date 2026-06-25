%% NGL02_postPhy. NGL Electrophysiology Spike Analysis (Stage 2 - spike path)
%
% PURPOSE:
%   Per-session SPIKE processing after Phy curation. Loads curated KS+Phy
%   clusters into the lab-standard spike struct, sorts spikes into trials,
%   computes firing rate, optionally runs population-dynamics analysis.
%   Multi-area aware: when input.Areas is set, builds nested
%   spike.<Area> / neurons.<Area> / fireRate.<Area> / neuralDynamics.<Area>.
%
% USAGE:
%   Do NOT run or edit this script directly. Configure and run from your
%   project's copy of NGL_SetAndRunMe.m. NGL_SetAndRunMe defines
%   'datadrive', 'studyname', 'subjects', 'dates', and 'opt' before
%   calling this script.
%
% REQUIRED WORKSPACE VARIABLES (set in NGL_SetAndRunMe):
%   datadrive, studyname, subjects, dates, opt - same as NGL02_LFP.
%
% REQUIRES on disk (per session, produced by NGL01_Main + Phy curation):
%   - Kilosort/Phy output: <preproc>/kilosort4/ (single-area) or
%                          <preproc>/<Area>/   (multi-area), each with
%                          params.py and cluster_info.tsv.
%   - <trialSorted>/trialdef.mat, events.mat, condition.mat.
%   The pre-flight checkNGL01Outputs (inside prepSession) errors up
%   front if any of these are missing for the current session.
%
% PIPELINE:
%   00.  NGL00_Prep            - parse inputs into input/opt structs
%   00b. Areas recovery        - pull input.Areas from master preprocInfo
%   01.  set_default           - validate opt, build paths
%   02.  findSessions          - discover session folders on disk
%     per session:
%   2.01. prepSession          - prepforsession + pre-flight + preprocInfo overlay
%   2.02. Offline video        - if opt.offlineTrack (project-specific)
%   2.03. SPIKE branch         - if opt.doSpikething
%           loadSpikes -> spike
%           sort2trials -> neurons
%           calculate_fireRate_general -> fireRate
%           calculate_population_dynamics -> neuralDynamics (if opt.popDyn.do)
%           social-tracking sub-step (if opt.useTrack, single-area only)
%   2.04. LFP branch warning   - if opt.doLFPthing reaches here, point
%                                user at NGL02_LFP.m sibling script
%
% OUTPUTS (per session, paths set by prepforsession):
%   <spikeSorted>/spike.mat           Curated clusters (flat or nested per area)
%   <analysis>/neurons.mat            Per-trial spike times per cluster
%   <analysis>/fireRate.mat           Binned firing rates
%   <analysis>/neuralDynamics.mat     (if opt.popDyn.do)
%   <analysis>/blob.mat               (if opt.offlineTrack)
%   <analysis>/plots/single_fr/       Per-cluster firing-rate heatmaps
%   <analysis>/plots/population_dynamics/  PCA trajectories (if popDyn.pca)
%
% DEPENDENCIES:
%   set_default, findSessions, prepSession, loadPreprocInfo, applyPreprocInfo,
%   checkNGL01Outputs, loadSpikes, sort2trials, calculate_fireRate_general,
%   calculate_population_dynamics (and the population-dynamics family),
%   processAndTrack_video, getSocialEvents.
%
% SEE ALSO:
%   NGL02_LFP (LFP path, runs independently of Phy),
%   NGL03_acrossSession (cross-session aggregator).
%
% Last modified 02.06.2026 (Jesus) - NGL01-style doc header (#9 K)

%% 00. Check current inputs.
% Check if input variable exist already. Parse values.
NGL00_Prep

%% 00b. Pre-flight: recover input.Areas from the NGL01 master snapshot.
% set_default needs input. Areas BEFORE it can build the area map. If the
% user did not re-set it in NGL_SetAndRunMe (common case when NGL02 runs in
% a fresh MATLAB session), pull it from analysisCode\preprocInfo_lastRun.mat.
% A user-supplied input.Areas always overrides.
if ~isfield(input,'Areas') || isempty(input.Areas)
    analysisCodePath = fullfile(input.datadrive, input.studyName, 'analysisCode');
    if ~contains(analysisCodePath, ':\') && ~isempty(input.datadrive)
        % datadrive might be a bare letter at this point; normalise.
        analysisCodePath = fullfile([input.datadrive(1) ':\'], input.studyName, 'analysisCode');
    end
    try
        masterInfo = loadPreprocInfo(analysisCodePath, 'master');
        if isfield(masterInfo,'Areas') && ~isempty(masterInfo.Areas)
            input.Areas = masterInfo.Areas;
            fprintf('NGL02: recovered input.Areas = {%s} from preprocInfo_lastRun.mat\n', ...
                    strjoin(input.Areas, ', '));
        end
    catch ME
        if strcmp(ME.identifier, 'NGL:loadPreprocInfo:notFound')
            % No master snapshot. Continue in single-area mode.
        else
            warning('NGL02:preflight', 'Could not read master preprocInfo: %s', ME.message);
        end
    end
end

% This call guarantees opt is complete, validated, and consistent.
% It will error early with a clear message if anything is wrong.
[input, opt] = set_default(input, opt);

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

%% 02. Proceed with data per session
for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        %% 2.01. Per-session scaffolding (prepforsession + pre-flight + preprocInfo overlay).
        input.run = [x y];
        [input, opt] = prepSession(input, opt);

        %% 2.02. Offline Video blob detector
        % Very specific for Social learning videos from central cenital camera. 
        if opt.offlineTrack
            [blob] = processAndTrack_video(opt);
            save(fullfile(opt.analysis, "blob.mat"), 'blob');
        end

        %% 2.03. SPIKE DATA
        if opt.doSpikething
            % Multi-area ready. When input.areaMap is set (fr. NGL01),
            % build nested spike.<Area>, neurons.<Area>, fireRate.<Area>, 
            % etc.<Area>.
            % Single-area runs keep the flat struct.
            isMultiArea = isfield(input, 'areaMap') && ~isempty(input.areaMap);
            if isMultiArea
                areaList = input.areaMap.uniqueAreas;
            else
                areaList = {'all'};  % only for log; opt.area stays 'all'
            end

            % param is left undefined by user setups that have moved
            % away from postPhy_param.m; guarantee an empty struct so
            % downstream functions can fall back to their inline defaults.
            if ~exist('param','var'), param = struct(); end

            %% spike: load cached file, or build from KS/Phy per area
            needsBuild = ~exist(fullfile(opt.spikeSorted, "spike.mat"), 'file');
            if ~needsBuild
                load(fullfile(opt.spikeSorted, "spike.mat"));  % loads 'spike'
                if isMultiArea
                    for a = 1:numel(areaList)
                        if ~isfield(spike, areaList{a})
                            warning('NGL02:spikeMatShapeMismatch', ...
                                ['Saved spike.mat is missing area "%s" (likely ', ...
                                 'produced by an older single-area run). Regenerating.'], ...
                                areaList{a});
                            needsBuild = true;
                            break
                        end
                    end
                end
            end

            % skipAreas accumulates names of areas that have no usable
            % clusters (loadSpikes errored or returned 0). Every per-area
            % loop below checks this list and skips those areas. Single-
            % area equivalent: skipSession = true means we bail out of
            % the whole doSpikething block for this session.
            skipAreas    = {};
            skipSession  = false;

            if needsBuild
                if isMultiArea
                    spike = struct();
                    for a = 1:numel(areaList)
                        areaName         = areaList{a};
                        optArea          = opt;
                        optArea.area     = areaName;
                        optArea.KSfolder = opt.KSfolders.(areaName);
                        fprintf('NGL02 doSpikething: loading clusters for area %s\n', areaName);
                        try
                            tmp = loadSpikes(optArea);
                            if isfield(tmp,'spike'), tmp = tmp.spike; end
                        catch ME
                            warning('NGL02:loadSpikesFailed', ...
                                'loadSpikes failed for area %s (%s). Marking area as skipped.', ...
                                areaName, ME.message);
                            noteSkippedArea(opt, areaName, ...
                                sprintf(['loadSpikes errored. Most likely cause: no ', ...
                                         'cluster_info.tsv (Phy curation not run yet) ', ...
                                         'or Kilosort produced no output for this area.\n', ...
                                         '\nMATLAB error:\n  %s'], ME.message));
                            skipAreas{end+1} = areaName;
                            continue
                        end
                        if ~isfield(tmp,'label') || isempty(tmp.label)
                            warning('NGL02:emptyArea', ...
                                'Area %s: 0 curated clusters. Marking area as skipped.', areaName);
                            noteSkippedArea(opt, areaName, ...
                                ['loadSpikes returned 0 clusters for this area. ', ...
                                 'Likely all clusters are Phy-labelled noise, or ', ...
                                 'opt.spparams.excludeNoise=true is filtering them out.']);
                            skipAreas{end+1} = areaName;
                            continue
                        end
                        if ~hasCuratedClusters(tmp)
                            warning('NGL02:noCuratedClusters', ...
                                ['Area %s: %d clusters loaded but none have HumanLabel ', ...
                                 'in {good, mua}. Marking area as skipped.'], ...
                                areaName, numel(tmp.label));
                            noteSkippedArea(opt, areaName, ...
                                ['Area has loaded clusters but none have Phy HumanLabel = ', ...
                                 '''good'' or ''mua''. Most likely Phy curation is incomplete: ', ...
                                 'open Phy and assign good / mua / noise to every cluster, then re-run NGL02.']);
                            skipAreas{end+1} = areaName;
                            continue
                        end
                        spike.(areaName) = tmp;
                    end
                    if numel(skipAreas) == numel(areaList)
                        warning('NGL02:allAreasEmpty', ...
                            'All areas have 0 clusters. Skipping doSpikething for this session.');
                        noteSkippedArea(opt, 'all', ...
                            'Every area in input.areaMap returned 0 clusters or errored in loadSpikes; see per-area *_skipped.txt files.');
                        skipSession = true;
                    end
                else
                    % Single-area: opt.area defaults to 'all' from set_default.
                    try
                        spike = loadSpikes(opt);
                        if isfield(spike,'spike'), spike = spike.spike; end
                    catch ME
                        warning('NGL02:loadSpikesFailed', ...
                            'loadSpikes failed (%s). Skipping doSpikething for this session.', ...
                            ME.message);
                        noteSkippedArea(opt, 'all', ...
                            sprintf(['loadSpikes errored. Most likely cause: no ', ...
                                     'cluster_info.tsv (Phy curation not run yet) ', ...
                                     'or Kilosort produced no output.\n', ...
                                     '\nMATLAB error:\n  %s'], ME.message));
                        skipSession = true;
                    end
                    if ~skipSession && (~isfield(spike,'label') || isempty(spike.label))
                        warning('NGL02:emptySession', ...
                            '0 curated clusters. Skipping doSpikething for this session.');
                        noteSkippedArea(opt, 'all', ...
                            ['loadSpikes returned 0 clusters. Likely all clusters are ', ...
                             'Phy-labelled noise, or opt.spparams.excludeNoise=true.']);
                        skipSession = true;
                    end
                    if ~skipSession && ~hasCuratedClusters(spike)
                        warning('NGL02:noCuratedClusters', ...
                            ['%d clusters loaded but none have HumanLabel in {good, mua}. ', ...
                             'Skipping doSpikething for this session.'], numel(spike.label));
                        noteSkippedArea(opt, 'all', ...
                            ['Session has loaded clusters but none have Phy HumanLabel = ', ...
                             '''good'' or ''mua''. Most likely Phy curation is incomplete: ', ...
                             'open Phy and assign good / mua / noise to every cluster, then re-run NGL02.']);
                        skipSession = true;
                    end
                end
                if ~skipSession
                    save(fullfile(opt.spikeSorted, "spike.mat"), 'spike', '-mat');
                end
            end
            if skipSession
                clear events trialdef condition blob neurons spike fireRate
                continue   % next session
            end

            %% neurons: sort spikes into trials, per area if multi
            if ~exist(fullfile(opt.analysis, "neurons.mat"),'file')
                if ~exist('trialdef','var'), load(fullfile(opt.trialSorted, "trialdef.mat")); end

                if isMultiArea
                    neurons = struct();
                    for a = 1:numel(areaList)
                        areaName = areaList{a};
                        if ismember(areaName, skipAreas), continue, end
                        [neurons.(areaName), ~] = sort2trials(spike.(areaName), trialdef, opt);
                    end
                else
                    [neurons, ~] = sort2trials(spike, trialdef, opt);
                end
                save(fullfile(opt.analysis, "neurons.mat"), 'neurons', '-mat')
            end

            %% fireRate: per area if multi
            % Auto-detect cached shape mismatches:
            %   #19: rows must equal numel(condition.aborted) (Ntotal).
            %   #26: columns must equal numel(opt.alignto) (Nalign).
            % If either is wrong (pre-#19 cache has Nremaining rows;
            % pre-#26 cache has 1 column even when multiple alignments
            % were requested), delete the stale cache so the regen path
            % runs with the current shape contract.
            fireRateCache = fullfile(opt.analysis, "fireRate.mat");
            if isfile(fireRateCache)
                if ~exist('condition','var')
                    load(fullfile(opt.trialSorted, "condition.mat"));
                    if ~exist('condition','var'), condition = conditions; clear conditions
                    end
                end
                cached = load(fireRateCache);
                cFr = cached.fireRate;
                % Pull the .sps cell array from either flat or nested layout.
                if isstruct(cFr) && isfield(cFr,'sps') && ~isempty(cFr.sps)
                    spsCells = cFr.sps;
                elseif isstruct(cFr)
                    fn = fieldnames(cFr);
                    spsCells = cFr.(fn{1}).sps;
                else
                    spsCells = {};
                end
                if ~isempty(spsCells)
                    sample = spsCells{1};
                    Ncol   = size(spsCells, 2);
                    needsRegen = false;
                    % Row-count check
                    if isfield(condition,'aborted') && ~isempty(sample) && ...
                            size(sample,1) ~= numel(condition.aborted)
                        warning('NGL02:fireRateShapeMismatch', ...
                            ['Cached fireRate.mat has %d rows but condition ', ...
                             'expects %d (pre-#19 cache). Deleting and regenerating.'], ...
                             size(sample,1), numel(condition.aborted));
                        needsRegen = true;
                    end
                    % Column-count check
                    if ~needsRegen && Ncol ~= numel(opt.alignto)
                        warning('NGL02:fireRateAlignMismatch', ...
                            ['Cached fireRate.mat has %d alignment column(s) but ', ...
                             'opt.alignto requests %d (pre-#26 cache). Deleting ', ...
                             'and regenerating.'], Ncol, numel(opt.alignto));
                        needsRegen = true;
                    end
                    if needsRegen, delete(fireRateCache); end
                end
                clear cached cFr spsCells sample Ncol needsRegen
            end

            if ~exist(fireRateCache,'file')
                if ~exist('neurons','var'),   load(fullfile(opt.analysis, "neurons.mat"));    end
                if ~exist('events','var'),    load(fullfile(opt.trialSorted, "events.mat"));   end
                if ~exist('condition','var')
                    load(fullfile(opt.trialSorted, "condition.mat"));
                    if ~exist('condition','var'), condition = conditions; clear conditions
                    end
                end

                % General function, no conditions: 'allInitiated' by default.
                if isMultiArea
                    fireRate = struct();
                    for a = 1:numel(areaList)
                        areaName = areaList{a};
                        if ismember(areaName, skipAreas), continue, end
                        % Non-block FR calculation, for block/treatment sessions,
                        % use 'calculate_fireRate_byBlock' instead.
                        % No 'events' needed here.
                        fireRate.(areaName) = calculate_fireRate_general( ...
                                                neurons.(areaName), [], condition, opt, param);
                    end
                else
                    % Non-block FR calculation, for block/treatment sessions,
                    % use 'calculate_fireRate_byBlock' instead. No 'events'
                    % needed here.
                    fireRate = calculate_fireRate_general(neurons, [], condition, opt, param);
                end
                save(fullfile(opt.analysis, "fireRate.mat"), 'fireRate', '-mat', '-v7.3')

                % Plotting decoupled from calculate_fireRate_general.
                % param.plot gates the helper; in multi-area
                % mode we call it once per area so titles/filenames
                % carry the right area tag via opt.area.
                if ~isfield(param,'plot') || param.plot
                    if isMultiArea
                        for a = 1:numel(areaList)
                            areaName     = areaList{a};
                            if ismember(areaName, skipAreas), continue, end
                            optArea      = opt;
                            optArea.area = areaName;
                            plot_fireRate_session(fireRate.(areaName), condition, param, optArea);
                        end
                    else
                        plot_fireRate_session(fireRate, condition, param, opt);
                    end
                end

                % % Project specific
                % param.IncludeFS = true; % NS and FS
                % param.trial2plot = 'allInitiated'; % for correct trials
                % fireRate = calculate_fireRate_extintion(neurons, events, condition, opt, param);
                % save(fullfile(opt.analysis, "fireRate_extintion.mat"), 'fireRate', '-mat')
            end

            %% population dynamics: per area if multi.
            % Dispatches to whichever methods are enabled under
            % opt.popDyn (pca / jPCA / GPFA / trialEmbed). jPCA and
            % GPFA are placeholders today; see calculate_neural_jpca.m
            % and calculate_neural_gpfa.m for the roadmap.
            if opt.popDyn.do
                if ~exist(fullfile(opt.analysis, "neuralDynamics.mat"),'file')
                    if ~exist('trialdef','var'), load(fullfile(opt.trialSorted, "trialdef.mat")); end
                    if ~exist('events','var'),    load(fullfile(opt.trialSorted, "events.mat"));   end
                    if ~exist('condition','var')
                        load(fullfile(opt.trialSorted, "condition.mat"));
                        if ~exist('condition','var'), condition = conditions; clear conditions
                        end
                    end
                    if ~exist('neurons','var'),   load(fullfile(opt.analysis, "neurons.mat"));    end
                    if ~exist('fireRate','var'),   load(fullfile(opt.analysis, "fireRate.mat"));    end

                    if isMultiArea
                        neuralDynamics = struct();
                        for a = 1:numel(areaList)
                            areaName         = areaList{a};
                            if ismember(areaName, skipAreas), continue, end
                            optArea          = opt;
                            optArea.area     = areaName;
                            neuralDynamics.(areaName) = calculate_population_dynamics( ...
                                neurons.(areaName), fireRate.(areaName), trialdef, condition, optArea);
                        end
                    else
                        neuralDynamics = calculate_population_dynamics(neurons, fireRate, trialdef, condition, opt);
                    end
                    save(fullfile(opt.analysis, "neuralDynamics.mat"), 'neuralDynamics', '-mat')
                end
            end

            %% Tracking in Social Arena (project-specific; single-area only for now)
            if opt.useTrack
                if isMultiArea
                    warning('NGL02:useTrackMultiArea', ...
                        ['opt.useTrack is not multi-area aware. Skipping ', ...
                         'social-tracking spike indexing for this session.']);
                else
                    % Spiking indexing for Social interactions. Checks blob
                    % interaction times (+-5s) and extracts spiking activity
                    % around them. Input needs to be 'blob.Merges', instead
                    % of a regular 'trialdef'.
                    if ~exist('blob','var'), load(fullfile(opt.analysis, "blob.mat")); end
                    if ~isfield(neurons,"interactions")
                        neurons.interactions = sort2trials_blob(spike, blob.Merges, opt);
                    end
                    save(fullfile(opt.analysis, "neurons.mat"), 'neurons', '-mat')

                    % Also, use video assessment excel files to extract the
                    % logical indexing of Social events.
                    if ~exist('events','var'), load(fullfile(opt.analysis, "events.mat")); end
                    if ~isfield(events,"social")
                        [events.social] = getSocialEvents(opt);
                        save(fullfile(opt.analysis, "events.mat"), 'events', '-mat')
                    end
                end
            end

        end
        
        %% Clean up to move on to next session
        % Drop everything loaded or built inside this iteration so it
        % cannot leak into the next session.
        clear events trialdef condition blob neurons spike fireRate
     
    end
end

