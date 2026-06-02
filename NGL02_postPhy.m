%% NGL02_postPhy
% To run after manual curation of desired sessions is completed.
% Reads KS results after manual curation, builds spike and LFP variables
% in lab-standard format, and performs trial-sorting, firing-rate
% calculation, and optional LFP time-frequency analysis.
%
% Requires: NGL_SetAndRunMe.m has been run, and Phy curation is complete.
%
% Jesus 02.06.2026

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

            if needsBuild
                if isMultiArea
                    spike = struct();
                    for a = 1:numel(areaList)
                        areaName         = areaList{a};
                        optArea          = opt;
                        optArea.area     = areaName;
                        optArea.KSfolder = opt.KSfolders.(areaName);
                        fprintf('NGL02 doSpikething: loading clusters for area %s\n', areaName);
                        tmp = loadSpikes(optArea);
                        if isfield(tmp,'spike'), tmp = tmp.spike; end
                        spike.(areaName) = tmp;
                    end
                else
                    % Single-area: opt.area defaults to 'all' from set_default.
                    spike = loadSpikes(opt);
                    if isfield(spike,'spike'), spike = spike.spike; end
                end
                save(fullfile(opt.spikeSorted, "spike.mat"), 'spike', '-mat');
            end

            %% neurons: sort spikes into trials, per area if multi
            if ~exist(fullfile(opt.analysis, "neurons.mat"),'file')
                if ~exist('trialdef','var'), load(fullfile(opt.trialSorted, "trialdef.mat")); end

                if isMultiArea
                    neurons = struct();
                    for a = 1:numel(areaList)
                        areaName = areaList{a};
                        [neurons.(areaName), ~] = sort2trials(spike.(areaName), trialdef, opt);
                    end
                else
                    [neurons, ~] = sort2trials(spike, trialdef, opt);
                end
                save(fullfile(opt.analysis, "neurons.mat"), 'neurons', '-mat')
            end

            %% fireRate: per area if multi
            if ~exist(fullfile(opt.analysis, "fireRate.mat"),'file')
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
                save(fullfile(opt.analysis, "fireRate.mat"), 'fireRate', '-mat')

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
                        [neurons.interactions, ~] = sort2trials(spike, blob.Merges, opt);
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

