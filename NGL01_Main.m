%% NGL01_Main. NGL Electrophysiology Preprocessing Pipeline (Stage 1)
%
% PURPOSE:
%   Main script. Processes one or more recording sessions from raw INTAN
%   or Deuteron data into spike-sorting (.bin) and LFP (FieldTrip .mat)
%   formats, extracts event codes and trial definitions, runs Kilosort 4
%   spike sorting, and optionally runs Bombcell and opens Phy for manual
%   curation.
%
% USAGE:
%   Do NOT run or edit this script directly. Instead, configure and run from
%   your project's copy of NGL_SetAndRunMe.m (stored in analysisCode\).
%   NGL_SetAndRunMe defines 'datadrive', 'studyname', 'subjects', 'dates',
%   and 'opt' before calling this script.
%
% REQUIRED WORKSPACE VARIABLES (set at NGL_SetAndRunMe):
%   datadrive   - char, drive letter (e.g. 'D')
%   studyname   - char, project folder name
%   subjects    - char 'all' or cell array of subject IDs (e.g. {'ABC'})
%   dates       - char 'all' or cell array of session dates (e.g. {'20260101'})
%   opt         - struct, user options (merged with non-explicit defaults)
%
% PIPELINE:
%   00. set_default       - validate opt, build paths, load dependencies
%   01. findSessions      - discover session folders on disk
%   02. prepforsession    - per-session path setup and format detection
%   03. processing        - INTAN or Deuteron Wrappers
%   04. master_kilosort4  - Kilosort 4 spike sorting (loop per area if multi-area)
%   05. Bombcell_Main     - automatic cluster quality (if opt.bombcell)
%   06. Phy               - manual curation GUI (if opt.phy; blocks MATLAB)
%
% OUTPUTS (per session, paths set in prepforsession):
%   <session>.bin       - flat int16 binary for Kilosort (preprocessing\)
%   kilosort\           - KS4 output folder, or <Area>\ in multi-area mode
%   EventRecord.mat     - raw event list (preprocessing\)
%   trialdef.mat        - trial boundaries in ms (trialSorted\)
%   events.mat          - trial-aligned event struct (trialSorted\)
%   *_FTcont.mat        - continuous FieldTrip LFP (trialSorted\, if FieldTrip)
%   *_<event>.mat       - trial-parsed FieldTrip LFP (trialSorted\, if FieldTrip)
%
% DEPENDENCIES:
%   set_default, findSessions, prepforsession, INTAN_PipelineWrapper,
%   Deuteron_PipelineWrapper, master_kilosort4, Bombcell_Main
%   All toolbox paths are added automatically by set_default.
%
% IKN folder standard:
%   gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure
%
% Last modified 18.06.2026 (Jesus) - regen-from-preprocessed mode:
%                                     opt.regenFrom.preproc=true skips all
%                                     raw-side steps (Kilosort, Bombcell,
%                                     NWB, format wrappers) and only
%                                     re-runs EventProcess so events /
%                                     trialdef / condition are rebuilt
%                                     from existing EventRecord.mat.

%% 00. Check current inputs.
NGL00_Prep
[input, opt] = set_default(input, opt);

%% 00b. Regen-from-preprocessed mode (LOCAL-PC PostPhy).
% When opt.regenFrom.preproc=true the user is re-running NGL01 against a
% data tree where the raw folder is intentionally empty (typically the
% curated outputs were moved to a different machine and the raw .dat /
% .rhd files were not transferred). Force the heavy raw-side stages
% OFF so we never try to read what isn't there, and rebuild
% events/trialdef/conditions from the existing EventRecord.mat instead.
regenMode = isfield(opt,'regenFrom') && isfield(opt.regenFrom,'preproc') ...
            && opt.regenFrom.preproc;
if regenMode
    fprintf('\nNGL01_Main: regenFrom.preproc=true (system=''%s''). Forcing kilosort / bombcell / phy / doNWB OFF; will rebuild events / trialdef / conditions only.\n\n', ...
            opt.regenFrom.system);
    opt.kilosort  = false;
    opt.bombcell  = false;
    opt.callBcGUI = false;
    opt.phy       = false;
    opt.doNWB     = false;
end

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

%% 01b. Save study-wide preprocessing snapshot (master copy in analysisCode\).
% Captures the resolved, post-set_default opt before any per-session paths
% are added. NGL02_postPhy can fall back to this if a per-session snapshot
% is missing. Overwritten on every NGL01 run.
savePreprocInfo(input, opt, 'master');

for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        %% 02. Prepare to proceed with a single session.
        input.run = [x y];
        [input, opt] = prepforsession(input, opt);

        %% 03. Proceed to appropriate pipeline.
        if regenMode
           % 03.0 Regen: raw is absent. Skip the format wrapper.
           % Loads an existing EventRecord.mat, re-runs and writes 
           % refreshed events.mat, trialdef.mat, condition.mat
           EventProcess(input, opt);
           savePreprocInfo(input, opt, 'session');
           clear FT_data INTANdata txt
           continue
        end

        switch input.sessions(input.run(1)).info.fileformat
            case {'DT2', 'DF1'}
               % 03.1 Deuteron Pipeline
               [input, opt] = Deuteron_PipelineWrapper(input, opt);

            case {'fileperch', 'filepertype', 'tradFormat'}
               % 03.2 INTAN Pipeline
               [input, opt] = INTAN_PipelineWrapper(input, opt);

            case {'FieldTrip'}
               % 03.3 FT Pipeline — session already preprocessed, events only.
               EventProcess(input, opt);
               disp('Session skipped because continuous FT file was found')
               continue

            otherwise
               warning('Something went wrong during format verification. Skipping Session');
               continue
        end

        %% 04. Kilosort
        if opt.kilosort
            if isfield(input, 'areaMap') && ~isempty(input.areaMap)
                % Multi-area: one KS run per unique area.
                % The shared .bin file is reused across runs; only the
                % connected mask (via chanMap) and output folder differ.
                for a = 1:numel(input.areaMap.uniqueAreas)
                    areaName              = input.areaMap.uniqueAreas{a};
                    fprintf('\n Kilosort: area %s \n', areaName);
                    optArea               = opt;
                    optArea.KSfolder      = opt.KSfolders.(areaName);
                    optArea.KSchanMapFile = input.areaMap.chanMapFiles{a};
                    master_kilosort4(input, optArea)
                end
            else
                % Single-area (default): existing behaviour unchanged.
                master_kilosort4(input, opt)
            end
            close all
        end

        %% 05. Bombcell
        if opt.bombcell
            if isfield(input, 'areaMap') && ~isempty(input.areaMap)
                % Multi-area: run Bombcell for each area's KS output folder.
                for a = 1:numel(input.areaMap.uniqueAreas)
                    areaName         = input.areaMap.uniqueAreas{a};
                    fprintf('\n Bombcell: area %s \n', areaName);
                    optArea          = opt;
                    optArea.KSfolder = opt.KSfolders.(areaName);
                    Bombcell_Main(input, optArea)
                end
            else
                % Single-area (default): existing behaviour unchanged.
                Bombcell_Main(input, opt)
            end
        end

        %% 06. Open Phy to manual curation or just inspection
        % This is BEST done manually, once all your sessions have been processed
        % by opening Phy one by one. This automatization after a session is
        % processed could be useful in specific cases at the time of parameter
        % optimization, or checking specific datasets one by one.
        if opt.phy
            % Will change to current session directory and open phy.
            % ! Keeps MATLAB busy until interface is closed.
            cd(opt.FolderProcDataMat)
            system('phy template-gui params.py');
        end

        %% 07. Save per-session preprocessing snapshot.
        % Authoritative record of the exact opt used for THIS session,
        % including the resolved session-specific paths. Lives next to
        % trialdef.mat / events.mat so NGL02_postPhy can pick it up
        % naturally at the top of each iteration.
        savePreprocInfo(input, opt, 'session');

        %% Clean up to move on to next session
        clear FT_data INTANdata txt

    end % sessions loop
end % subjects loop
