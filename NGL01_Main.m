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
%   00. set_default     - validate opt, build paths, load dependencies
%   01. findSessions    - discover session folders on disk
%   02. prepforsession  - per-session path setup and format detection
%   03. processing      - INTAN or Deuteron Wrappers
%   04. master_kilosort4  - Kilosort 4 spike sorting
%   05. Bombcell_Main   - automatic cluster (if opt.bombcell)
%   06. Phy             - manual curation GUI (if opt.phy; blocks MATLAB)
%
% OUTPUTS (per session, paths set in prepforsession):
%   <session>.bin       - flat int16 binary for Kilosort (preprocessing\)
%   kilosort\           - KS4 output folder (preprocessing\)
%   EventRecord.mat     - raw event list (preprocessing\)
%   trialdef.mat        - trial boundaries definitions, in ms (trialSorted\)
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
% Last modified 06.05.2026 (Jesus)

%% 00. Check current inputs.
% Check if input variable exist already. Parse values.
if ~exist("input","var")
    input = struct( 'datadrive' , datadrive , ...   % force char array
                    'studyName' , studyname , ...   % force char array
                    'subjects'  , [], ...           % do NOT force char array
                    'dates'     , []        );      % do NOT force char array
    input.dates     = dates;    % place as it comes
    input.subjects  = subjects; % place as it comes
end

% This single call guarantees opt is complete, validated, and consistent.
% It will break here if anything is wrong.
[input, opt] = set_default(input, opt);

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        %% 02. Prepare to proceed with a single session.
        input.run = [x y]; % Current run, to pass to functions.
        [input, opt] = prepforsession(input, opt);

        %% 03. Proceed to appropiated pipeline.
        switch input.sessions(input.run(1)).info.fileformat
            case {'DT2', 'DF1'} 
               % 03.1 Deuteron Pipeline
               [input, opt] = Deuteron_PipelineWrapper(input, opt); 
    
            case {'fileperch', 'filepertype', 'tradFormat'}
               % 03.2 INTAN Pipeline
               [input, opt] = INTAN_PipelineWrapper(input, opt);

            case {'FieldTrip'}   
               % 03.3 FT Pipeline
               % Check for events, neurons and spike variables.
               [events, trialdef, EventRecord] = EventProcess(input, opt);

               % so far, reaching this point means there was no raw data,
               % and under analysis there is FT formatted data, so prob
               % this has already been preprocessed but we only have the
               % minimal data here for analysis.
               disp('Session skipped because continuous FT file was found')
               continue
                
            otherwise
               warning('Something went wrong during format verification. Skipping Session');
               continue
        end 
        
        %% 04. Kilosort
        if opt.kilosort
            % Kilosort 4 will run without GUI.
            master_kilosort4(input, opt)
            close all
        end

        %% 05. Bombcell
        if opt.bombcell
            % Kilosort will run without GUI.
            Bombcell_Main(input, opt) 
        end

        %% 06. Open Phy to manual curation or just inspection
        if opt.phy
            % Will change to current session directory and open phy.
            % ! Keeps MATLAB busy until interface is closed.
            cd(opt.FolderProcDataMat)
            system('phy template-gui params.py');
        end

        %% Clean up to move on to next session
        clear FT_data INTANdata txt events EventRecord trialdef

    end % sessions loop
end % subjects loop