%% Pipeline process INTAN and Deuteron continous data.
% Will read and process INTAN, DEUTERON (or ALLEGO) data, from selected sessions for a given animal.
% The main pipeline will be: INTAN/DEUTERON raw formats to be located, then
% converted to .bin files (spike sorting), and Fieldtrip .mat structures
% (for LFP). Once sorted, spike data will be attached to the FieldTrip
% structure. For Arena experiments, motion sensor data will be extracted and
% interpreted. Data will be trial-parsed using EventCodes.
%
% The hard disk data structure SHOULD fit the IKN standard published at:
% gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure
%
% DEPENDENCIES:
% Requires that all pipeline dependencies are properly located. 
% I suggest to include the 'mainfolder' in Matlab's permanent path system.
% The function 'set_default' will take care of the rest of folders on each run.
%
% OUTPUTS:
% For one single session or for a batch of sessions, from one single animal:
%       Fieldtrip (.mat), binary (.bin) and/or HDF5 (.h5) and .nwb files from
%           1. Deuteron .DT2 or .DF1 data.
%           2. INTAN file-per-type and file-per-channel format data.
%           3. (ALLEGO data?)
%       EventRecord.mat file, from Deuteron session.
%       MotionData.mat file, From Deuteron sensors.
%       Plots snippets of time- and frequency-domain data, from FieldTrip
%       
% Last modified 05.03.2024 (Jesus)

% TODO LIST
% Order channels as incremental ordinals. (in 'Deuteron_ExtractEvents' ~136)
% If Deuteron2Kilosort(opt) filter for DF1 format works, set filter out of format cases (generalize)
% Continue with 'Deuteron_GetDigInEvents' when we get a recording with EVENTS
% Check for Deuteron_GetDigInEvents(EventRecord) status.
% Check for FT trial-parsing using EventRecord with MAT2FieldTrip(data, opt, varargin)
%    Create a 'trial-parsed' stream in 'mat2FieldTrip' VS. add post-hoc parsing
% Extract nChannels from EventsRecord. Find first 'File started' then use
%    'strsplit(EventRecord(50).Details,{';','='})' and find the 6th cell
% There seems to be an ERROR on 2nd and following runs of the NWB functionalities.
%    Figure out what's going on with the NWB/H5 DLLs that block either when the other has been performed...

% Version 28.05.2024 (Jesus)

%% 00. Check current inputs.
% Check if input variable exist already. Parse values.
if ~exist("input","var")
    input = struct( 'datadrive' , datadrive , ...   % force char array
                    'studyName' , studyname , ...   % force char array
                    'toolbox'   , toolbox   , ...   % force char array
                    'subjects'  , [], ...           % do NOT force char array
                    'dates'     , []        );      % do NOT force char array
    input.dates     = dates;    % place as it comes
    input.subjects  = subjects; % place as it comes
end

% Set default inputs and dependencies.
input = set_default(input, opt);

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        %% 02. Prepare to proceed with a single session.
        input.run = [x y]; % Current run, to pass to functions.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);

        %% 03. Proceed to appropiated pipeline.
        switch input.sessions(input.run(1)).info.fileformat
            case {'DT2', 'DF1'} 
               % 03.1 Deuteron Pipeline
               opt = Deuteron_PipelineWrapper(input, opt);
    
            case {'fileperch', 'filepertype'}
               % 03.2 INTAN Pipeline
               input = INTAN_PipelineWrapper(input, opt);

            case {'FieldTrip'}   
               % 03.3 FT Pipeline
               % Check for events, neurons and spike variables.
               [events, trialdef, EventRecord] = EventProcess(opt);

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
        if opt.kilosort == 2
            % Kilosort 2 will run without GUI.
            master_kilosort(input, opt)
        elseif opt.kilosort == 4
            % Kilosort 4 will run without GUI.
            master_kilosort4(input, opt)
        end
        close all
        
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