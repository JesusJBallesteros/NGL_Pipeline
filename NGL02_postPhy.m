%% NGL02_postPhy (in progress)
% To run after manual curation of desired sessions is completed. Will read
% the resulting KS results after manual curation.
%
% Jesus 12.06.2024
if ~isfield(opt, 'doSpikething') || isempty(opt.doSpikething),  opt.doSpikething = true;    end
if ~isfield(opt, 'doLFPthing') || isempty(opt.doLFPthing),      opt.doLFPthing   = true;    end
if ~isfield(opt, 'offlineTrack') || isempty(opt.offlineTrack),  opt.offlineTrack = false;   end
if ~isfield(opt, 'FLIP') || isempty(opt.FLIP),                  opt.FLIP         = false;   end
% if exist('regions','var'),                                      opt.multregion   = true;    end

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
else
    disp('Using INPUTS from NGL01_MAIN.')
end

% Set default inputs and dependencies. In case NGL01 did not before.
input = set_default(input, opt);
if ~isfield(opt,'useTrack'), opt.useTrack = false; end

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

%% 02. Proceed with data per session
for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
            %% 2.01. Prepare to proceed with a single session.
            input.run = [x y]; % Current run, to pass to functions.
            [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);           

            %% 2.02. Offline Video blob detector
            % Very specific for Social learning videos from central cenital camera. 
            if opt.offlineTrack
                [blob] = processAndTrack_video(opt);
                save(fullfile(opt.analysis, "blob.mat"), 'blob');
            end

            %% 2.03. SPIKE DATA
            if opt.doSpikething
                % Extract preprocessed 'spike' and recover 'events' data.
                if exist(fullfile(opt.spikeSorted, "spike.mat"),'file')
                    load(fullfile(opt.spikeSorted, "spike.mat"));
                else
                    % Spike clusters after sorting and curation.
                    spike = loadSpikes(opt);
                    if isfield(spike,"spike"), spike = spike.spike; end % Simplify loaded structure if needed
                            
                    % Save output to \spikesorted
                    save(fullfile(opt.spikeSorted, "spike.mat"), 'spike', '-mat');
                end
                
                % Sort 'spike' into 'trialdef' to create 'neurons'
                if ~exist(fullfile(opt.analysis, "neurons.mat"),'file')
                    % Iterate trough all units and sort them into trials.
                    % Recover trial definitions created after event extraction and processing. 
                    % Can have as many variations as requested at that time.
                    % To create new alignments, it would have to be ran again.
                    if ~exist('trialdef','var'), load(fullfile(opt.trialSorted, "trialdef.mat")); end
                    
                    % Outputs are saved to data\analysis. 
                    % TODO fix Fieldtrip extraction
                    [neurons, ~] = sort2trials(spike, trialdef, opt);

                    % Save output to \analysis
                    save(fullfile(opt.analysis, "neurons.mat"), 'neurons', '-mat')
                end

                % Calculate fire rate and normalized fire rate
                if ~exist(fullfile(opt.analysis, "fireRate.mat"),'file')
                    if ~exist('neurons','var'), load(fullfile(opt.analysis, "neurons.mat")); end
                    if ~exist('events','var'), load(fullfile(opt.analysis, "events.mat")); end
                    if ~exist('condition','var'), load(fullfile(opt.analysis, "condition.mat")); end
    
                    % General function, no conditions: 'allInitiated' by default
                    fireRate = calculate_fireRate_general(neurons, [], conditions, opt, param);
                        
                    save(fullfile(opt.analysis, "fireRate.mat"), 'fireRate', '-mat')

                    % % Project specific    
                    % param.IncludeFS = true; % NS and FS
                    % param.trial2plot = 'allInitiated'; % for correct trials
                    % fireRate = calculate_fireRate_extintion(neurons, events, conditions, opt, param);
                    % save(fullfile(opt.analysis, "fireRate_extintion.mat"), 'fireRate', '-mat')
                end
                
                % calculate dynamics
                %
                %

                % Tracking in Social Arena
                if opt.useTrack
                    % Spiking indexing for Social intereactions. Checks blob
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
                    % Read excel file, num trials
                    if ~exist('events','var'), load(fullfile(opt.analysis, "events.mat")); end
                    if ~isfield(events,"social")
                        [events.social] = getSocialEvents(opt);
                        save(fullfile(opt.analysis, "events.mat"), 'events', '-mat')
                    end
                end

            end

            %% 2.04. Continuous LFP DATA. UNDER DEVELOPMENT
            if opt.doLFPthing
                % Extract preprocessed FTcont file.
                if opt.trialparsed
                   disp('Loading FT trial parsed file...')
                   load(fullfile(opt.trialSorted, [opt.SavFileName '_stimOn2.mat']), "-mat", 'FT_data');
                else, disp('Loading FT continuous file...')
                      load(fullfile(opt.trialSorted, [opt.SavFileName '_FTcont.mat']), "-mat", 'FT_data');
                      FT_data.cfg.continuous = 'yes';
                end
                
                if isfield(FT_data,"FT_data"), FT_data = FT_data.FT_data; end   % Simplify loaded structure if needed

                % Obtain or create trial definition to pass to FT
                if isfile(fullfile(opt.trialSorted, 'trialdef.mat')), load(fullfile(opt.trialSorted, "trialdef.mat")); end
                if isfile(fullfile(input.analysis, 'data_all.mat'))
                    load(fullfile(input.analysis, "data_all.mat"), 'allconditions');
                    conditions = allconditions{x,y};
                end

                if ~exist('trialdef','var')
                    load(fullfile(opt.analysis, "events.mat"));
                    if exist('events','var')
                        [~, trialdef, ~] = trialdefGen(events, opt, 1);
                        save(fullfile(opt.trialSorted, "trialdef.mat"), 'trialdef');
                    else, warning('Neither trial definitions or events found for this session.')
                    end
                end

                % % Specific for chgDtctPCue
                % trialdef{2,1} = ceil(trialdef{2,1}/32);                
                % % Outputs are saved to data\analysis. 
                % MAT2FieldTrip(FT_data, opt, trialdef); 
                % clear FT_data events trialdef

                % Artifact detection and rejection. IN PROGRESS
                if opt.artifdet
                    FT_data = artifact_detRej_lfp(FT_data, condition, opt, param);
                end

                % Time-frequency analisys. IN PROGRESS
                if opt.spectrogram
                    if strcmp(FT_data.cfg.continuous, 'yes')
                        allTFR_continuous{x,y} = continous_MTspectrogram(FT_data, conditions, param, opt);

                        % Save once all iterations are done
                        if y == input.sessions(x).nsessions
                            save(fullfile(input.analysis,'allTFR_continuous.mat'), 'allTFR_continuous', 'param', 'opt', '-mat');
                        end
                    else
                        % Specific parameters
                        param.testname = 'ASL_Clean_Final_Correct';

                        % pass the analysis folder info too % TO optimize
                        opt.analysisCode = input.analysisCode;

                        [allTFR_trialparsed{x,y}, TFRcgf{x,y}] = trialparsed_MTspectrogram(FT_data, conditions, param, opt);
                        
                        % Save once all iterations are done
                        if y == input.sessions(x).nsessions
                            save(fullfile(input.analysis,'allTFR_chbych_trialparsed.mat'), 'allTFR_trialparsed', 'param', 'opt', '-mat');
                        end
                    end
                    
                end
    
                % % vFLIP Analysis
                % if opt.FLIP
                %      laminaraxis = 0:0.05:1.55;
                %      freqaxis = 1:150;
                % 
                %     % Specific for chgDtctPCue
                %     [FLIP, relpow, ~] = vFLIP_NGL(FT_data, laminaraxis, freqaxis, 0);
                % 
                % end

            end

        %% Clean up to move on to next session
        clear events EventRecord trialdef blob neurons spike
     
    end
end

