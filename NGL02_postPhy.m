%% NGL02_postPhy (in progress)
% To run after manual curation of desired sessions is completed. Will read
% the resulting KS results after manual curation.
%
% Jesus 12.06.2024
if ~isfield(opt, 'doSpikething') || isempty(opt.doSpikething),  opt.doSpikething = true;    end
if ~isfield(opt, 'doLFPthing') || isempty(opt.doLFPthing),      opt.doLFPthing   = true;    end
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

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

%% 02. Load Post_Phy parameter file (under '/analysisCode')
try run(fullfile(input.analysisCode, 'postPhy_param.m'));
catch, error('No script found with parameters for the Post-phy pipeline. Find it and locate it into your analysisCode folder.')
end

%% 03. Proceed with data per session
for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
            input.run = [x y]; % Current run, to pass to functions.
            
            %% 04. Prepare to proceed with a single session.
            [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);           

            %% 05. Offline Video blob detector
            % Very specific for Social learning videos from central cenital camera. 
            if opt.offlineTrack
                [blob] = processAndTrack_video(opt);
                save(fullfile(opt.behavFiles, "blob.mat"), 'blob');
            end

            %% 06. SPIKE DATA
            if opt.doSpikething
                % 04.1 Extract preprocessed spikes and recover event data.
                % Spike clusters after sorting and curation.
                spike = loadSpikes(opt);
                if isfield(spike,"spike"), spike = spike.spike; end % Simplify loaded structure if needed
                            
                % Save output to \spikesorted
                save(fullfile(opt.spikeSorted, "spike.mat"), 'spike', '-mat');
                
                % 04.2 Iterate trough all units and sort them into trials.
                % Recover trial definitions created after event extraction and processing. 
                % Can have as many variations as requested at that time.
                % To create new alignments, it would have to be ran again.
                if ~exist('trialdef','var'), load(fullfile(opt.trialSorted, "trialdef.mat")); end
                
                % Outputs are saved to data\analysis. 
                % TODO fix Fieldtrip extraction
                [neurons, ~] = sort2trials(spike, trialdef, opt);
            
                % 04.3 Save output to \analysis
                save(fullfile(opt.analysis, "neurons.mat"), 'neurons', '-mat')
            end

            %% 07. Continuous LFP DATA. UNDER DEVELOPMENT
            if opt.doLFPthing
                % 05.1 Extract preprocessed FTcont file.
                disp('Loading FT continuous file...')
                load(fullfile(opt.analysis, input.sessions(input.run(1)).info.files.name));
                if isfield(FT_data,"FT_data"), FT_data = FT_data.FT_data; end   % Simplify loaded structure if needed

                % 05.2 Obtain or create trial definition to pass to FT
                if isfile('trialdef.mat'), load(fullfile(opt.trialSorted, "trialdef.mat")); end
                
                if ~exist('trialdef','var')
                    load(fullfile(opt.analysis, "events.mat"));

                    if exist('events','var')
                        [~, trialdef, ~] = trialdefGen(events, opt, 1);
                        save(fullfile(opt.trialSorted, "trialdef.mat"), 'trialdef');
                    else, warning('Neither trial definitions or events found for this session.')
                    end
                end

                % Specific for chgDtctPCue
                trialdef{2,1} = ceil(trialdef{2,1}/32);                
                % Outputs are saved to data\analysis. 
                MAT2FieldTrip(FT_data, opt, trialdef); 
                clear FT_data events trialdef

                % Artifact detection and rejection
                load(fullfile(opt.analysis, [input.sessions(input.run(1)).list{input.run(2)}, '_startON.mat']), "-mat", 'FT_data');
                
                cfg = [];
                cfg.trl         = FT_data.cfg.trl;
                cfg.continuous  = 'no';
                cfg.artfctdef.zvalue.channel    = 'all';
                cfg.artfctdef.zvalue.cutoff     = 20;
                cfg.artfctdef.zvalue.trlpadding = 0;
                cfg.artfctdef.zvalue.fltpadding = 0;
                cfg.artfctdef.zvalue.artpadding = 0;                
                
                % The optional configuration settings (see below) are:
                  cfg.artfctdef.zvalue.artfctpeak       = 'no';
                  cfg.artfctdef.zvalue.interactive      = 'no';
                  cfg.artfctdef.zvalue.zscore           = 'yes';
                
                [~, artifact] = ft_artifact_zvalue(cfg, FT_data);
    
                % The following configuration options are supported
                cfg = [];
                  cfg.artfctdef.reject          = 'partial';
                  cfg.artfctdef.zvalue.artifact = artifact;
    %               cfg.artfctdef.minaccepttim    = when using partial rejection, minimum length
    %                                               in seconds of remaining trial (default = 0.1)
                [FT_data_art] = ft_rejectartifact(cfg, FT_data);
    
                % vFLIP Analysis
                if opt.FLIP
                     laminaraxis = 0:0.05:1.55;
                     freqaxis = 1:150;
    
                    % Specific for chgDtctPCue
                    [FLIP, relpow, ~] = vFLIP_NGL(FT_data, laminaraxis, freqaxis, 0);
    
                end
            end
            
    end
end