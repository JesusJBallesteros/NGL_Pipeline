%% NGLXX_postPhy (in progress)
% To run after manual curation of desired sessions is completed. Will read
% the resulting KS results after manual curation.

%% Input storage drive, project name and toolbox folder:
input.datadrive     = 'F:\';
input.studyName     = 'ephysLabComparison'; %'ephysTestBundleWires';
input.toolbox       = 'C:\Code\ephys-data-pipeline'; % Default: 'C:\Code\Scripts\ephys-data-pipeline'

% To run the script on all subjects and sessions included in your project,
% just leave both as 'all'. For a session-to-session process, explicit the
% subject and session/s to process. 
input.subjects       = 'all'; % char array 'all', or a single subject denomination e.g. 'DOE'
input.dates          = 'all'; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}

%% General Options. What you want to obtain:
% Those used for all sessions. Specific options can be set below or defaulted in the functions.
opt = struct();
    opt.plotdrift      = false; % logic to trigger plot
    opt.plotAmpDepth   = false; % logic to trigger plot
    opt.excludeNoise   = true; % param for plots                 
    opt.loadPCs        = false; % param for plots                 
    
%% 00. Check inputs, set defaults and dependencies.
set_default(input);
opt.postPhy = true; % Always true, to differentiate from NGL01

for s = 1:input.nsubjects
    % Read requested sessions from specified animal folder.
    sessions = findSessions(input, opt);

    %% 01. Loop subjects and sessions to process.
    for ss = 1:sessions(s).nsessions
        results{ss,s} = struct;

        % Navigate to session's raw data folder.
        cd(fullfile(sessions(s).folder,sessions(s).list{ss}));
                
        % Progress report.
        fprintf(sprintf('\n --> Subject %s, session %d out of %d: %s \n', ...
                     input.subjects(s).name, ss, sessions(s).nsessions, sessions(s).list{ss}));

        % Determine where processed data will be saved, done for every session.
        opt.PathRaw           = pwd;
        opt.FolderProcDataMat = fullfile(input.sorted, input.subjects(s).name, sessions(s).list{ss});
        opt.SavFileName       = sessions(s).list{ss}; 
        
        % Report and create folder.
        disp(strcat('Processed data will be saved to: >', opt.FolderProcDataMat));
        if exist(opt.FolderProcDataMat,"dir") == 0
            mkdir(opt.FolderProcDataMat);
        end

        %% 02. Proceed with reading data from preprocessed files
        % Load information as stored post-Phy curation.
        [results{ss,s}.spike] = loadKSdir(opt.PathRaw); 

        %% 03? Get 'events.mat' and generate 'conditions' struct
        % TODO: trials needs to be readed from an event mat file created
        % after a proper experiment. Not valid for testing recordings. 
        
        % If 'trials' does not exists, defaults here to 1 (single, long trial)
          % smth like... if exists trials then use it, otherwise, trials = 1 
          trials = 1;

        % TODO: create conditions
        % This variable will be saved under '...\data\spikeSorted\...' for further access

        %% 03. Here we can create the 'neurons' cell variable, according to the IKN standard.
        % This variable will be saved under '...\data\spikeSorted\...' for further access
        getneurons(results{ss,s});

        %% 04. For now, we can create some plots using those created by 'Cortex-Lab' at ULC, for instance.
        if ~isempty(results{ss,s}.spike.spikeTemplates)
            % Will update some spike parameter and create the template structure.
            [results{ss,s}.spike, results{ss,s}.template] = plot_KSresults(results{ss,s}.spike, opt);
        end

    end
end
