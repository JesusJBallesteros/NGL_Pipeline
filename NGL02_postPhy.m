%% NGLXX_postPhy (in progress)
% To run after manual curation of desired sessions is completed. Will read
% the resulting KS results after manual curation.

%% Input storage drive, project name and toolbox folder:
input.datadrive     = 'D:\';
input.studyName     = 'ephysTestATLAS'; % For SPP people: 'Dorian\SPP'
input.toolbox       = 'C:\Code\Scripts\ephys-data-pipeline'; % Default: 'C:\Code\Scripts\ephys-data-pipeline'

% To run the script on all subjects and sessions included in your project,
% just leave both as 'all'. For a session-to-session process, explicit the
% subject and session/s to process. 
input.subjects       = {'478'}; % char array 'all', or a single subject denomination e.g. 'DOE'
input.dates          = {'20230515' '20230516' '20230517'}; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}

%% General Options. What you want to obtain:
% Those used for all sessions. Specific options can be set below or defaulted in the functions.
opt = struct();
    opt.postPhy        = true; % Always true, to differentiate from NGL01

    % Optionals
    opt.plotdrift      = true; % testing
    opt.plotAmpDepth   = true; % testing
        
%% 00. Check inputs, set defaults and dependencies.
set_default(input);

for s = 1:input.nsubjects
    % Read requested sessions from specified animal folder.
    sessions = findSessions(input, opt);

    %% 02. Loop subjects and sessions to process.
    for ss = 1:sessions(s).nsessions
        % Navigate to session's raw data folder.
        cd(fullfile(sessions(s).folder,sessions(s).list{ss}));
                
        % Progress report.
        txt = sprintf('\n --> Subject %s, session %d out of %d: %s \n', ...
                     input.subjects(s).name, ss, sessions(s).nsessions, sessions(s).list{ss});
        fprintf(txt);

        % Determine where processed data will be saved, done for every session.
        opt.PathRaw           = pwd;
        opt.FolderProcDataMat = fullfile(input.sorted, input.subjects(s).name, sessions(s).list{ss});
        opt.SavFileName       = sessions(s).list{ss}; 
        
        % Report and create folder.
        disp(strcat('Processed data will be saved to: >', opt.FolderProcDataMat));
        mkdir(opt.FolderProcDataMat);

        %% 03. Proceed with reading data from preprocessed files
        [spike, template] = read_KSresults(opt);


    end
end
