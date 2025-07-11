function sessions = findSessions(input, varargin)
% Finds and list all sessions requested, no matter the input format.

% Version 28.05.2024 (Jesus)

if ~exist("sessions","var")
    % Goes over every subject's folder and reads existing sessions
    for s = 1:input.nsubjects
        if exist(fullfile(input.datafolder, string(input.subjects(s).name)),"dir")
            cd(fullfile(input.datafolder, string(input.subjects(s).name)))
        else
            cd(fullfile(input.analysis, string(input.subjects(s).name)))
        end
        
        % Select sessions
        ss = dir(); % List all content in folder
        dirFlags = [ss.isdir]; % Set to keep only directories
        ss = ss(dirFlags); % Keep only directories
        ss(ismember({ss.name}, {'.', '..'})) = []; % Remove spurious '.' and '..'
    
        % Check 'all' vs explicit sessions request
        if iscell(input.dates) % Request is a subset
            nameFlags = ismember({ss.name}, input.dates); % Index subset
            ss = ss(nameFlags); % Keep indexed subset
        end
    
        % Get number of sessions added and its folder.
        sessions(s).nsessions = length(ss);
        sessions(s).folder = ss(1).folder;
        
        % Create a definitive list of sessions
        for d = 1:sessions(s).nsessions
            sessions(s).list{d} = ss(d).name;
        end
    end
end
end