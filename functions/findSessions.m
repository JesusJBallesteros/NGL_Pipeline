function sessions = findSessions(input)
% findSessions. Discover and list session directories for all requested subjects.
%
% PURPOSE:
%   Scans the raw data and returns a struct array describing the sessions
%   that match the subjects and dates requested in 'input'.
%   Works for all recording formats (INTAN, Deuteron, FieldTrip).
%
% USAGE:
%   sessions = findSessions(input)
%
% INPUT:
%   input  - struct built by set_default, must contain:
%              .datafolder   (char) root of data/raw/
%              .analysis     (char) fallback folder if raw does not exist
%              .subjects     (dir-struct array) subjects to process
%              .nsubjects    (scalar) number of subjects
%              .dates        (char 'all' | cell of 'YYYYMMDD' strings)
%
% OUTPUT:
%   sessions - (1 × nsubjects) struct array, each element with:
%                .nsessions  (scalar) number of sessions found
%                .folder     (char)   parent folder of those sessions
%                .list       (1 × nsessions cell of char) session names
%
% NOTES:
%   - If the subject folder does not exist under datafolder, falls back to
%     the analysis folder (supports partially-preprocessed data trees).
%   - Session names must match 'YYYYMMDD' format.
%
% Last modified 06.05.2026 (Jesus)

if ~exist("sessions", "var")
    % Goes over every subject's folder and reads existing sessions
    for s = 1:input.nsubjects
        if exist(fullfile(input.datafolder, string(input.subjects(s).name)), "dir")
            cd(fullfile(input.datafolder, string(input.subjects(s).name)))
        else
            cd(fullfile(input.analysis, string(input.subjects(s).name)))
        end
        
        % Select sessions
        ss = dir(); % List all content in folder
        dirFlags = [ss.isdir]; % Set to keep only directories
        ss = ss(dirFlags); % Keep only directories
        ss(ismember({ss.name}, {'.', '..'})) = []; % Remove '.' and '..' from the list
    
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

disp('Subjects and sessions found and listed.')
end