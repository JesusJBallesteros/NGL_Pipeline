function sessions = findSessions(input)
% Finds and list all sessions requested, no matter the input format.
%
% Version 01.03.2023 Jesus

if iscell(input.dates) % input is cell array of dates
    input.dates = input.datafolder + input.animal + "\" + input.animal + "_" + input.dates(:);
    [sessions.folder,sessions.name,~] = fileparts(input.dates);
    sessions.folder = unique(sessions.folder);

    % Get and Count sessions
    cd(sessions.folder)
    sessions.nSessions = length(sessions.name);
    for s = 1:sessions.nSessions
        sessions.list(s) = string(sessions.name(s));
    end
    clear s

elseif strcmp(input.dates, 'all') % input is 'all'
    sessions.folder     = fullfile(input.datafolder, input.animal);
    sessions.list       = dir(sessions.folder + '\' + input.animal + '*');
    sessions.nSessions  = length(sessions.list);
    
    cd(sessions.folder)
end

end
