function sessions = findSessions(input)
% UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if iscell(input.dates) % input is cell array of dates
    input.dates = input.datafolder + input.animal + "\" + input.animal + "_" + input.dates(:) + '*';
    [sessions.folder,sessions.name,~] = fileparts(input.dates);
    sessions.folder = unique(sessions.folder);

    % Get and Count sessions
    cd(sessions.folder)
    sessions.nSessions = length(sessions.name);
    for s = 1:sessions.nSessions
        sessions.list(s) = dir(sessions.name(s));
    end
    clear s

elseif strcmp(input.dates, 'all') % input is 'all'
    sessions.folder     = fullfile(input.datafolder, input.animal);
    sessions.list       = dir(sessions.folder + '\' + input.animal + '*');
    sessions.nSessions  = length(sessions.list);
    
    cd(sessions.folder)
end

end