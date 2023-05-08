function sessions = findSessions(input)
% Finds and list all sessions requested, no matter the input format.
%
% Version 04.04.2023 Jesus

% Goes over every subject's folder and reads the sessions    
for s = 1:input.nsubjects
    cd(fullfile(input.datafolder, string(input.subjects(s).name)))
    ss = dir();
    dirFlags = [ss.isdir];
    ss = ss(dirFlags);
    ss(ismember({ss.name}, {'.', '..'})) = [];

        % Check 'all' vs explicit sessions request
        if iscell(input.dates)
            nameFlags = ismember({ss.name}, input.dates);
            ss = ss(nameFlags);
        end

    sessions(s).nsessions = length(ss);
    sessions(s).folder = ss(1).folder;
    
    for d = 1:sessions(s).nsessions
        sessions(s).list{d} = ss(d).name;
    end
end

end