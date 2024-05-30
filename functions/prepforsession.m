function [info, opt] = prepforsession(input, opt)
% Check system and version. Determine where processed session data will be saved.

% Version 28.05.2024 (Jesus)

% Extract subject and session 
subject = input.subjects(input.run(1)).name;
session = input.sessions(input.run(1)).list{input.run(2)};

% Navigate to session's raw data folder.
cd(fullfile(input.sessions(input.run(1)).folder, session));
        
% Report.
txt = sprintf('\n --> Subject %s, session %d out of %d: %s \n', ...
             subject, input.run(2), input.sessions(input.run(1)).nsessions, session);
fprintf(txt);

% Check system and version.
info = chckV();

% Determine where processed session data will be saved.
opt.PathRaw             = pwd;
opt.SavFileName         = session; 
opt.FolderProcDataMat   = fullfile(input.processed, subject, session);
opt.behavFiles          = fullfile(input.bhvfolder, subject, session);
opt.KSfolder            = [opt.FolderProcDataMat, '\kilosort', int2str(opt.kilosort)];
opt.spikeSorted         = fullfile(input.spikeSorted, subject, session);
opt.trialSorted         = fullfile(input.trialSorted, subject, session);

% Create folder.
mkdir(opt.FolderProcDataMat);

end