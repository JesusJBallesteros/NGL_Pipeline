function [info, opt] = prepforsession(input, opt)
% Check system and version. Determine where processed session data will be saved.

% Version 12.06.2024 (Jesus)

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

% Collect data to create paths.
opt.PathRaw             = pwd;
opt.SavFileName         = session; 

% Create paths to session-specific folders
opt.FolderProcDataMat   = fullfile(input.processed, subject, session);
opt.KSfolder            = [opt.FolderProcDataMat, '\kilosort', int2str(opt.kilosort)];
opt.behavFiles          = fullfile(input.bhvfolder, subject, session);
opt.spikeSorted         = fullfile(input.spikeSorted, subject, session);
opt.trialSorted         = fullfile(input.trialSorted, subject, session);
opt.analysis            = fullfile(input.analysis, subject, session);

% Create session-specific folders.
mkdir(opt.FolderProcDataMat);
mkdir(opt.behavFiles);
mkdir(opt.spikeSorted);
mkdir(opt.trialSorted);
mkdir(opt.analysis);

end