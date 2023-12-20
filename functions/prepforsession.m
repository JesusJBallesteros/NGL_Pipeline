function [info, opt] = prepforsession(input, sessions, opt, ss)

% Navigate to session's raw data folder.
cd(fullfile(sessions.folder,sessions.list{ss}));
        
% Report.
txt = sprintf('\n --> Subject %s, session %d out of %d: %s \n', ...
             input.subjects(ss).name, ss, sessions.nsessions, sessions.list{ss});
fprintf(txt);

% Check system and version.
info = chckV();

% Determine where processed session data will be saved.
opt.PathRaw           = pwd;
opt.FolderProcDataMat = fullfile(input.processed, input.subjects(ss).name, sessions.list{ss});
opt.behavFiles        = fullfile(input.bhvfolder, input.subjects(ss).name, sessions.list{ss});
opt.SavFileName       = sessions.list{ss}; 

% Report and create folder.
disp(strcat('Processed data will be saved to: >', opt.FolderProcDataMat));
mkdir(opt.FolderProcDataMat);

end