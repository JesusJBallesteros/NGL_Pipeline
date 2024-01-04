function [info, opt] = prepforsession(input, opt)
% 
% Version 02.01.2024 (Jesus)

% Navigate to session's raw data folder.
cd(fullfile(input.sessions(input.run(1)).folder, input.sessions(input.run(1)).list{input.run(2)}));
        
% Report.
txt = sprintf('\n --> Subject %s, session %d out of %d: %s \n', ...
             input.subjects(input.run(1)).name, input.run(2), input.sessions(input.run(1)).nsessions, input.sessions(input.run(1)).list{input.run(2)});
fprintf(txt);

% Check system and version.
info = chckV();

% Determine where processed session data will be saved.
opt.PathRaw           = pwd;
opt.FolderProcDataMat = fullfile(input.processed, input.subjects(input.run(1)).name, input.sessions(input.run(1)).list{input.run(2)});
opt.behavFiles        = fullfile(input.bhvfolder, input.subjects(input.run(1)).name, input.sessions(input.run(1)).list{input.run(2)});
opt.SavFileName       = input.sessions(input.run(1)).list{input.run(2)}; 

% Report and create folder.
disp(strcat('Processed data will be saved to: >', opt.FolderProcDataMat));
mkdir(opt.FolderProcDataMat);

end