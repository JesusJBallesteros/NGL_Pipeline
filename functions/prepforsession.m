function [info, opt] = prepforsession(input, sessions, opt, ss)
% This function prepares to run a single session dynamically.
% INPUTS:
%       input: struct with basic info about paths, project and current run
%       sessions: struct with listed sessions per subject.
%       opt: struct with options defined by user or defaulted.
%       ss: integer as ordinal subject to run.
% OUTPUT:
%       info: struct generated with single session values. Ideally gets integrated
%             into the input variable 'sessions'
%       opt: struct updated with especific single-session paths.
%
% Jesus. 21.12.2023

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