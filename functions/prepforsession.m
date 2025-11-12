function [info, opt] = prepforsession(input, opt)
% Check system and version. Determine where processed session data will be saved.
if ~isfield(opt, 'kilosort') || isempty(opt.kilosort),   opt.kilosort = 4; end
% Version 10.10.2025 (Winston)

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
if ~exist(fullfile(opt.FolderProcDataMat),"dir"), mkdir(opt.FolderProcDataMat); end
if ~exist(fullfile(opt.behavFiles),"dir"), mkdir(opt.behavFiles); end
if ~exist(fullfile(opt.spikeSorted),"dir"), mkdir(opt.spikeSorted); end
if ~exist(fullfile(opt.trialSorted),"dir"), mkdir(opt.trialSorted); end
if ~exist(fullfile(opt.analysis),"dir"), mkdir(opt.analysis); end

% Check number of expected channels vs number of raw files. Create a 
% reduced channel map if mismatched, and save in preprocessing output dir.
if ~isempty(dir('amp*.dat')) %skip if Deuteron is the recording system 
    if length(dir('amp*.dat')) > opt.numChannels
        error('More INTAN files than number of channels specified in NGL_SetAndRunMe.m')
    elseif length(dir('amp*.dat')) < opt.numChannels
        % Create reduced probe map .mat file for kilosort
        opt = reduceChanMap(input,opt); % also reassigns opt.KSchanMapFile to new map
        opt.numChannels = length(dir('amp*.dat'));
    end
end