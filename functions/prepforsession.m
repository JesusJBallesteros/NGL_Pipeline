function [input, opt] = prepforsession(input, opt)
% prepforsession. Per-session setup: format detection, path assignment, folder creation.
%
% PURPOSE:
%   Called at the start of each session iteration in NGL01_Main.
%   Navigates to the raw data folder, checks the file format, builds
%   all session-specific output paths in 'opt', creates missing output
%   directories, reads the INTAN header (if INTAN), and generates a
%   reduced channel map if nfiles is different from opt.numChannels.
%
% USAGE:
%   [input, opt] = prepforsession(input, opt)
%   Called after input.run = [subject_index, session_index] ONLY
%
% INPUTS:
%   input  - struct from set_default, with input.run = [x y] indicating the
%            current subject (x) and session (y)
%   opt    - options struct (complete, from set_default)
%
% OUTPUTS:
%   input  - updated with input.sessions(x).info populated with format metadata
%   opt    - updated with session-specific paths:
%              .PathRaw          raw data folder (current dir)
%              .SavFileName      session name string
%              .FolderProcDataMat  preprocessing output folder
%              .KSfolder         Kilosort output folder
%              .behavFiles       behaviour folder
%              .spikeSorted      spike-sorted output folder
%              .trialSorted      trial-sorted output folder
%              .analysis         analysis output folder
%            and, if channel mismatch detected:
%              .KSchanMapFile    updated to reduced map filename
%              .numChannels      updated to actual file count
%
% CALLS:
%   chckV, findSetting, reduceChanMap
%
% Last modified 06.05.2026 (Jesus)

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
input.sessions(input.run(1)).info = chckV();

% Collect data to create paths.
opt.PathRaw             = pwd;
opt.SavFileName         = session; 

% Create paths to session-specific folders
opt.FolderProcDataMat   = fullfile(input.processed, subject, session);
opt.KSfolder            = fullfile(opt.FolderProcDataMat, 'kilosort', int2str(opt.kilosort));
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

% Brought here from Intan wrapper, so all header info is available already
if contains(input.sessions(input.run(1)).info.fileformat,'fileper')
    % Find out INTAN settings and header file. Extract info.
    %  Uses a modified Intan function, to make the basic information
    %  available at 'info{ss}' and a more detailed info at
    %  the '.INTAN_hdr' sub-structure.
    input = findSetting(input);
end

% Check number of expected channels vs number of raw files. Create a 
% reduced channel map if mismatched, and save in preprocessing output dir.
% Only for INTAN (07.01.2026, Winston)
if contains(input.sessions(input.run(1)).info.fileformat,'fileper')
    if length(dir('amp*.dat')) > opt.numChannels
        error('More INTAN files than number of channels specified in NGL_SetAndRunMe.m')
    elseif length(dir('amp*.dat')) < opt.numChannels
            % Create reduced probe map .mat file for kilosort
            opt = reduceChanMap(input,opt); % also reassigns opt.KSchanMapFile to new map
            opt.numChannels = length(dir('amp*.dat'));
    end
end

disp('Single session paths created, and their formats and settings extracted.')
end