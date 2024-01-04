% Input file where there is no access to any of the running code, making
% this the only file that needs to be modified, and that could call all
% pipelines as a sequence of easily swichable runs by simply commenting 
% lines. Could be ran line-by-line (F9) or all at once (F5).
%
% Once the folder system is created, should this file be located INSIDE?
% i.e. at '...\studyName\data\analysis' ??
%
% Jesus 04.01.2024

%% INPUT
% A) CRITICAL
% Specify drive where data is located, Project Name and the toolbox folder.
% RECOMMENDED to add the toolbox folder to MATLAB folder system, but not necessary.
input.datadrive     = 'F';  % Define the LETTER of the drive where the data structure will be created or already exists 
input.studyName     = 'Pilot_SocialLearning';  % Name the Study or Project to be used.
input.toolbox       = 'C:\Code\ephys-data-pipeline'; % Absolute address to the toolbox.

% B) SUBJECTS AND SESSIONS
% To run the script on all subjects and sessions, or as session-to-session process.
input.subjects      = 'all'; % {'257'};      % char array 'all', or cell with a single subject denomination e.g. {'DOE'} or {'042'}
input.dates         = 'all'; % {'20231031'}; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}

% C) GENERAL Options.
opt = struct();
    opt.RetrieveEvents      = true;     % Retrieve event log.
        opt.useexe          = false;    % Eventually, only option for Deuteron recordings (TODO)
        opt.usepar          = true;     % temporarily use of .par files from Juan's behavior paradigm
    opt.GetMotionSensors    = false;    % Retrieve data from motion sensors in Deuteron. (TODO)
    opt.kilosort            = true;     % Call to kilosort processing. !! NEEDS configfile saved under 'studyName\analysisCode\'
        opt.spkTh           = -2;       % Default: -4.5.
        opt.KSchanMapFile   = '';       % Empty to use simple non-mapped linear probe, or e.g.'chanMapPoly3Deut', 'chanMapATLASTri'. 
    opt.bombcell            = true;     % Run bombcell on the KS output, previously to manual curation
        opt.rerun           = true;     % To overwrite previous runs of BombCell
        opt.nRawSpikesToExtract = 1000; % Parameter for bombcell run
    opt.phy                 = false;     % Calls phy for manual inspection or curation. !! It PUTS MATLAB on HOLD!
    opt.FieldTrip           = true;     % Creation of FieldTrip ready .mat file.
    opt.parsetrial          = false;    % Define and split data into trials

% D) INCLUDE README.TXT
% It contains details about the project. File can also be modified later.
readmecontent = ["Study name: StudyName"                                , ...
                 "Readme date: DD/MM/YYYY"                              , ...
                 "Person (1) responsible for data repository: X "       , ...
                 "Person(s) responsible for study: Y, Z "               , ...
                 "Hardware used: Deuteron/Intan "                       , ...
                 "Related Publication(s): None."                        , ...
                 "Short description of study: Bird stuff "                   ];

%% RUN pipelines
%% Start with Project folder system preparation. Commonly to run only once,
% before any data exists, since the folder for the raw data is created
% here. If it already exists, nothing will change.

NGL00_Prep

%% Continue with the main script, which most of times will locate sessions, determine
% formats, extarct EventCodes and Motion data, convert to Kilosort and FieldTtrip formats, 
% and perform Kilosort and BombCell analysis. 
% Additionally it could launch Phy for manual curation, use EventCodes to
% trial parse the data and others.

NGL01_Main

%% The following script is still on the works. Most likely will include
% statistical treatments and plots, and NO FURTHER processing (currently
% some is done)

% NGL02_postPhy

%% Other scripts to come.

% NGLXX_something
