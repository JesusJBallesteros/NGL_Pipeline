%% 0) READ. 
% SET AND RUN
% Input file where there is no access to any of the running code, making
% this the only file that needs to be modified, and that could call all
% pipelines as a sequence of easily swichable runs by simply commenting 
% lines. Could be ran line-by-line (F9) or all at once (F5).
%
% Once the folder system is created, should this file be located INSIDE?
% i.e. at '...\studyName\data\analysis' ??

%% 1) SET.
clear all
% A) SYSTEM
% RECOMMENDED to add the toolbox folder to MATLAB folder system, but not necessary.
datadrive   = 'F';                   % The LETTER of the drive where the data structure is/will be created.
studyname   = 'projectName';  % Name of the study to be used (main folder for the data)
toolbox     = 'C:\Code\ephys-data-pipeline'; % Absolute path to the toolbox.

% B) SUBJECTS AND SESSIONS
% To run the script on all subjects and sessions, or as session-to-session process.
subjects    = {'123'}; % 'all';      % char array 'all', or cell with a single subject denomination e.g. {'DOE'} or {'042'}.
dates       = {'YYYYMMDD'}; % 'all'; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}.

% C) OPTIONS.
opt = struct();
    % General options for NGL01_Main
    opt.numChannels             = 32;       % For now, explicit 32. (Deprecate?)
    opt.CAR                     = 1;        % Default: 1. To remove fast-ample transients and other noise for .bin file. If == 2 also CAR for lowpass (not recommended)
    opt.linefilter              = 0;        % If > 0, filter noise at given value +-2 (Hz)
    opt.bin                     = true;     % Create a .bin file with the high-pass data, input to Kilosort for spike sorting.
        opt.highpass            = 400;      % Low frequency boundary (Hz)
    opt.FieldTrip               = true;     % Create a FieldTrip ready .mat file with the low-pass data, either continuous, trial-parsed or both. 
        opt.lowpass             = 150;      % High frequency boundary (Hz)
    opt.GetMotionSensors        = true;     % Retrieve data from motion sensors in Deuteron. NEEDS IMPROVEMENT. Not yet for Intan's accelerometer
    opt.RetrieveEvents          = true;     % Retrieve events.
        opt.alignto             = {'itiOn', 'stimOn1'}; % cell array e.g. {'itiOn', 'rwd'}. 'itiON' should be the very minimum to align to.
        opt.trEvents            = {'tr1'};  % Explicit Out-trial events, i.e. events at ITI like treatments, tutors, block or phase changes.
    opt.kilosort                = 4;        % 2/4 for KS2/KS4 !! NEEDS corresponding config files saved under 'studyName\analysisCode\'
        opt.KSchanMapFile       = 'chanMapE32-S2_DeutSN11.mat';  % Empty '' to use linearly increasing array. Or e.g. 'chanMapXXX.mat' for maps saved under 'studyName\analysisCode\'
        opt.spkTh               = -4.5;     % For KS2, usually a single value. If array, it runs iterations with threshold -X, -Y ... -Z.
    opt.bombcell                = false;    % Run bombcell on the KS output, as previous step to manual curation. !! NEEDS revision after KS updates.
%          opt.rerun              = false;    % To overwrite previous runs of BombCell.
    opt.phy                     = false;    % Open phy for manual inspection or curation. !! It puts MATLAB on HOLD!

% D) README.TXT
% It contains details about the project. File can also be modified later.
readmecontent = ["Study name: Generic Readme project", ...
                 "Readme date: 06/06/0006"                          , ...
                 "Person (1) responsible for data repository: J DOE", ...
                 "Person(s) responsible for study: Mr. White"    , ...
                 "Hardware used: Deuteron. 2x ATLAS E-32-S2"                  , ...
                 "Related Publication(s): None."                    , ...
                 "Short description of study: Template" ];

%% 2) RUN.
%% 2.0 Start with project folder system preparation. Commonly to be ran only ONCE,
% before any data exists, since the folder for the raw data is created
% here. If it already exists, nothing will change.
NGL00_Prep

%% 2.1 Continue with the Main script, which locate sessions, determine formats, extract
% EventCodes and Motion data, convert to Kilosort and FieldTtrip formats, 
% and perform Kilosort automatic sorting. 
% Additionally it can launch Phy for manual curation after each sessions, or first
% run Bombcell to semi-automatize this porcess (only once appropiate
% parameters are known) and then launch Phy.
NGL01_Main

%% 2.2 Proceed with post-Phy processing. Once data is curated.
% Uses events and trial definitions obtained before to trial-parse the spike or
% LFP data, creating the variables into the lab standard.
NGL02_postPhy

%% 2.3 Plotting. (Some plotting happens before this. Rearrange?)
% Some basic plots are provided for exploratory-descriptive plotting, 
% either for the day-to-day data check or for the whole of sessions
% plotting, to obtain examples of clusters, effects, etc.)
% Othert elaborated or dedicated plots could be added on a personal basis,
% or implemented as default if decided as standard.
NGL03_plotting 

%% 2.etc Statistics?
% etc