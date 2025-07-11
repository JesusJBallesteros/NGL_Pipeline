%% 0) READ. 
% SET AND RUN
% Input file where there is no access to any of the running code, making
% this the only file that needs to be modified, and that could call all
% pipelines as a sequence of easily swichable runs by simply commenting 
% lines. Could be ran line-by-line (F9) or all at once (F5).

%% 1) SET.
clear all
% A) SYSTEM
% RECOMMENDED to add the toolbox folder to MATLAB folder system, but not necessary.
datadrive   = 'D';                   % The LETTER of the drive where the data structure is/will be created.
studyname   = 'Default';  % Name of the study to be used (main folder for the data)
toolbox     = 'C:\Code\ephys-data-pipeline'; % Absolute path to the toolbox.

% B) SUBJECTS AND SESSIONS
% To run the script on all subjects and sessions, or as session-to-session process.
subjects    = {'666'}; % char array 'all', or cell with a single subject denomination e.g. {'DOE'} or {'042'}.
dates       = {'19970829'}; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}.

% C) OPTIONS.
opt = struct();
    % General options for NGL01_Main
    opt.numChannels             = 32;       % For now, explicit 32 if not SpikeLog-64C was used (Deuteron). INTAN: comment.
    opt.CAR                     = 1;        % Default: 1. CAR to remove fast-ample transients and other noise for .bin file. If == 2 also CAR for lowpass (not recommended)
    opt.linefilter              = 0;        % If not 0, filter line noise at given value +-2 (Hz)
    opt.doNWB                   = true;     % TESTING INTAN-NEUROCONV (python) with a Matlab wrapping for no python-user interaction
    opt.bin                     = true;     % Create a .bin file with the high-pass data, to be passed to Kilosort for spike sorting.
        opt.highpass            = 400;      % Give as low boundary frequency value. (High boundary is set at recording time)
    opt.FieldTrip               = true;     % Create a FieldTrip ready .mat file with the low-pass data, either continuous, trial-parsed or both. 
        opt.lowpass             = 250;      % Give as high boundary frequency value.
    opt.GetMotionSensors        = true;    % Retrieve data from motion sensors in Deuteron. NEEDS IMPROVEMENT on head direction interpretation.
    opt.RetrieveEvents          = true;     % Retrieve event log. If not further options defaulted to Deuteron txt log extraction.
        opt.alignto             = {'itiOn', 'stimOn1', 'rwd'};  % single char array e.g. 'itiOn', or cell array e.g. {'itiOn', 'rwd'}. 'itiON' should be the very least to align to.
        opt.trEvents            = {'na3'};  % Event definition of 'special events' i.e. events at the ITI like treatments, tutors, etc...
        opt.addtime             = 1500;     % Expands the trial definition start/end by X ms in both directions. 
    opt.kilosort                = 4;        % Kilosort processing. == 2 for KS2, == 4 for KS4 !! KS2 NEEDS configfile saved under 'studyName\analysisCode\'
        opt.KSchanMapFile       = 'chanMapE32-S2_DeutSN11.mat';  % Empty '' to use non-mapped, linear array. Or e.g.'chanMapXXX.mat' for custom maps saved under 'studyName\analysisCode\'
        opt.spkTh               = -6;       % Only for KS2. Usually a single value. If multiple [-X -Y ... -Z], cycle runs with thresholds -X, -Y ... -Z each.
    opt.bombcell                = false;    % Run bombcell on the KS output. =2 (KS2) or =4 (KS4). Previous step to manual curation.
         opt.rerun              = false;    % To overwrite previous runs of BombCell.
    opt.phy                     = false;    % Open phy for manual inspection or curation. !! It puts MATLAB on HOLD! Needs bin file in same folder.
                        
    % General options for NGL02_postPhy
    postPhy_param();

% D) README.TXT
% It contains details about the project. File can also be modified later.
readmecontent = ["Study name: DefaultName", ...
                 "Readme date: 29/08/1997"                          , ...
                 "Person (1) responsible for data repository: Main researcher", ...
                 "Person(s) responsible for study: Main R and associates"    , ...
                 "Hardware used: ACQ SYS. Probe type. Arena/Box"                  , ...
                 "Related Publication(s): None."                    , ...
                 "Short description of study: A more elaborated description of project paradigm, goals, etc ." ];

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

% TODO. separate data from different ports at this level to
%       effective CAR use on different brain regions
NGL01_Main

%% 2.2 Proceed with post-Phy processing. Once data is curated.
% Includes steps towards spike/trial sorting of the curated data. Uses
% events and trial definitions obtained before to trial-parse the spike or
% LFP data, creating the variables into the lab standard.

% TODO
NGL02_postPhy

%% 2.3 Plotting.
% Having all necessary variables ('neurons', 'events', 'conditions',
% 'spike', 'trialdef', etc...) proceed to plot data. 
%
% Some basic plots are provided for exploratory-descriptive plotting, 
% either for the day-to-day data check or for the whole of sessions
% plotting, to obtain examples of clusters, effects, etc.)
% Othert elaborated or dedicated plots could be added on a personal basis,
% or implemented as default if decided as standard.

% TODO
% NGL03_plotting01
% NGL03_plotting02

%% 2.4 Aggregating. (IN PROGRESS)
% Get data from all specified animals and sessions and aggregate them into
% single variables.

NGL04_aggregate % So far, only living inside 'SocialLearning'

%% 2.XX More...
% NGLXX_something
