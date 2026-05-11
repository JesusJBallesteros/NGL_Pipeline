%% NGL_SetAndRunMe  — User Configuration and Run Script (TEMPLATE)
%
% PURPOSE:
%   The single file a user must edit to drive the entire NGL pipeline.
%   Copy this template from the toolbox root into your project's
%   analysisCode/ folder and edit it there. Never run it from the toolbox
%   root itself.
%
% USAGE:
%   Run section-by-section (F9) or all at once (F5) from MATLAB.
%   All subsequent pipeline scripts (NGL01_Main, NGL02_postPhy, etc.) are
%   called from within this file — they should never be opened directly.
%
% SECTIONS:
%   1) PREPARE   — set study metadata, drive and project name, run NGL00_Prep
%   2) SET       — choose subjects/sessions and configure opt struct
%   3) RUN       — call NGL01_Main, NGL02_postPhy, and optional stages
%
% REQUIRED CONFIG FILES (place in analysisCode/ alongside this script):
%   NGL_machineConfig.m   — machine-specific Python env and toolbox paths
%   eventDefinitions.m    — project event code definitions
%   conditions_script.m   — condition grouping logic
%   chanMapXXX.mat        — custom Kilosort channel map (if applicable)
%
% NOTES:
%   - opt fields not listed here receive safe defaults from default_opt.m.
%   - See wiki_NGL01_pipeline.md for a full opt field reference.
%   - NGL_machineConfig.m is machine-specific; do not commit it to git.
%
% Last modified 06.05.2026 (Jesus)

%% 1) PREPARE.
clear all
% A) README.TXT
readmecontent = ["Study name: DefaultName", ...
                 "Readme date: 29/08/1997"                          , ...
                 "Person (1) responsible for data repository: Main researcher", ...
                 "Person(s) responsible for study: Main R and associates"    , ...
                 "Hardware used: ACQ SYS. Probe type. Arena/Box"                  , ...
                 "Related Publication(s): None."                    , ...
                 "Short description of study: A more elaborated description of project paradigm, goals, etc ." ];

% B) DATA LOCATION
% ADD the toolbox folder to MATLAB folder system !!
datadrive   = 'D';                   % The LETTER of the drive where the data structure is/will be created.
studyname   = 'ProjectName';  % Name of the study to be used (main folder for the data)

% Before any data exists, the folder for the raw data is created here.
% If it already exists, nothing will change.
NGL00_Prep

%% 2) SET.
% A) SUBJECTS AND SESSIONS
% To run the script on all subjects and sessions, or as session-to-session process.
subjects    = {'XXX'}; % char array 'all', or cell with a single subject denomination e.g. {'DOE'} or {'042'}.
dates       = {'YYYYMMDD', 'YYYYMMDD'}; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}.

% A2) MULTI-AREA (optional).
% Uncomment and fill in if your probe spans more than one brain area.
% One label per kcoords group in the chanMap (kcoords==1 → Areas{1}, etc.).
% Repeated labels indicate shanks from the same area (processed together).
% Kilosort will run once per unique area; results go to preprocessing\<session>\<Area>\.
% Leave commented out for standard single-area behaviour.
%
%   input.Areas = {'NCL', 'NCL', 'STR'};   % 2 NCL shanks (kcoords 1-2) + 1 STR shank (kcoords 3)

% B) OPTIONS.
opt = struct();
    % NECESSARY options for NGL01_Main
    opt.numChannels         = 32;      % For now, explicit 32 if not SpikeLog-64C was used (Deuteron). INTAN: comment.
    opt.KSchanMapFile       = 'chanMap_XXX.mat';  % Empty '' to use non-mapped, linear array. Or e.g.'chanMapXXX.mat' for custom maps saved under 'studyName\analysisCode\'
    opt.RetrieveEvents      = true;    % Retrieve event log. If not further options defaulted to Deuteron txt log extraction.
        opt.alignto         = {'itiOn', 'stimOn1', 'rwd'};  % single char array e.g. 'itiOn', or cell array e.g. {'itiOn', 'rwd'}. 'itiON' should be the very least to align to.
    opt.GetMotionSensors    = false;   % Retrieve data from motion sensors in Deuteron. NEEDS IMPROVEMENT on head direction interpretation.
    opt.FieldTrip           = false;   % Create a FieldTrip ready .mat file with the low-pass data, either continuous, trial-parsed or both. 
    opt.bombcell            = true;    % Run bombcell on the KS output. Previous step to manual curation.
    opt.lowpass             = 10000;   % If < 9500, high boundary frequency value for low-pass.

    % Change only with good reasons.
    opt.addtime             = 0;       % Expands the trial definition around start/end by X ms in both directions. 
    opt.trEvents            = {};      % Add inter-trial events, if any, to delimit e.g. block changes 
    opt.phy                 = false;   % Open phy for manual inspection or curation. !! It puts MATLAB on HOLD! Needs bin file in same folder.
        % !! Realize that manual curation via PHY must be PERFORMED, to use Post-Phy scripts.
        % But it does NOT need to be IMMEDIATELY after KS-BC automatic job.

    % opt.noise               = [];
    % opt.doNWB               = false;   % TESTING INTAN-NEUROCONV (python) with a Matlab wrapping for no python-user interaction
    % opt.CAR                 = 0;       % If not 0, removes fast-ample transients and other noise. (KS4 should do this)
    % opt.linefilter          = 0;       % If not 0, filter line noise at given value +-2 (Hz)
    % opt.lowpassFT           = 250;     % Give as high boundary frequency value for FT.
                        
%% 3) RUN.
% 3.1 Continue with the Main script, which locate sessions, determine formats, extract
% EventCodes and Motion data, convert to Kilosort and FieldTtrip formats, 
% and perform Kilosort automatic sorting. 
% Additionally it can launch Phy for manual curation after each sessions, or first
% run Bombcell to semi-automatize this porcess (only once appropiate
% parameters are known) and then launch Phy.
NGL01_Main

%% 3.2 Proceed with post-Phy processing. Once data is curated.
% Includes steps towards spike/trial sorting of the curated data. Uses
% events and trial definitions obtained before to trial-parse the spike or
% LFP data, creating the variables into the lab standard.

% General options for NGL02_postPhy
cd(input.analysisCode)
postPhy_param();

% *IMPORTANT*: phy2 must have been run beforehand, so a key file exists to
% extract information from
cd(toolbox)
NGL02_postPhy

%% 3.3 Plotting.
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

%% 3.4 Aggregating. (IN PROGRESS)
% Get data from all specified animals and sessions and aggregate them into
% single variables.

NGL04_aggregate % So far, only living inside 'SocialLearning'

%% 2.XX More...
% NGLXX_something
