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
studyname   = 'Pre_SocialLearning';  % Name of the study to be used (main folder for the data)
toolbox     = 'C:\Code\ephys-data-pipeline'; % Absolute path to the toolbox.

% B) SUBJECTS AND SESSIONS
% To run the script on all subjects and sessions, or as session-to-session process.
subjects    = {'913'};      % char array 'all', or cell with a single subject denomination e.g. {'DOE'} or {'042'}.
dates       = {'20240219'}; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}.

% C) OPTIONS.
opt = struct();
    % opt.cooking = false;    % Temporary option to run or not things under development
    opt.numChannels             = 32;       % For now, explicit 32 if not SpikeLog-64C was used (Deuteron). INTAN: comment.
    opt.bin                     = true;     % Create a .bin file with the high-pass data, to be passed to Kilosort for spike sorting.
        opt.highpass            = [450 7000]; % Give as [low high] frequency values.
        opt.CAR                 = 1;        % Default: 1. Common Average Referencing (median) to remove fast-ample transients and other noise. 
    opt.FieldTrip               = false;     % Create a FieldTrip ready .mat file with the low-pass data, either continuous, trial-parsed or both. 
        opt.lowpass             = [0 250];  % Give as [low high] frequency values.
    opt.GetMotionSensors        = false;    % Retrieve data from motion sensors in Deuteron. NEEDS IMPROVEMENT on head direction interpretation.
    opt.RetrieveEvents          = false;    % Retrieve event log. If not further options defaulted to Deuteron txt log extraction.
        % opt.useexe              = false;   % Use Deuteron's executable software. Prob to be discontinued.
        % opt.usepar              = false;   % temporarily use '_par' files from Juan's behavior paradigm to extract events. Prob to be discontinued.
    % opt.parsetrial              = false;   % Define and split data into trials. Prob to discontinue as will be assumed true when 'RetrieveEvents' = true
        % opt.eventdef            = [];      % To pass non-standard event descriptions. If missing, use NGL standad.
    opt.kilosort                = 4;        % Kilosort processing. == 2 for KS2, == 4 for KS4 !! KS2 NEEDS configfile saved under 'studyName\analysisCode\'
        opt.KSchanMapFile       = 'chanMapE32-S2_linearized_DeutSN11.mat';  % Empty '' to use non-mapped, linear array. Or e.g.'chanMapXXX.mat' for custom maps saved under 'studyName\analysisCode\'
        opt.spkTh               = -4.5;     % Usually a single value. If multiple [-X -Y ... -Z], cycle runs with thresholds -X, -Y ... -Z each.
    opt.bombcell                = false;     % Run bombcell on the KS output, as previous step to manual curation. TODO: go over several KS outputs if existing.
         opt.rerun              = false;    % To overwrite previous runs of BombCell.
    opt.phy                     = false;    % Open phy for manual inspection or curation. !! It puts MATLAB on HOLD! Needs bin file in same folder.
                        
%    % Post-processing (after manual curation) TODO
%     opt.postPhy                  = false;     % Would habilitate the postPhy processes. Prob not important.
%     opt.KS2spkmat                = false;   % We need the KS output, once curated, to be read and saved as MATLAB structure. This could be analyzed independently.
%     opt.spkmat2FT                = false;   % We can add the spike data to the Fieldtrip LFP data for combined analysis.

% D) README.TXT
% It contains details about the project. File can also be modified later.
readmecontent = ["Study name: Preparation for social learning paradigms", ...
                 "Readme date: 22/02/2024"                          , ...
                 "Person (1) responsible for data repository: Jesus", ...
                 "Person(s) responsible for study: Juan, Jesus "    , ...
                 "Hardware used: Deuteron. E32-S2"                  , ...
                 "Related Publication(s): None."                    , ...
                 "Short description of study: Autoshaping with social learning." ];

%% 2) RUN.
%% 2.0 Start with project folder system preparation. Commonly to be ran only ONCE,
% before any data exists, since the folder for the raw data is created
% here. If it already exists, nothing will change.
NGL00_Prep

%% 2.1 Continue with the Main script, which locate sessions, determine formats, extract
% EventCodes and Motion data, convert to Kilosort and FieldTtrip formats, 
% and perform Kilosort and BombCell analysis. 
% Additionally it could launch Phy for manual curation, use EventCodes to
% trial parse the data and others.
NGL01_Main

%% 2.2 Still on the works. Most likely will include
% statistical treatments and plots, and NO FURTHER processing (currently
% some is done)
% NGL02_postPhy

%% 2.X Other scripts to come.

% NGLXX_something
