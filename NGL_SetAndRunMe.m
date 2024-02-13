%% SET AND RUN
% Input file where there is no access to any of the running code, making
% this the only file that needs to be modified, and that could call all
% pipelines as a sequence of easily swichable runs by simply commenting 
% lines. Could be ran line-by-line (F9) or all at once (F5).
%
% Once the folder system is created, should this file be located INSIDE?
% i.e. at '...\studyName\data\analysis' ??

%% 1) SET
% A) SYSTEM 
% Specify drive where data is located, Project Name and the toolbox folder.
% RECOMMENDED to add the toolbox folder to MATLAB folder system, but not necessary.
datadrive   = 'F';                   % Define the LETTER of the drive where the data structure will be created or already exists.
studyname   = 'MotionSensing';       % Name the Study or Project to be used.
toolbox     = 'C:\Code\ephys-data-pipeline';% Absolute address to the toolbox.

% B) SUBJECTS AND SESSIONS
% To run the script on all subjects and sessions, or as session-to-session process.
subjects    = {'079'};      % char array 'all', or cell with a single subject denomination e.g. {'DOE'} or {'042'}.
dates       = {'20240123'}; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}.

% C) OPTIONS.
opt = struct();
    opt.cooking = false;    % Temporary option to run or not things under development
    % Data Extraction and Pre-processing
    opt.bin                     = true;    % Create a .bin file with the high-pass data, usually to be passed to Kilosort for spike sorting.
    opt.FieldTrip               = true;    % Create a FieldTrip ready .mat file with the low-pass data, either continuous, trial-parsed or both. 
    opt.GetMotionSensors        = true;    % Retrieve data from motion sensors in Deuteron. NEEDS IMPROVEMENT on head direction interpretation.
    opt.RetrieveEvents          = true;    % Retrieve event log. So far, for Deuteron.
        % opt.eventdef            = [];      % To pass non-standard event descriptions. If missing, use NGL standad.

        % opt.parsetrial          = false;   % Define and split data into trials. Prob to discontinue as will be assumed true when 'RetrieveEvents' = true
        % opt.useexe              = false;   % Use Deuteron's executable software. Prob to be discontinued.
        % opt.usepar              = false;   % temporarily use '_par' files from Juan's behavior paradigm to extract events. Prob to be discontinued.
    
    % Data Processing
    opt.kilosort                = true;     % Call to kilosort processing. !! NEEDS configfile saved under 'studyName\analysisCode\'
        opt.spkTh               = -4.5;     % A single value. Default: -4.5. Possible to introduce an array of negative values to run several thresholds. Results are saved in subfolders and bin file needs to be copied manually inside to check sorter results.
        opt.KSchanMapFile       = '';       % Empty to use simple, non-mapped, linear array. Or e.g.'chanMapPoly3Deut', 'chanMapATLASTri' for custom maps saved under 'studyName\analysisCode\'
    opt.bombcell                = false;    % Run bombcell on the KS output, as previous step to manual curation. TODO: go over several KS outputs if existing.
%         opt.rerun               = false;    % To overwrite previous runs of BombCell.
    opt.phy                     = false;    % Open phy for manual inspection or curation. !! It puts MATLAB on HOLD!

%    % Post-processing (after manual curation) TODO
%     opt.postPhy                  = false;     % Would habilitate the postPhy processes. Prob not important.
%     opt.KS2spkmat                = false;   % We need the KS output, once curated, to be read and saved as MATLAB structure. This could be analyzed independently.
%     opt.spkmat2FT                = false;   % We can add the spike data to the Fieldtrip LFP data for combined analysis.
%  

% D) README.TXT
% It contains details about the project. File can also be modified later.
readmecontent = ["Study name: MotionSensing"                 , ...
                 "Readme date: DD/MM/YYYY"                          , ...
                 "Person (1) responsible for data repository: X "   , ...
                 "Person(s) responsible for study: Y, Z "           , ...
                 "Hardware used: Deuteron "                         , ...
                 "Related Publication(s): None."                    , ...
                 "Short description of study: Bird stuff "          ];

%% 2) RUN
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
