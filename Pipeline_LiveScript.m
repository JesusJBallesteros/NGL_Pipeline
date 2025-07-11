%[text] %[text:anchor:T_5E89BFB6] # Ephys Data Pipeline, LiveScript version
%[text] General explanation and Instructions go here
%[text:tableOfContents]{"heading":"**Table of Contents**"}
%%
%[text] %[text:anchor:T_5782] # Main inputs
%[text] %[text:anchor:H_320a] ## Paths to toolbox and data
%[text] %[text:anchor:H_6d7b] Absolute path to the toolbox and project folders.
toolbox = "C:\Code\ephys-data-pipeline"; %[control:filebrowser:6d69]{"position":[11,40]}
Project = "D:\S3_ExperimentArena"; %[control:filebrowser:8d93]{"position":[11,34]}
%[text] %[text:anchor:H_9913] ## Data selection
%[text] Input a string 'all' to use all valid directories at project, or explicit a subset of them
subjects = 'all'; % Input as 'all' or {'DOE', ..., '666'} %[control:editfield:77cf]{"position":[12,17]}
dates    = 'all'; % Input as 'all' or {'YYYYMMDD', ..., 'YYYYMMDD'} %[control:editfield:8f23]{"position":[12,17]}
%[text] Get into toolbox's given path and set default parameters.
cd(toolbox); addpath(toolbox);
default_opt;
%[text] %[text:anchor:T_8f81] # User options
%[text] %[text:anchor:H_4A35002E] Defaults will be used if no changes are made. See 'default\_opt'. Alternatively, you can modify the  options below, or add your own.
%[text] Number of channels:
opt.numChannels =32; % channels %[control:slider:0798]{"position":[18,20]}
%[text] Common Average Reference (remove fast-ample transients and other noise):
opt.CAR = true; % The method actually uses the median, not the average %[control:checkbox:5e20]{"position":[11,15]}
%[text] Band-stop Filter at value ($\\pm${"editStyle":"visual"}2 Hz):
opt.linefilter = 0; % Hz %[control:dropdown:138b]{"position":[18,19]}
%[text] Create .bin file for Kilosort to spike sorting:
opt.bin = true; %[control:checkbox:1ac4]{"position":[11,15]}
%[text] High-pass filter's lower boundary frequency:
opt.highpass = 400; % Hz %[control:slider:0bf3]{"position":[16,19]}
%[text]  Create .mat file with FieldTrip format (continuous, trial-parsed or both): 
opt.FieldTrip = false; %[control:checkbox:835c]{"position":[17,22]}
%[text] Low-pass filter's higher boundary frequency:
opt.lowpass = 250; % Hz %[control:slider:65af]{"position":[15,18]}
%[text] Retrieve data from motion sensors: (Deuteron only, so far)
opt.GetMotionSensors = false; %[control:checkbox:9b7b]{"position":[24,29]}
%[text] Retrieve Event Codes and timestamps:
opt.RetrieveEvents = true; %[control:checkbox:5fc5]{"position":[22,26]}
%[text] Event names to align data to:
opt.alignto = {'itiOn','stimOn1', 'rwd'}; % keep 'itiON' %[control:editfield:188c]{"position":[24,40]}
%[text] If any, ITI Event names to locate (e.g stimulus type, block change, tutor in/out):
opt.trEvents = {'na1', 'na2'}; %[control:editfield:954b]{"position":[17,29]}
%[text] Miliseconds of data to retrieve from start/end of each trial, in both directions:
opt.addtime = 1500; % miliseconds %[control:slider:3fad]{"position":[15,19]}
%[text] Kilosort processing via v.2 or v.4. A **configuration file is needed** at '...\\analysisCode\\'
opt.kilosort = 4; % to DEPRECATE and leave 4 only  %[control:dropdown:760b]{"position":[16,17]}
%[text] Channel map file for sorting. The **file must be at** '...\\analysisCode\\' **with name** 'chanMapXXX.mat'
opt.KSchanMapFile = "chanMapXXX.mat"; % empty for a non-mapped, linear array %[control:editfield:061a]{"position":[21,37]}
%[text] Threshold for spike detection (only for KS2):
opt.spkTh = -6; % (std dev) % to DEPRECATE %[control:slider:13e2]{"position":[13,15]}
%[text] Run BombCell automatic curation:
opt.bombcell = 0; %[control:dropdown:4406]{"position":[16,17]}
%[text] Open Phy automatically, for manual curation. 
opt.phy = false; % Puts MatLab on HOLD until Phy is closed. %[control:checkbox:68a4]{"position":[11,16]}
%[text] Jump directly to plot functions:
opt.jump2plot = false; % if preprocessing was already done %[control:checkbox:4be0]{"position":[17,22]}
%[text] Create the to-plot list:
opt.plot.Raster_trialbytrial = false; %[control:checkbox:6bbb]{"position":[32,37]}
opt.plot.Raster_Aligned      = false; %[control:checkbox:4cf3]{"position":[32,37]}
opt.plot.Normalized_firerate = false; %[control:checkbox:6711]{"position":[32,37]}
opt.plot.YYY = false; %[control:checkbox:24ce]{"position":[16,21]}
opt.plot.ZZZ = false; %[control:checkbox:3610]{"position":[16,21]}
%[text] Some <u>General Parameters</u> for plots are set by default as 'param' in 'default\_opt'. Add any specific parameters here:
param.trial2plot = 'allInitiated'; % 'correct', 'incorrect', 'omission', 'allInitiated'
param.FS2plot   = true; 
param.nBlocks   = 8;
param.interval  = [-2000 10000];
param.baseline  = 2000;
%[text] %[text:anchor:H_790c] ### <u>**IMPORTANT**</u>: Add Project-specific configurations, maps, parameters.
%[text] All should be located inside '...\\analysisCode\\' folder.
%[text] <u>**Once all is ready, uncheck**</u> the box <u>**and click**</u> Run Section.
StopGo = true;   %[control:checkbox:88e9]{"position":[10,14]} %[control:button:97e9]{"position":[16,17]}
%[text] %[text:anchor:H_8179] Summary of options:
opt %[output:50f8331b]
%%
%[text] %[text:anchor:T_5663] # <u>**00**</u>. Check folder system, inputs and dependencies.
%[text] No user input is needed. <u>**Click to Continue.**</u>
  %[control:button:8b0b]{"position":[1,2]}
%[text] %[text:anchor:H_6a33] #### Folder syster preparation.
%[text] Get the Drive where the data structure is/will be created and the desired Project Name.
[datadrive, studyname] = fileparts(Project);
%[text] Prepares the project folder system if not existing, and copies all default scripts and configurations to the appropiate folder. If everything is already set from a previous run, nothing will change.
NGL00_Prep %[text:anchor:H_666E4CBA] %[output:5fd23422]
%[text] %[text:anchor:H_1dab] #### Set Default global parameters and dependencies.
%[text] Inputs are collected and structured for convenience.
input = struct('datadrive', datadrive, 'studyName', studyname, 'toolbox', toolbox, ... % force char array
               'subjects', [], 'dates', []); % do NOT force char array
input.dates    = dates;    % place as it comes
input.subjects = subjects; % place as it comes
%[text] Add other default inputs and set necessary dependencies.
input = set_default(input, opt);
%[text] %[text:anchor:H_5447] #### Locate requested sessions.
input.sessions = findSessions(input);
%[text] Summary of data to be processed:
input.sessions
%%
%[text] %[text:anchor:T_78a6] # <u>**01**</u>. Find and run the data preprocessing.
%[text] No user input is needed. <u>**Click to Continue.**</u>
  %[control:button:3bbe]{"position":[1,2]}
%[text] Continue with locating sessions, determine formats, extract EventCodes and Motion data, convert to Kilosort and FieldTtrip formats, performs Kilosort automatic sorting. Additionally it can launch Phy for manual curation after each sessions, or first run Bombcell to semi-automatize this porcess (only once appropiate parameters are known) and then launch Phy.
if StopGo, error('Did you put all your customized files in its place? Did you check the box?'); end %[output:296e2cb5]
%[text] If all checks passed, then we start preprocessing in a subject and session basis.
if ~opt.jump2plot % If not expicitly skipped
 for x = 1:input.nsubjects % Go over subjects
 for y = 1:input.sessions(x).nsessions % Go over sessions
%[text] Prepare for session preprocessing.
 input.run = [x y]; 
 [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);
 input.sessions(input.run(1)).info.fileformat
%[text] Direct script towards the appropiate pipeline based on data format:
 switch input.sessions(input.run(1)).info.fileformat
%[text] %[text:anchor:H_7efe] ###  <u>Deuteron Pipeline</u>
    case {'DT2', 'DF1'} 
%[text]             Aggregate a few parameters for convenience:
    opt.myFiles = input.sessions(input.run(1)).info.files; % Prob can go out
    opt.ext     = input.sessions(input.run(1)).info.fileformat; % Prob can go out            
    opt.numberOfAdcBits   = input.sessions(input.run(1)).info.numADCBits;
    opt.sampleRate        = input.sessions(input.run(1)).info.amplifier_sample_rate;
    opt.voltageResolution = input.sessions(input.run(1)).info.voltageRes;
    opt.offset            = 2^(opt.numberOfAdcBits-1);
%[text] %[text:anchor:H_54d7] ####           Extract Events (General method)
    [events, trialdef, EventRecord] = EventProcess(opt); % Needs trialdef existence check inside)
%[text] %[text:anchor:H_5689] ####            Convert Deuteron highpass data to .bin format, if requested.
    if opt.bin, disp('Converting Deuteron files to .bin format.')
       Deuteron2Kilosort(opt);
    end
%[text] %[text:anchor:H_4171] ####           Convert Deuteron lowpass data to FieldTrip format, if requested.
    if opt.FieldTrip, disp('Converting Deuteron files into a pseudo-FT file.')
       FT_data = Deuteron2Fieldtrip(opt);
%[text] %[text:anchor:H_0632] ####           Giving proper FieldTrip format (General method)
       disp('Giving proper FieldTrip format.')
       MAT2FieldTrip(FT_data, opt, trialdef, 1); 
    end
%[text] %[text:anchor:H_7fdc] ####           Get Motion Data from Deuteron data, if requested.
    if opt.GetMotionSensors, disp('Extracting Motion Sensor data from Deuteron.')
       Deuteron_GetMotionSensors(opt);
    end
%[text] %[text:anchor:H_5e8f] ###  <u>INTAN Pipeline</u>
    case {'fileperch', 'filepertype'}
%[text] %[text:anchor:H_2ac1] ####           Get INTAN metadata
    input.sessions(input.run(1)) = findSetting(input.sessions(input.run(1)));
%[text] %[text:anchor:H_17bd] ####           Extract Events (General method)
    [events, trialdef, EventRecord] = EventProcess(opt); % Needs trialdef existence check inside)
%[text] %[text:anchor:H_3590] ####           NWB file format creation (INTAN specific)
    if input.useNWB, intan2NWB_wrapper(input, opt); end
%[text] %[text:anchor:H_528f] ####           INTAN to Kilosort method (INTAN specific)
    if opt.bin && ~isfile(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']))
       Intan2Kilosort_wrapper(input.sessions(input.run(1)), opt);
    end
%[text] %[text:anchor:H_9e9e] ####           Convert INTAN lowpass data to FieldTrip format, if requested.
    if opt.FieldTrip && ~isfile(fullfile(opt.analysis,[opt.SavFileName '_FTcont.mat']))
       FT_data = intan2MAT_wrapper(input.sessions(input.run(1)), opt);
%[text] %[text:anchor:H_290e] ####           Giving proper FieldTrip format (General method)
       disp('Giving proper FieldTrip format.')
       MAT2FieldTrip(FT_data, opt, trialdef, 1);
    end
%[text] %[text:anchor:H_133e] ####           Get Motion Data from INTAN data, if requested. (TODO)
    if opt.GetMotionSensors, disp('Extracting Motion Sensor data from INTAN.')
       % not_a_function_yet(opt);
    end
%[text] %[text:anchor:H_1e13] ###  <u>FieldTrip Pipeline</u>
    case {'FieldTrip'}
%[text] %[text:anchor:H_317f] ####            Extract Events (General method)
        [events, trialdef, EventRecord] = EventProcess(opt);
        disp('Session skipped because continuous FT file was found'), continue
%[text]             If the format is not any of the previous cases, throw a warning and progress to next  
 otherwise, warning('Something went wrong during format verification. Skipping Session'), continue
 end % format switch
%[text] %[text:anchor:H_8b54] ### <u>Kilosort Processing</u>
%[text] Any Kilosort instance will run withour GUI. All options are given explicitly as an specific config file.
 if opt.kilosort == 2,     master_kilosort(input, opt)
 elseif opt.kilosort == 4, master_kilosort4(input, opt)
 end, close all
%[text] %[text:anchor:H_5229] ### <u>Bombcell Processing</u>
%[text] Any Bombcell options are given explicitly as an specific config file.
 if opt.bombcell, Bombcell_Main(input, opt), end
%[text] %[text:anchor:H_100a] ### <u>Phy Programatic Call</u>
%[text] Parameters for proper Phy function are created at the time of kilosort processing. Check them!
 if opt.phy, cd(opt.FolderProcDataMat)
   system('phy template-gui params.py'); % This will keep Matlab on HOLD until the Phy instance is closed!
 end
%[text] %[text:anchor:H_3fef] ### <u>Blob tracking</u> 
%[text] (Social Learning-specific Method). Optimiced for videos from central cenital camera. 
 if opt.offlineTrack, blob = processAndTrack_video(opt);
  save(fullfile(opt.analysis, "blob.mat"), 'blob');
 end
%[text] %[text:anchor:T_6418] # <u>02.</u> Beging Post-phy processing.
%[text] %[text:anchor:H_2366] ### <u>Spike Data</u> Processing
 if opt.doSpikething
%[text] %[text:anchor:H_3be1] ####      Recover Spike clusters after sorting and curation. (General Method)
  spike = loadSpikes(opt);
  if isfield(spike, "spike"), spike = spike.spike; end
  save(fullfile(opt.spikeSorted, "spike.mat"), 'spike', '-mat');
%[text] %[text:anchor:H_3662] ####      Parse cluster's spikes into trials, using time definitions. (General Method)
  if ~exist('trialdef', 'var'), load(fullfile(opt.trialSorted, "trialdef.mat")); end
  [neurons, ~] = sort2trials(spike, trialdef, opt); % TODO fix Fieldtrip extraction
%[text] %[text:anchor:H_068a] ####      Spiking indexing for Social intereactions.
  if opt.useTrack
    if ~isfield(neurons, "interactions")
       [neurons.interactions, ~] = sort2trials(spike, blob.Merges, opt); end
%[text] %[text:anchor:H_0325] ####      Use video human-assesment to extract Social Event indexing.
    if ~isfield(events,"social"), [events.social] = getSocialEvents(opt);
       save(fullfile(opt.analysis, "events.mat"), 'events', '-mat'); end
  end
%[text]      Save output to \\analysis
  save(fullfile(opt.analysis, "neurons.mat"), 'neurons', '-mat')
 end
%[text] %[text:anchor:H_3f92] ### <u>LFP Data</u> Processing Pipeline.
 if opt.doLFPthing % IN PROGRESS) 
     LFP_Fieldtrip(neurons, spike, trialdef, input, opt)
 end
%[text] %[text:anchor:H_6e82] ### <u>Aggregate</u> all Spike data as cell arrays with session/subject indexing.
    allneurons{x,y}     = neurons;      % Collect Subjects/session data
    allevents{x,y}      = events;       % idem
    allconditions{x,y}  = conditions;   % idem
    allspike{x,y}       = spike;        % idem
    if opt.plot_SocLear, allblobs{x,y} = blob; end % Include blob tracking if Social
    clear neurons events conditions spike blob % Clean up Subjects/Session results after collection
 end % Sessions loop
 end % Subjects loop
%[text]  Save all preprocessed data.
  save(fullfile(input.analysis, "data_all.mat"),'allspike','allconditions','allevents','allneurons');
  if exist('allblobs','var'), save(fullfile(input.analysis, "data_all.mat"), 'allblobs','-append');end
end % End preprocessing line.
%%
%[text] %[text:anchor:T_93c6] # <u>03</u>. Plotting. 
%[text] %[text:anchor:T_610e] (TODO: Check what is not pure plotting and take it out of plotting scripts)
%[text] %[text:anchor:H_14b0] #### Create a plot directory
if ~exist(fullfile(opt.analysis,'plots'),"dir"), mkdir(fullfile(opt.analysis,'plots')), end
%[text] %[text:anchor:H_994e] #### Recover collected data
%[text] %[text:anchor:H_820b] if we skipped the preprocessing.
if opt.jump2plot % If we skipped preprocessing, we need to recoved data
    [allneurons, allevents, allconditions, allspike, allblobs] = collect_data(input, opt)
end
%%
%[text] %[text:anchor:H_47f8] ## Add your plotting scripts here. 
%[text] Make sure they are in accordance with the to-plot list set above.
% NGL_plotting01
% NGL_plotting02

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":33.1}
%---
%[control:filebrowser:6d69]
%   data: {"browserType":"Folder","defaultValue":"\"C:\\Code\\ephys-data-pipeline\"","label":"Toolbox path","run":"Nothing"}
%---
%[control:filebrowser:8d93]
%   data: {"browserType":"Folder","defaultValue":"\"D:\\\"","label":"Project Location","run":"Nothing"}
%---
%[control:editfield:77cf]
%   data: {"defaultValue":"'all';","label":"subjects","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:editfield:8f23]
%   data: {"defaultValue":"'all';","label":"dates","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:slider:0798]
%   data: {"defaultValue":32,"label":"numChannels","max":512,"min":32,"run":"Nothing","runOn":"ValueChanged","step":32}
%---
%[control:checkbox:5e20]
%   data: {"defaultValue":true,"label":"CAR","run":"Nothing"}
%---
%[control:dropdown:138b]
%   data: {"defaultValue":"0","itemLabels":["OFF","50","60"],"items":["0","50","60"],"label":"linefilter","run":"Nothing"}
%---
%[control:checkbox:1ac4]
%   data: {"defaultValue":true,"label":"bin","run":"Nothing"}
%---
%[control:slider:0bf3]
%   data: {"defaultValue":400,"label":"Highpass","max":500,"min":150,"run":"Nothing","runOn":"ValueChanged","step":50}
%---
%[control:checkbox:835c]
%   data: {"defaultValue":false,"label":"Fieldtrip","run":"Nothing"}
%---
%[control:slider:65af]
%   data: {"defaultValue":250,"label":"Highpass","max":300,"min":50,"run":"Nothing","runOn":"ValueChanged","step":25}
%---
%[control:checkbox:9b7b]
%   data: {"defaultValue":false,"label":"GetMotionSensors","run":"Nothing"}
%---
%[control:checkbox:5fc5]
%   data: {"defaultValue":true,"label":"RetrieveEvents","run":"Nothing"}
%---
%[control:editfield:188c]
%   data: {"defaultValue":"{'itiOn', 'stimOn1', 'rwd'}","label":"alignto","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:editfield:954b]
%   data: {"defaultValue":"{'itiOn', 'stimOn1', 'rwd'}","label":"alignto","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:slider:3fad]
%   data: {"defaultValue":1500,"label":"addtime","max":5000,"min":0,"run":"Nothing","runOn":"ValueChanged","step":500}
%---
%[control:dropdown:760b]
%   data: {"defaultValue":"4","itemLabels":["2","4"],"items":["2","4"],"label":"kilosort","run":"Nothing"}
%---
%[control:editfield:061a]
%   data: {"defaultValue":"\"chanMapXXX.mat\"","label":"KSchanMapFile","run":"Nothing","valueType":"String"}
%---
%[control:slider:13e2]
%   data: {"defaultValue":-6,"label":"spkTh","max":-4,"min":-10,"run":"Nothing","runOn":"ValueChanged","step":1}
%---
%[control:dropdown:4406]
%   data: {"defaultValue":"0","itemLabels":["0","2","4"],"items":["0","2","4"],"label":"kilosort","run":"Nothing"}
%---
%[control:checkbox:68a4]
%   data: {"defaultValue":false,"label":"phy","run":"Nothing"}
%---
%[control:checkbox:4be0]
%   data: {"defaultValue":false,"label":"","run":"SectionToEnd"}
%---
%[control:checkbox:6bbb]
%   data: {"defaultValue":false,"label":"Raster Trial by trial","run":"Nothing"}
%---
%[control:checkbox:4cf3]
%   data: {"defaultValue":false,"label":"Raster Aligned","run":"Nothing"}
%---
%[control:checkbox:6711]
%   data: {"defaultValue":true,"label":"Normalized Fire Rate","run":"Nothing"}
%---
%[control:checkbox:24ce]
%   data: {"defaultValue":false,"label":"YYY","run":"Nothing"}
%---
%[control:checkbox:3610]
%   data: {"defaultValue":false,"label":"ZZZ","run":"Nothing"}
%---
%[control:checkbox:88e9]
%   data: {"defaultValue":false,"label":"Ready!","run":"Nothing"}
%---
%[control:button:97e9]
%   data: {"label":"Run Section","run":"SectionAndStaleSectionsAbove"}
%---
%[control:button:8b0b]
%   data: {"label":"Continue","run":"SectionAndStaleSectionsAbove"}
%---
%[control:button:3bbe]
%   data: {"label":"Continue","run":"SectionAndStaleSectionsAbove"}
%---
%[output:50f8331b]
%   data: {"dataType":"textualVariable","outputData":{"header":"struct with fields:","name":"opt","value":"         numChannels: 32\n                 CAR: 1\n          linefilter: 0\n                 bin: 1\n            highpass: 400\n           FieldTrip: 0\n             lowpass: 250\n    GetMotionSensors: 0\n      RetrieveEvents: 1\n             alignto: {'itiOn'  'stimOn1'  'rwd'}\n            trEvents: {'na1'  'na2'}\n             addtime: 1500\n            kilosort: 4\n       KSchanMapFile: \"chanMapXXX.mat\"\n               spkTh: -6\n            bombcell: 0\n               rerun: 1\n                 phy: 0\n               StpSz: 1000000\n          parsetrial: 0\n              useNWB: 0\n        doSpikething: 1\n          doLFPthing: 1\n        offlineTrack: 0\n                FLIP: 0\n               getwF: 1\n           gwfparams: [1×1 struct]\n          plot_trial: 0\n          plot_align: 0\n       plot_sessions: 0\n       plot_fireRate: 1\n           jump2plot: 0\n                plot: [1×1 struct]\n"}}
%---
%[output:5fd23422]
%   data: {"dataType":"text","outputData":{"text":"Folder system for project \"S3_ExperimentArena\" already exists. Nothing changed. \n","truncated":false}}
%---
%[output:296e2cb5]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Did you put all your customized files in its place? Did you check the box?"}}
%---
