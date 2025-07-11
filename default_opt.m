%%  Create the structure 'opt' with default values. IN PROGRESS
% DO NOT MODIFY here, use input values in script,
% whose will run after this is executed.
opt = struct();
%% General
opt.numChannels      = 32;
opt.CAR              = 1;
opt.linefilter       = 0;
opt.bin              = true;
opt.highpass         = 400;
opt.FieldTrip        = false;
opt.lowpass          = 250;
opt.GetMotionSensors = false;
opt.RetrieveEvents   = true;
opt.alignto          = {'itiOn', 'stimOn1', 'rwd'};
opt.trEvents         = {'na1'};
opt.addtime          = 1500;
opt.kilosort         = 4;
opt.KSchanMapFile    = 'chanMapE32-S2_DeutSN11.mat';
opt.spkTh            = -6;
opt.bombcell         = false;
opt.rerun            = true;
opt.phy              = false;
%% 01. Deuteron pipeline
opt.StpSz            = 1000000;
opt.parsetrial       = false;
%% 01. INTAN pipeline 
opt.useNWB           = false;
%% 02. postPhy
opt.doSpikething     = true;
opt.doLFPthing       = true;
opt.offlineTrack     = false;
opt.FLIP             = false;
% opt.multregion     = true;
opt.getwF            = true;	
opt.gwfparams.nWf       = 2000;     % maximum # of waveforms to extract per cluster. Affects computing time.
opt.gwfparams.dataType  = 'int16';  % Data type of .dat file
opt.gwfparams.nCh       = 32;       % Number of channels that were streamed in .dat file
opt.gwfparams.wfWin     = [-20 41]; % Number of samples around spiketime to include in waveform
%% Plotting
param = struct('visible',      'off', ... % figure visibility at plotting
		       'Resolution',    300, ...  % Quality
               'size',          'adaptive', ... % 'adaptive' for screen's full resolution, or [W H]
		       'treatment',     false, ... % plot different levels due to treatment/s
		       'baseline',      1500, ... % data before event-aligned 0 (msec)
		       'post',          2500, ... % data after event-aligned 0 (msec)
		       'plotevent',     [2,3,7], ... % draw these events in raster trial plots  
               ... % Rasters
               'plotStyle',     'lines', ... % lines/dots
               'spkWidth',      .5, ... % line/dot width
               'lineLength',    1, ... % line length
		       ... % PSH
		       'smpRate',       1000, ... % To perform calculations on PSH plots (Hz)
               'binSize',       500, ...  % binning window (msec)
		       'stepSz',        50 ...    % window steps (msec)
               );     
% These parameters affect the plotting functions used after unit sorting.
opt.plot_trial          = false;   % in general, to plot full trial rasters/PSHS
opt.plot_align          = false;   % in general, to plot event-aligned pieces of trial
opt.plot_sessions       = false;  % in general, to plot data collected across session
opt.plot_fireRate       = true;
%% README.TXT 
% It contains details about the project. File can also be modified later.
readmecontent = ["Study name: Preparation for social learning paradigms", ...
                 "Readme date: 22/02/2024"                          , ...
                 "Person (1) responsible for data repository: Jesus", ...
                 "Person(s) responsible for study: Juan, Jesus "    , ...
                 "Hardware used: Deuteron. E32-S2"                  , ...
                 "Related Publication(s): None."                    , ...
                 "Short description of study: Autoshaping with social learning." ];
% end