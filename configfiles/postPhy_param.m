%% A config file to collect possible options for the Post-phy pipeline.
% In order to not crowd the main user script, and given the fact that 
% at this point the interests for individual projects could start to 
% divert from others. Add this file to the other config/parameter files
% under your '\analysisCode' folder

%% BRAIN AREA - MAP key
% Set an area pointer, according to your k-coords shank definition
opt.mapkey = struct();
    opt.mapkey.NCL = [1,2];
    opt.mapkey.STR = 3;

%% Spike analysis and plots 
% Proceed to some analysis and plots for clustered units obtained from
% KS-Phy processing.
% % Working.
opt.doSpikething = true;
	opt.getwF    = true;
	opt.useTrack = false;
    
% On development
opt.neurDyn = struct();
    opt.neurDyn.do      = true;
    opt.neurDyn.binsize = 0.1;
    opt.neurDyn.method  = "tSNE";
    opt.neurDyn.synth   = 1;
	
%% General options to extract single waveforms 
% from the clustered units. 
% Suboptions are probably to held fix for everyone.
% % Working
opt.gwfparams.nWf       	= 2000;     % maximum # of waveforms to extract per cluster. Affects computing time.
opt.gwfparams.dataType      = 'int16';  % Data type of .dat file
opt.gwfparams.nCh           = 32;       % Number of channels that were streamed in .dat file
opt.gwfparams.wfWin         = [-20 41]; % Number of samples around spiketime to include in waveform
opt.gwfparams.nWf           = 1;        % Proportion of total waveforms per unit to extract

%% LFP analysis and plots 
% from data obtained via Fieldtrip pathway.
% % Not totally functional yet. 
opt.doLFPthing = false;

%% Run blob ID and tracking on behavioral videos.
% So far, optimized for Social learning paradigm in half plus arena.
% Tracks simultaneous pigeons in the arena , measuring the distance between
% blobs and categorizing as 'interactions' those when the blob merge into a
% single one. 
% % Working.
opt.offlineTrack = false;

%% FLIP tests.
% First approach to the 'Beta-gamma cortical motif' paper from Miller lab,
% without much results yet but also no deeply investigated.
% % On development
opt.FLIP = false;

%% Plotting
% These parameters affect the plotting functions used after unit sorting.
param   = struct('visible',      'off', ... % figure visibility at plotting
	            'Resolution',    300, ...
		        'treatment',     false, ... % plot different levels due to treatment/s
		        'baseline',      1500, ...  % for event-aligned plots, msec before event
		        'post',          2000, ... % for event-aligned plots, msec after event
		        'plotevent',     3, ...    % mark these trial events in rasters
		        ... % PSH plotting
		        'binSize',       200, ...  % binning window, msec. Only for aligments other than itiOn
		        'stepSz',        20);      % window steps, msec. Only for aligments other than itiOn
		        % ... %
		        % 'raster', struct(), ...
		        % ... %
		        % 'psh', struct(), ...
		        % ... %
		        % 'poolraster', struct() ...
		        % );
		        % %  ... %

% % % Examples for other projects
% % Extintion specific. 
% param.trial2plot = 'allInitiated'; % 'correct', 'incorrect', 'omission', 'allInitiated'
% param.IncludeFS = false;
% param.FS2plot   = true; 
% param.binSize   = 500;
% param.stepSz    = 50;
% param.nBlocks   = 8;
% param.interval  = [-2000 10000];
% param.baseline  = 0 - param.interval(1);

% %
%
%
%
