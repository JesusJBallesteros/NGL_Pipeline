%% A config file to collect possible options for the Post-phy pipeline.
% In order to not crowd the main user script, and given the fact that 
% at this point the interests for individual projects could start to 
% divert from others. Add this file to the other config/parameter files
% under your '\analysisCode' folder

%% Run blob ID and tracking on behavioral videos.
% So far, optimized for Social learning paradigm in half plus arena.
% Tracks simultaneous pigeons in the arena , measuring the distance between
% blobs and categorizing as 'interactions' those when the blob merge into a
% single one. 
% % Working.
opt.offlineTrack = false;

%% Spike analysis and plots 
% Proceed to some analysis and plots for clustered units obtained from
% KS-Phy processing.
% % Working.
opt.doSpikething = true;
opt.getwF        = true;	
	
%% General options to extract single waveforms 
% from the clustered units. 
% Suboptions are probably to held fix for everyone.
% % Working
opt.gwfparams.nWf       	= 2000;     % maximum # of waveforms to extract per cluster. Affects computing time.
opt.gwfparams.dataType      = 'int16';  % Data type of .dat file
opt.gwfparams.nCh           = 32;       % Number of channels that were streamed in .dat file
opt.gwfparams.wfWin         = [-20 41]; % Number of samples around spiketime to include in waveform

%% LFP analysis and plots 
% from data obtained via Fieldtrip pathway.
% % Not totally functional yet. 
opt.doLFPthing = false;

%% FLIP tests.
% First approach to the 'Beta-gamma cortical motif' paper from Miller lab,
% without much results yet but also no deeply investigated.
% % On development
opt.FLIP = false;

%% Plotting
% These parameters affect the plotting functions used after unit sorting.
opt.plot_trial              = false;   % in general, to plot full trial rasters/PSHS
opt.plot_align              = false;   % in general, to plot event-aligned pieces of trial
opt.plot_sessions           = false;  % in general, to plot data collected across session
opt.plot_fireRate           = true;

%% General figure-related parameters
param = struct('visible',      'off', ... % figure visibility at plotting
		       'Resolution',    300, ...  %
               'size',          'adaptive', ... %
		       'treatment',     true, ... % plot different levels due to treatment/s
		       'baseline',      1500, ... % for event-aligned plots, msec before event
		       'post',          2500, ... % for event-aligned plots, msec after event
		       'plotevent',     [2,3,7], ... % mark these trial events in rasters
               'smpRate',       1000, ... %
		       ... % PSH plotting
		       'binSize',       500, ...  % binning window, msec. Only for full trial
		       'stepSz',        50 ...    % window steps, msec. Only for full trial
               );     

% Extintion specific
param.trial2plot = 'allInitiated'; % 'correct', 'incorrect', 'omission', 'allInitiated'
param.IncludeFS = false;
param.FS2plot   = true; 
param.binSize   = 500;
param.stepSz    = 50;
param.nBlocks   = 8;
param.interval  = [-2000 10000];
param.baseline  = 0 - param.interval(1);

