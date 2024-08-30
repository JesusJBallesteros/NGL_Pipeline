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

%% LFP analysis and plots 
% from data obtained via Fieldtrip pathway.
% % Not totally functional yet. 
opt.doLFPthing = false;

%% FLIP tests.
% First approach to the 'Beta-gamma cortical motif' paper from Miller lab,
% without much results yet but also no deeply investigated.
% % On development
opt.FLIP = false;

%% General options to extract single waveforms 
% from the clustered units. 
% Suboptions are probably to held fix for everyone.
% % Working
opt.getwF                   = true;
    opt.gwfparams.dataType      = 'int16';  % Data type of .dat file
    opt.gwfparams.nCh           = 32;       % Number of channels that were streamed in .dat file
    opt.gwfparams.wfWin         = [-20 41]; % Number of samples around spiketime to include in waveform
    opt.gwfparams.nWf           = 1;        % Proportion of total waveforms per unit to extract

%% Plotting
% These parameters affect the plotting functions used after unit sorting.
% % Working
param = struct('res',           true, ...  % To deprecate % load 'res' variable from Juan's
               'visible',       'off', ... % figure visibility at plotting
               'treatment',     true, ... % plot different levels due to treatment/block/phase
               'genstats',      true, ... % do raster plots
               'pooledstats',   true, ... % Plot a pool of all clusters together
               ... % raster plotting
               'plotcol',       [.4 .4 .4; 0.6350 0.0780 0.1840; 0 0 0], ... % levels coloring
               'plotStyle',     'square', ...  % plot marker
               'spkWidth',      3, ... % marker size
               'lineLength',    1, ...  % trial line width
               'timelim',       [-500 2000], ... % limits on msec around cero, to plot
               ... % PSH plotting
               'binSize',       100, ... % binning window, msec
               'stepSz',        10, ... % window running step, msec
               'smpRate',       1000, ... % samples per second in 'neurons'
               'interval',      [0 2500]); % time interval to calculate histogram bins
% These sub-structures could contain the options that currently are inside the plot functions, for further customization
%                     ... %
%                     'raster', struct(), ... 
%                     ... %
%                     'psh', struct(), ...
%                     ... %
%                     'poolraster', struct() ...
%                     );
%                   %  ... %