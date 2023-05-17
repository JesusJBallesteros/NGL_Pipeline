function [spike, template] = read_KSresults(varargin)
% Adaptation of 'spikes' toolbox from 'Cortex-Lab' at ULC
% (https://github.com/cortex-lab/spikes)
% These functions make it easy to import spioke data after kilosort/phy
% processing.
%
% INPUT: 
%
%
% OUTPUT:   spike: structure containing:
%               .st, are sp times in seconds.
%               .clu, are cluster identities.
%               .cgs, are cluster groups, the labels given during manual sorting in phy (1=MUA, 2=Good, 3=Unsorted).
%               .cids, cluster-sp.cgs entry index. e.g sp.cgs(sp.cids==9) gives the cluster group for cluster 9
%               .
%               .
%               .
%
%           template: structure containing:
%               .
%               .
%               .
%               .
%               .
%               .
%               .

if nargin < 1, opt = struct();
elseif nargin == 1, opt = varargin{1};
end

%% Defaults
if ~isfield(opt,'UseEvents') || isempty(opt.UseEvents),                 opt.UseEvents      = false;                  end
if ~isfield(opt,'plotdrift') || isempty(opt.plotdrift),                 opt.plotdrift      = false;                   end
if ~isfield(opt,'plotAmpDepth') || isempty(opt.plotAmpDepth),           opt.plotAmpDepth   = false;                   end

%% Loading data from kilosort/phy easily
% spikes from clusters labeled "noise" will be omitted
spike = loadKSdir(opt.PathRaw);
template = [];

%% Computing some useful details about spikes/neurons
%
[spike.Amps, ~, template.Ypos, template.Amps, template.UnW, template.Dur, template.PeakWF] = ...
    templatePositionsAmplitudes(spike.temps, spike.winv, spike.ycoords, spike.spikeTemplates, spike.tempScalingAmps);

%
[spike.Times, spike.Amps, spike.Depths, spike.Sites] = ksDriftmap(opt.PathRaw);

%% Plot Drift. 
if opt.plotdrift
    % To observe whether there was drift 'ksDriftmap' is useful. y-axis is depth.
    figure; 
    plotDriftmap(spike.Times, spike.Amps, spike.Depths);
end

%% Amplitudes and Depth
if opt.plotAmpDepth
    % Where spikes of different amplitudes were recorded. Colormap of the 
    % spike distribution across depth and amplitude.
    ampBins = 0:30:min(max(spike.Amps),800);
    
    % This depends on the probe
    depthBins = 0:15:190; % check for this


    [spike.AmpPDFS, spike.FrCDFS] = computeWFampsOverDepth(spike.Amps, spike.Depths, ampBins, depthBins, spike.st(end));
    plotWFampCDFs(spike.AmpPDFS, spike.FrCDFS, ampBins, depthBins);
end

end