function [spike, template] = plot_KSresults(spike, varargin)
% Adaptation of 'spikes' toolbox from 'Cortex-Lab' at ULC. https://github.com/cortex-lab/spikes
% These functions make it easy to import spike data after kilosort/phy processing.
%
% INPUT:    opt: with the different options including folders and tasks to
%             perform.
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

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Defaults
if ~isfield(opt,'UseEvents') || isempty(opt.UseEvents),                 opt.UseEvents      = false;                 end
if ~isfield(opt,'plotdrift') || isempty(opt.plotdrift),                 opt.plotdrift      = false;                 end
if ~isfield(opt,'plotAmpDepth') || isempty(opt.plotAmpDepth),           opt.plotAmpDepth   = false;                 end
if ~isfield(opt,'excludeNoise') || isempty(opt.excludeNoise),           opt.excludeNoise   = true;                  end
if ~isfield(opt,'loadPCs') || isempty(opt.loadPCs),                     opt.loadPCs        = false;                 end
% if ~isfield(opt,'optionalfield') || isempty(opt.optionalfield),           opt.optionalfield = false;                 end
% if ~isfield(opt,'optionalfield') || isempty(opt.optionalfield),           opt.optionalfield = false;                 end

% Parameters, optional. We set them up outside as 'opt'
params = struct();
    params.excludeNoise = opt.excludeNoise;  % if true, Spikes from clusters labeled "noise" will be omitted
    params.loadPCs      = opt.loadPCs; % if true, PC features will be ommited.
    params.useGood      = true;
    params.FRthresh     = 0.1;
    params.ampThresh    = 20;

%% Computing some basic details for spikes and templates
[spike.Amps, ~, template.Ypos, template.Amps, template.UnW, template.Dur, template.PeakWF] = ...
    templatePositionsAmplitudes(spike.temps, spike.winv, spike.ycoords, spike.spikeTemplates, spike.tempScalingAmps);

%% Plot Drift. 
if opt.plotdrift
    [spike.Times, ~, spike.Depths, spike.Sites] = ksDriftmap(opt.PathRaw);

    % To observe whether there was drift 'ksDriftmap' is useful. y-axis is depth.
    toplot = 'mark'; % Can be '' (empty), 'mark' or 'show'
    
    f1 = figure;
    plotDriftmap(spike.Times, spike.Amps, spike.Depths, toplot);
    title('Drift map')
    
        saveas(f1, fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_DepthAmpl']), 'png');
        saveas(f1, fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_DepthAmpl']), 'epsc');
end

%% Amplitudes and Depth
if opt.plotAmpDepth
    % Where spikes of different amplitudes were recorded. Colormap of the 
    % spike distribution across depth and amplitude.
    
    % This depends on the probe
    ampBins = 0:20:150; %min(max(spike.Amps),800);   
    depthBins = 0:25:250; % check for this

    [spike.AmpPDFS, spike.FrCDFS] = computeWFampsOverDepth(spike.Amps, spike.Depths, ampBins, depthBins, spike.st(end));
    
    f2 = plotWFampCDFs(spike.AmpPDFS, spike.FrCDFS, ampBins, depthBins);
        saveas(f2, fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_Amplpdf']), 'png');
        saveas(f2, fullfile(opt.FolderProcDataMat, [opt.SavFileName, '_Amplpdf']), 'epsc');

end

%% Plot Unit Scatter

    % f3 = plotUnitScatter(spike, template, params);

end