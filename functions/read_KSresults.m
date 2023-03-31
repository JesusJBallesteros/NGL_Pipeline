function [sp] = read_KSresults(varargin)
% Adaptation of 'spikes' toolbox from 'Cortex-Lab' at ULC
% (https://github.com/cortex-lab/spikes)
% These functions make it easy to import spioke data after kilosort/phy
% processing.
%
% INPUT: 
%
%
% OUTPUT:   sp: structure containing:
%               sp.st, are sp times in seconds.
%               sp.clu, are cluster identities.
%               sp.cgs, are cluster groups, the labels given during manual sorting in phy (1=MUA, 2=Good, 3=Unsorted).
%               sp.cids, cluster-sp.cgs entry index. e.g sp.cgs(sp.cids==9) gives the cluster group for cluster 9

if nargin < 1, opt = struct();
elseif nargin == 1, opt = varargin{1};
end

%% Defaults
if ~isfield(opt,'FolderProcDataMat'), opt.FolderProcDataMat  = pwd + "\processed"; end
if ~isfield(opt,'UseEvents'),         opt.UseEvents      = false;                  end
if ~isfield(opt,'drift'),             opt.drift          = true;                   end
if ~isfield(opt,'amplitude'),         opt.amplitude      = true;                   end
if ~isfield(opt,'psth'),              opt.psth           = false;                  end

%% Loading data from kilosort/phy easily
% spikes from clusters labeled "noise" will be omitted
sp = loadKSdir(opt.FolderProcDataMat);

%% Load synchronization data
if opt.UseEvents
    % Eventually, will have events to align to.
    % EventTimes = load('C:\...\data\EventTimes.mat'); % a vector of times in seconds of some event to align to

    % eventTimes = spikeGLXdigitalParse(syncDat, lfpFs);
    
    % - eventTimes{1} contains the sync events from digital channel 1, as three cells: 
    % - eventTimes{1}{1} is the times of all events
    % - eventTimes{1}{2} is the times the digital bit went from off to on
    % - eventTimes{1}{2} is the times the digital bit went from on to off
    
    % To make a timebase conversion, e.g. between two probes:
    % [~,b] = makeCorrection(syncTimesProbe1, syncTimesProbe2, false);
    
    % and to apply it:
    % correctedSpikeTimes = applyCorrection(spikeTimesProbe2, b);
end

%% Some (perhaps) convenient checks:
% Drift. 
if opt.drift
    % To observe whether there was drift 'ksDriftmap' is useful. y-axis is depth.
    [sp.Times, sp.Amps, sp.Depths, sp.Sites] = ksDriftmap(opt.FolderProcDataMat);
    
    figure; 
    plotDriftmap(sp.Times, sp.Amps, sp.Depths);
end

% sp amplitudes.
if opt.amplitude
    % Where spikes of different amplitudes were recorded. Colormap of the sp distribution across depth and amplitude.
    depthBins = 0:15:250; % check for this
    ampBins = 0:30:min(max(sp.Amps),800);
    sp.st = sp.st(end);
    
    [sp.pdfs, sp.cdfs] = computeWFampsOverDepth(sp.Amps, sp.Depths, ampBins, depthBins, sp.st);
    plotWFampCDFs(sp.pdfs, sp.cdfs, ampBins, depthBins);
end

% PSTHs
if opt.psth && opt.UseEvents
    % With the relevant 'eventTimes' (not the cell array as above, just a vector):
    window = [-0.3 1]; % look at sp times from 0.3 sec before each event to 1 sec after
    
    % If your events come in different types, like different orientations of a
    % visual stimulus, then you can provide those values as "trial groups",
    % which will be used to construct a tuning curve. Here we just give a
    % vector of all ones. 
    trialGroups = ones(size(eventTimes)); 
    
    % use left/right arrows to page through the clusters
    psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups);

% PSTHs as function of depth
    depthBinSize = 80;      % in units of the channel coordinates, in this case µm
    timeBinSize = 0.01;     % seconds
    bslWin = [-0.2 -0.05];  % window in which to compute "baseline" rates for normalization
    psthType = 'norm';      % show the normalized version
    eventName = 'stimulus onset'; % for figure labeling
    
    [timeBins, depthBins, allP, ~] = psthByDepth(spikeTimes, spikeDepths, ...
                                                 depthBinSize, timeBinSize, ...
                                                 eventTimes, window, bslWin);
    
    figure;
    plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);
    
end

end