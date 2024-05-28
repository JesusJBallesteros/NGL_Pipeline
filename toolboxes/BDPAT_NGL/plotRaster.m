function trialCounter = plotRaster(spikes,trialCounter,plotCol,varargin)
%%function trialCounter = plotRaster(spikes,trialCounter,plotCol,varargin)
%
% This function creates a raster plot for the given spike train.
%
%INPUTS
%  * 'spikes'          : cell containing aligned spike times per trial
%  * 'trialCounter'    : arbitrary trial number to start plotting at
%  * 'plotCol'         : RGB color triplet for plot color
%
%OPTIONAL INPUTS
%  * 'spkWidth'        : size of the marker for plotting (default is 5)
%  * 'plotStyle'       : string of marker style for plotting (default is
%                        'lines')
%
%OUTPUTS
%   * 'trialCounter'   : arbitrary trial number until which spikes were
%                        plotted

% VERSION HISTORY:
% Author:         Lukas Hahn
% Version:        1.1.1
% Last Change:    11.12.2023
%
% 15.07.2019, Lukas: v1.0.0 release version
% 28.11.2023, Lukas: v1.1.0 added plot as square markers option
% 11.12.2023, Lukas: v1.1.1 updated documentation
%%
plotColBkp = plotCol;
if nargin==3
    spkWidth = 5;
    plotStyle = 'lines';
else
    spkWidth = varargin{1};
    plotStyle = varargin{2};
end
for trial=1:size(spikes,1) %for all trials
    if numel(spikes{trial,1})<4
        plotCorrection = [NaN; NaN; NaN; NaN]; %to correct for line bugs
    else
        plotCorrection = [];
        plotCol = plotColBkp;
    end
    trialCounter = trialCounter+2;
    spikeTrains = [plotCorrection; spikes{trial,1}];

    %helper = randperm(size(spikeTrains,1));
    %reduces raster to 20 %
    %spikeTrains = spikeTrains(helper(1:ceil(size(helper,1)/5)+2));
    if strcmp(plotStyle,'lines')
        line([spikeTrains spikeTrains],[trialCounter trialCounter+2],...
            'Color',plotCol,'LineWidth',spkWidth)
    else
        plot(spikeTrains,trialCounter,'s','MarkerEdgeColor',plotCol,...
            'MarkerFaceColor',plotCol,'MarkerSize',spkWidth)
        hold on
    end
end
end