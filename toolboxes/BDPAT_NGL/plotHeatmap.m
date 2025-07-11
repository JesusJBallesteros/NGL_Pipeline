function trialCounter = plotHeatmap(spikes,trialCounter,varargin)
%% Under development
%
% This function creates a raster plot for the given spike train.
%
%INPUTS
%  * 'spikes'          : cell containing aligned spike times per trial
%  * 'trialCounter'    : arbitrary trial number to start plotting at
%
%OPTIONAL INPUTS
% %  * 'plotCol'         : RGB color triplet for plot color, single row
% %                        vector (results in all spike indicators having the
% %                        same color, default), or a matrix of RGB triplets
% %                        (with each triplet corresponding to a trial,
% %                        resulting in trial unique colors)
% %  * 'spkWidth'        : size of the marker for plotting (default is 3)
%  * 'plotStyle'       : string of marker style for plotting (default is
%                        'lines')
% %  * 'lineLength'      : vertical length of line indicator (default is 1)
%  * 'timelim'         : int array, plot spikes only between [-timelim(1),
%                        timelim(2)]. Default is empty.
%
%OUTPUTS
%   * 'trialCounter'   : arbitrary trial number until which spikes were
%                        plotted

% VERSION HISTORY:
% Author:         Jesus Ballesteros
% Version:        0.1
% Last Change:    08.07.2024
%
% 08.07.2024, Jesus: 

%%
%default values
% plotCol = zeros(1,3);
% spkWidth = 3;
plotStyle = 'lines';
% lineLength = 1;
timelim = [];

%user specified values
if nargin>2
    if nargin==5 %ensures backward compatibility to use as
        % plotRaster(spikes,trialCounter,plotColor,spikeWidth,'plotStyle')
%         plotCol = varargin{1};
%         spkWidth = varargin{2};
        plotStyle = varargin{3};
    else
        for i=1:length(varargin)
            if isa(varargin{i},'char') || isa(varargin{i},'string')
                switch lower(varargin{i})
                    case 'plotcol'
                        plotCol = varargin{i+1};
                    case 'spkwidth'
                        spkWidth = varargin{i+1};
                    case 'plotstyle'
                        plotStyle = varargin{i+1};
                    case 'linelength'
                        lineLength = varargin{i+1};
                    case 'timelim' % Jesus, added on 04.06.2024
                        timelim = varargin{i+1};
                    otherwise
                        %next input
                end
            else
                %next input
            end
        end
    end
end
%%
%ensure that every trial has a corresponding color
if size(plotCol,1)<size(spikes,1)
    plotCol = repmat(plotCol(1,:),size(spikes,1),1);
end
%%
for trial = 1:size(spikes,1) % for all trials
    % Jesus, added timelim optional input, used here
    if ~isempty(timelim)
        spikestoplot = spikes{trial,1}(spikes{trial,1}>(timelim(1)) & spikes{trial,1}<(timelim(2)));
    else
        spikestoplot = spikes{trial,1};
    end

    if numel(spikestoplot)<4
        plotCorrection = [NaN; NaN; NaN; NaN]; %to correct for line bugs
    else
        plotCorrection = [];
    end
    spikeTrains = [plotCorrection; spikestoplot];

    %helper = randperm(size(spikeTrains,1));
    %reduces raster to 20 %
    %spikeTrains = spikeTrains(helper(1:ceil(size(helper,1)/5)+2));
    if strcmp(plotStyle,'lines')
        line([spikeTrains spikeTrains],...
            [trialCounter trialCounter+lineLength],...
            'Color', plotCol(trial,:), 'LineWidth', spkWidth)
    else
        plot(spikeTrains,trialCounter, 's', 'MarkerEdgeColor',...
            plotCol(trial,:), 'MarkerFaceColor', plotCol(trial,:),...
            'MarkerSize', spkWidth)
        hold on
    end

    trialCounter = trialCounter+1;
end

end