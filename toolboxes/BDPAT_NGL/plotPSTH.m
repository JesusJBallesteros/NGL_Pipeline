function upperY = plotPSTH(spikes,stepSz,binSize,interval,smpRate,varargin)
%%
%upperY = plotPSTH(spikes,stepSz,binSize,interval,smpRate,varargin)
%
% Use this function to calculate the firing rate of a neuron over the time
% course of the trial.
%
%INPUTS
%  * 'spikes'          : cell containing aligned spike times per trial
%  * 'stepSz'          : step size of for the walk through the interval of
%                        interest (in samples)
%  * 'binSize'         : size of the bins in which the spike train will be
%                        divided (in samples)
%  * 'interval'        : matrix containing the start and the end point(s)
%                        (in samples) of the interval of interest, relative
%                        to the alignment
%  * 'smpRate'         : sampling rate of the recording (in Hz)
%
% OPTIONAL INPUTS
% varargin parameter input pairs
%  * 'plotCol'         : RGB color triplet for plot color
%  * 'meanline'        : string for mean linestyle (e.g. '-', '--', '.-')
%  * 'smoothPlot'      : 1 - yes, or 0 - no, determines if plot is smoothed
%                        by a standard window of 5
%OUTPUTS
%  * 'upperY'          : largest value of
%                        nanmean(fireRate)+nanstd(fireRate)

% VERSION HISTORY:
% Author:         Lukas Hahn
% Version:        1.1.3
% Last Change:    25.01.2024
%
% 15.07.2019, Lukas: v1.0.0 release version
% 12.08.2019, Lukas: v1.1.0 added handling for no spikes (returns a
%                           2-by-interval length matrix of nan for
%                           plotting)
% 27.06.2020, Lukas: v1.1.1 added output to determine upper bound of ylim
% 11.12.2023, Lukas: v1.1.2 updated documentation
% 25.01.2024, Lukas: v1.1.3 added smooth option to input, made plotCol,
%                           meanline and smoothplot optional inputs
%% check inputs
%set defaults
plotCol = [0 0 0];
meanline = '-';
smoothPlot = true;
if nargin>5
    for i=1:2:length(varargin)
        if isa(varargin{i},'char') || isa(varargin{i},'string')
            switch lower(varargin{i})
                case 'plotcol'
                    plotCol = varargin{i+1};
                case 'meanline'
                    meanline = varargin{i+1};
                case 'smoothplot'
                    %set to ~ because nanMeanSterrHistogram uses dontsmooth
                    smoothPlot = ~varargin{i+1};
                otherwise
                    error(['unknown input parameter: ' varargin{i+1}])
            end
        else
            error(['unknown input parameter: ' varargin{i+1}])
        end
    end
    %%
    fireRate = calcFireRate(spikes,stepSz,binSize,interval,smpRate);
    if size(fireRate{1,1},1)<2 % if less than two trials have spikes
        fireRate{1,1} = nan(2,sum(abs(interval))/stepSz);
        warning(...
            'No trial with spikes detected, firing rate was set to NaN.')
    end
    for intervalNo=1:size(interval,1)
        [trlmean, trlsterr, ~] = nanMeanSterrHistogram(...
            fireRate{intervalNo},...
            'errcol',plotCol,'meancol',plotCol,...
            'nomean',false,'erralpha',0.7,'meansize',3,...
            'dontsmooth',smoothPlot,'meanline',meanline);
    end
    upperY = ceil(max(trlmean+trlsterr))+1;
end