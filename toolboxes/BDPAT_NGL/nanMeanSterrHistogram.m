function [trlmean, trlsterr, nanWarn] = nanMeanSterrHistogram(dat,varargin)
%% Plots histograms as mean with standard error. Uses only not NAN entries
% based om meanSterrHistogram.m
%
% Uses nanmean and nanste to ignote NaNs. If means / ste still contain NaN,
%   will be set to zero and a warning added to the plot.
%
% Inputs:
%   - dat       : matrix of binned data : dat(trials,bins)
% Inputs-optional:
%   - varargin  : parameter-value pairs
%               : 'Winsze'      size of the sliding window in bins
%               : 'meancol'     color of the histogram
%               : 'meansize'    linesize of the histogram
%               : 'errstyle'    1: error-lines 0: error as shades
%               : 'errcol'      color of the errorplot
%               : 'errline'     linestyle of errorlines
%               : 'nomean'      plot only sterr, no mean
%               : 'nanwarn'     warn if NaNs are replaced by zeros
%               : 'suprPlot'    supress plotting functions
% Returns:
%   - trlmean   : mean accross all trials, each bin
%   - trlsterr  : standard error accross all trial, eachbin
%   - nanWarn   : a flag, true if NaN were encoutered in mean/ ste

% Author:         Jonas
% Version:        1.6
% Last Change:    03/03/25

% 1.0 : first draft
% 1.2 : added error aplpha
% 1.3 : added ouput
% 1.4 : catching NaN, warning of corretion
% 1.5 : changed dimensions of trlmean / nanmean if dontsmooth is true
% 1.6 : added linestyle choice for mean 'meanline'(Lukas)
% 1.7: Jesus, added plot supression to keep calculations as done before but
       % without the need of graphical output, saving time

%%  BUGS / TODO:

%% Default values, use optional string arguments to set
winsze      = 5;                % size of the moving average-widow in bins
errcol      = [.8 .8 .8];       % 'k' for lines grey shade color
errAlpha    = 1;       			% default error SHADE has alpha 1
errstyle    = 0;                % plot sterr as lines (1) or as shades (0)
errline     = ':';
meanline    = '-';              %addition here
meancol     = 'b';
meansize    = 1;
nomean      = false;
smoother    = 'moving';
nosmooth    = false;
nanWarn     = false;
doNanWarn   = false;
suprPlot    = false; % Added by Jesus to supress plot and maintain calculations

%% get the variable input arguments
if nargin > 1
    for i=1:length(varargin)
        if isa(varargin{i},'char')
            switch lower(varargin{i})
                case 'winsze'
                    winsze      = varargin{i+1};
                case 'meancol'
                    meancol     = varargin{i+1};
                case 'meansize'
                    meansize    = varargin{i+1};
                case 'errstyle'
                    errstyle    = varargin{i+1};
                case 'errcol'
                    errcol      = varargin{i+1};
                case 'erralpha'
                    errAlpha     = varargin{i+1};
                case 'errline'
                    errline     = varargin{i+1};
                case 'meanline'                     %addition here
                    meanline    = varargin{i+1};
                case 'nomean'
                    nomean      = varargin{i+1};
                case 'smoother'
                    smoother    = varargin{i+1};
                case 'dontsmooth'
                    nosmooth    = varargin{i+1};
                case 'nanwarn'
                    doNanWarn   = varargin{i+1};
                case 'suprplot'
                    suprPlot   = varargin{i+1};
            end
        end
    end
end

% % smoothen trial by trial (DEPRECATED)
% trlLen = size(dat,2);
% for t=1:size(dat,1)
%     % smooth data, including a bit of iti in the end to avoid edge-effects
%     tmp         = smooth([dat(t,:) dat(t,1:winsze)],winsze,smoother);
%     dat(t,:)    = tmp(1:trlLen);
% end
%
% trlmean     = mean(dat);
% trlsterr    = ste(dat);

%% Mean/StErr calculation
if ~nosmooth
    % calculate mean, sterr and smoothen
    trlmean     = smooth(nanmean(dat),winsze,smoother)';
    trlsterr    = smooth(nanste(dat),winsze,smoother)';
else
    % calculate mean, sterr
    trlmean     = nanmean(dat);
    trlsterr    = nanste(dat); %deleted ' to fix problem (Lukas 19.12.2016)
end

%% Plotting
if ~suprPlot
    % % make axes current
    % axes(axHand);
    hold on;
    
    % plot standard error
    if errstyle
        plot(trlmean-trlsterr,'Color',errcol,'LineStyle',errline);
        plot(trlmean+trlsterr,'Color',errcol,'LineStyle',errline);
    else
        back = trlmean-trlsterr;
        back = back(end:-1:1);
        if any(isnan(back)) || any(isnan(trlmean)) || any(isnan(trlsterr))
            back(isnan(back))           = 0;
            trlmean(isnan(trlmean))     = 0;
            trlsterr(isnan(trlsterr))   = 0;
            nanWarn = true;
        end
        hand = fill([1:length(trlmean), length(trlmean):-1:1],...
            [trlmean+trlsterr, back],errcol,...
            'Linestyle','none','edgecolor','none');
        alpha(hand,errAlpha);
        if doNanWarn && nanWarn
            text(1, max(trlmean+trlsterr)- max(trlmean+trlsterr)*.1,...
                'set NaN to zero','color',[1 0 0]);
        end
    end
    
    if ~nomean
        % plot the mean-histogram
        plot(trlmean,'Color',meancol,'Linewidth',meansize,'LineStyle',meanline); %change here
    end
end

end