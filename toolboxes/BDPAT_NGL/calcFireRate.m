function [fireRate, normFireRate, meanNormFireRate] = calcFireRate(spikes, opt, param)
% Use this function to calculate the firing rate of a neuron over the time
% course of the trial.
%
% INPUTS
% 'spikes'   : cell containing aligned spike times per trial
%
% 'opt'     : contains options to plot and save graphs
%
% 'param' structure with fields:
%  '.stepSz'          : step size of for the walk through the interval of
%                        interest (in samples)
%  '.binSize'         : size of the bins in which the spike train will be
%                        divided (in samples)
%  '.interval'        : matrix containing the start and the end point(s)
%                        (in samples) of the interval of interest, relative
%                        to the alignment
%  '.smpRate'         : sampling rate of the recording (in Hz)
%
% 'varargin' OPTIONAL include pairwise inpus as:
%   'baseline'        : followed by the number of miliseconds to use as
%                       baseline period
%   'plot'            : true or false
%
% OUTPUTS
% 'fireRate'    : cell containing absolute counts of spiking in each
%                  trial (for the selected step through the interval),
%                  1 cell per interval.
% 'normFireRate': cell containing dimensionless value of spiking in each
%                  trial, 1 cell per trial.
% VERSION HISTORY:
% Author:         Lukas Hahn
%
% 15.07.2019, Lukas: v1.0.0 release version
% 12.12.2023, Lukas: v1.0.1 updated documentation
% 10.04.2024, Lukas: v1.0.2 bug fix: loop index 'int' is now based on 
%                           size of correct windowBorder dimension (2)
% 03.03.2024, Jesus: v1.0.3 modified input for simplicity. Added normalization
%                           over baseline optional input.
% 14.03.2024, Jesus: v1.1   Added possibility to plot fr and Norm fr

%% Defaults
if isempty(param)
    param.stepSz  = 200;
    param.binSize = 400;
    param.interval= [-1000 3000];
    param.smpRate = 1000;  
    param.baseline  = [];
    param.plot  = false;
    param.blockchange = [];
end

% To find bin corresponding to alignment (time 0)
param.xtick = (0:(0-param.interval(1)):diff([param.interval(1) param.interval(2)]))/param.stepSz; % time ticks
param.xticklabels = {((param.xtick*param.stepSz)+param.interval(1))/1000};
param.inibin = param.xtick(find(param.xticklabels{1}==0));

%% Define interval borders
windowBorder = cell(1,size(param.interval,1));
for i=1:size(param.interval,1) %for all intervals
    windowBorder{1,i} = param.interval(i,1):param.stepSz:param.interval(i,2);
end

%% Firing rate per trial with moving window
fireRate = cell(size(spikes,2),1);
for int = 1:size(windowBorder,2)
    for trl = 1:size(spikes,1) %for all trials
        % for all steps within the borders of the selected interval
        if any(isnan(spikes{trl}))
            fireRate{int}(trl,:) = nan(1, length(windowBorder{1,int})-1);
        else
            for bin = 1:length(windowBorder{1,int})-1
                % absolute number of spikes within the selected portion
                fireRate{int}(trl,bin) = sum(...
                    spikes{trl}>=windowBorder{1,int}(bin) &...
                    spikes{trl}<(windowBorder{1,int}(bin)+param.binSize))...
                    *(param.smpRate/param.binSize); %transformation to per second
            end
        end
    end
end

%% Normalized firing rate for the whole period, if baseline specified
normFireRate = [];
if ~isempty(param.baseline)
    baseBins = param.baseline(1)/param.stepSz-1;
    for int = 1:size(windowBorder,2)
        if ~isempty(fireRate{int})
            
            % Use of single block/phase to calculate baseline
            if size(param.baseline, 2) > 1
                if param.baseline(2) == 1, basetrials = 1:param.blockchange(param.baseline(2));
                else,                      basetrials = param.blockchange(param.baseline(2)-1):param.blockchange(param.baseline(2));
                end
            else,  basetrials = 1:length(fireRate{int});
            end
            
            % Compute the baseline for normalization 
            [~, C, S] = normalize(fireRate{int}(basetrials,param.inibin-baseBins:floor(param.inibin)), 1, "zscore", "std");
            
            % Normalize the entire time series
            normFireRate = normalize(fireRate{int}, "center", mean(C), "scale", mean(S));
            
            % Get average normalize firing rate
            meanNormFireRate = mean(normFireRate,1,"omitnan");
        else
            meanNormFireRate = nan(length(windowBorder{1,int})-1,1);
        end
    end
end

%% If plot requested, plot per cluster and per session
if param.plot
    plot_single_fireRate(fireRate{1,1}, normFireRate, param, opt)
end

end