function [fireRate, normFireRate] = calcFireRate(spikes, par, varargin)
% Use this function to calculate the firing rate of a neuron over the time
% course of the trial.
%
% INPUTS
% 'spikes'   : cell containing aligned spike times per trial
%
% 'par' structure with fields:
%  '.stepSz'          : step size of for the walk through the interval of
%                        interest (in samples)
%  '.binSize'         : size of the bins in which the spike train will be
%                        divided (in samples)
%  '.interval'        : matrix containing the start and the end point(s)
%                        (in samples) of the interval of interest, relative
%                        to the alignment
%  '.smpRate'         : sampling rate of the recording (in Hz)
%
% OUTPUTS
% 'fireRate'    : cell containing absolute counts of spiking in each
%                  trial (for the selected step through the interval),
%                  1 cell per interval.
% 'normFireRate': cell containing dimensionless value of spiking in each
%                  trial, 1 cell per trial.
% VERSION HISTORY:
% Author:         Lukas Hahn
% Version:        1.0.3
%
% 15.07.2019, Lukas: v1.0.0 release version
% 12.12.2023, Lukas: v1.0.1 updated documentation
% 10.04.2024, Lukas: v1.0.2 bug fix: loop index 'int' is now based on 
%                           size of correct windowBorder dimension (2)
% 03.03.2024, Jesus: v1.0.3 modified input for simplicity. Added
%                           normalization over baseline optional input, to obtail trial-long
%                           dimensionless fire rate change.

%% Defaults
if isempty(par)
    par.stepSz  = 50;
    par.binSize = 500;
    par.interval= [-1000 2000];
    par.smpRate = 1000;  
end
baseBins = 0;

%% varargin
if nargin > 2
   for i = 1:2:length(varargin)
       if isa(varargin{i},'char') || isa(varargin{i},'string')
          switch lower(varargin{i})
                case 'baseline'
                    baseBins = varargin{i+1}/par.stepSz;
                otherwise
                    error(['unknown input parameter: ' varargin{i+1}])
          end
       else
           error(['unknown input parameter: ' varargin{i+1}])
       end
   end
   else
      %use defaults
end

%% Define interval borders
windowBorder = cell(1,size(par.interval,1));
for i=1:size(par.interval,1) %for all intervals
    windowBorder{1,i} = par.interval(i,1):par.stepSz:par.interval(i,2);
end

%% Firing rate per trial with moving window
fireRate = cell(size(spikes,2),1);
for int = 1:size(windowBorder,2)
    for trl = 1:size(spikes,1) %for all trials
        % for all steps within the borders of the selected interval
        for bin = 1:length(windowBorder{1,int})-1
            % absolute number of spikes within the selected portion
            fireRate{int}(trl,bin) = sum(...
                spikes{trl}>=windowBorder{1,int}(bin) &...
                spikes{trl}<(windowBorder{1,int}(bin)+par.binSize))...
                *(par.smpRate/par.binSize); %transformation to per second
        end
    end
end

% Normalized firing rate for the whole period
normFireRate = [];
if baseBins > 0
    for int = 1:size(windowBorder,2)
        if ~isempty(fireRate{int})
            % Compute the baseline (2:baseBins) normalization 
            [~, C, S] = normalize(mean(fireRate{int}(:,1:baseBins),1));
            
            % Normalize the entire time series
            normFireRate = normalize(fireRate{int}, "center", C, "scale", S);
            normFireRate = mean(normFireRate,1);
        else
            normFireRate = nan(length(windowBorder{1,int})-1,1);
        end
    end
end

end