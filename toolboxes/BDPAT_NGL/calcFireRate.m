function fireRate = calcFireRate(alignedSpikes,stepSz,binSize,interval,...
    smpRate)
%%
%function fireRate = calcFireRate(alignedSpikes,stepSz,binSize,interval,...
% smpRate)
%
% Use this function to calculate the firing rate of a neuron over the time
% course of the trial.
%
%INPUTS
%  * 'alignedSpikes'   : cell containing aligned spike times per trial
%  * 'stepSz'          : step size of for the walk through the interval of
%                        interest (in samples)
%  * 'binSize'         : size of the bins in which the spike train will be
%                        divided (in samples)
%  * 'interval'        : matrix containing the start and the end point(s)
%                        (in samples) of the interval of interest, relative
%                        to the alignment
%  * 'smpRate'         : sampling rate of the recording (in Hz)
%
%OUTPUTS
%   * 'fireRate'   : cell containing absolute counts of spiking in each
%                    trial (for the selected step through the interval),
%                    1 cell per interval.

% VERSION HISTORY:
% Author:         Lukas Hahn
% Version:        1.0.2
% Last Change:    10.04.2024
%
% 15.07.2019, Lukas: v1.0.0 release version
% 12.12.2023, Lukas: v1.0.1 updated documentation
% 10.04.2024, Lukas: v1.0.2 bug fix: loop index 'int' is now based on 
%                           size of correct windowBorder dimension (2)
%%
windowBorder = cell(1,size(interval,1));
for i=1:size(interval,1) %for all intervals
    windowBorder{1,i} = interval(i,1):stepSz:interval(i,2);
end
%firing rate during the selected interval
fireRate = cell(size(alignedSpikes,2),1);
for int=1:size(windowBorder,2)
    for trl=1:size(alignedSpikes,1) %for all trials
        %for all steps within the borders of the selected interval
        for bin=1:length(windowBorder{1,int})-1
            %absolute number of spikes within the selected portion
            fireRate{int}(trl,bin) = sum(...
                alignedSpikes{trl}>=windowBorder{1,int}(bin) &...
                alignedSpikes{trl}<(windowBorder{1,int}(bin)+binSize))...
                *(smpRate/binSize); %transformation to per second
        end
    end
end