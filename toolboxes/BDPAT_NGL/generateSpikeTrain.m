function spkPerTrial = generateSpikeTrain(spksPerSec,binSize,varargin)
%%spkPerTrial = generateSpikeTrain(spksPerSec,binSize,varargin)
%
% This function creates a spike train for an arbitrary amount of bins and
% trials. The function automatically assumes a 2 ms refractory period
% around each spike, i.e. two consecutive spikes in the resulting spike
% train are separated by at least 3 ms.
%
%INPUTS
%  * 'spksPerSec'         : matrix containing number of spikes per bin
%                           (bins rows, repetitions columns)
%  * 'binSize'            : temporal interval of a bin, in ms
%
%OPTIONAL INPUTS
% varargin                : parameter-value pairs
%  * 'conditionindex'     : index of condition
%  * 'stimcondition'      : stimulation conditions (vector with true or
%                           false per condition)
%  * 'stimFactor'         : factor of stimulation (vector with one value
%                           per condition)
%  * 'stimStart'          : first bin of stimulation (vector with one value
%                           per condition)
%  * 'stimend'            : last bin of stimulation (vector with one value
%                           per condition)
%  * 'stimMin'            : start of stimulation interval within a bin in
%                           ms (vector with one value per condition)
%  * 'stimMax'            : end of stimulation interval within a bin in ms
%                           (vector with one value per condition)
%
%OUTPUTS
%   * 'spkPerTrial'       : cell (row repetitions) containing spike trains

% VERSION HISTORY:
% Author:         Lukas Hahn
% Version:        1.0.0
% Last Change:    25.01.2024
%
% 25.01.2024, Lukas: v1.0.0 release version
%% check inputs
%set working default values
cI = 1;
stimCond = 0;
stimulationFactor = 0;
stimBinS = 0;
stimBinE = 0;
stimMin = 0;
stimMax = binSize;
if nargin>2 %take provided values
    for a=1:size(varargin,2)
        if ischar(varargin{a})
            switch lower(varargin{a})
                case 'conditionindex'
                    cI = varargin{a+1};
                case 'stimcondition'
                    stimCond = varargin{a+1};
                case 'stimfactor'
                    stimulationFactor = varargin{a+1};
                case 'stimstart'
                    stimBinS = varargin{a+1};
                case 'stimend'
                    stimBinE = varargin{a+1};
                case 'stimmin'
                    stimMin = varargin{a+1};
                case 'stimmax'
                    stimMax = varargin{a+1};
                otherwise
                    error(['unknown parameter: ' varargin{a}])
            end
        end
    end
end
%% run spike train generation
rng('shuffle')
spkPerTrial = cell(size(spksPerSec,2),1);
for trlNo=1:size(spksPerSec,2) %for repetitions (trials)
    spkNo = 0;
    for binNo=1:size(spksPerSec,1) %for bins within repetition
        if spksPerSec(binNo,trlNo)>0
            %determines if a spike is a 'baseline spike'
            % or a 'stimulation spike',
            %relevant to achieve distinct intervals of
            % higher spiking while maintaining regular
            %spiking in surrounding intervals
            spkGrp = rand(spksPerSec(binNo,trlNo),1);
            for spk=1:spksPerSec(binNo,trlNo) %for all spikes
                spkNo = spkNo+1;

                %spike timepoint in ms

                %to account for refractory time
                % (+/- 2 ms around spike event)
                % determine blocked time points
                blockedSpkTimes = [...
                    spkPerTrial{trlNo,1}; ...
                    spkPerTrial{trlNo,1}-1; ...
                    spkPerTrial{trlNo,1}-2; ...
                    spkPerTrial{trlNo,1}+1; ...
                    spkPerTrial{trlNo,1}+2];

                %draw up random order of spike times within
                %interval
                allSpkTimes = randperm(binSize)+...
                    ((binNo-1)*binSize);

                if spkGrp(spk)>=(1/stimulationFactor(cI)) && ...
                        ismember(binNo,...
                        stimBinS(cI):stimBinE(cI)) && ...
                        stimCond(cI)==1
                    % if stimulation

                    %limit spike times to pre-defined interval
                    allSpkTimes(~ismember(allSpkTimes,...
                        (stimMin(cI):stimMax(cI))+...
                        ((binNo-1)*binSize))) = [];
                else %normal spiking
                    % no changes required
                end

                %remove blocked spike times
                allSpkTimes(ismember(...
                    allSpkTimes,blockedSpkTimes)) = [];

                if ~isempty(allSpkTimes)
                    %set spike time
                    spkPerTrial{trlNo,1}(spkNo,1) = allSpkTimes(1);
                else
                    %no possible time point left to place
                    % another spike into the bin
                    warning(['bin is full, no timepoints left to place '...
                        'a spike'])
                end
            end
        end
    end
end
end