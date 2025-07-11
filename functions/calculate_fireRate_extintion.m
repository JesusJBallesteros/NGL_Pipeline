function fireRate = calculate_fireRate_extintion(neurons, events, conditions, opt, param)
% Jesus. 06.03.2025.
% Description
% INPUTS
% OUTPUTS

%% Default options.
if ~isfield(param,'nBlocks'),       param.nBlocks        = 8;            end
if ~isfield(param,'IncludeFS'),     param.IncludeFS      = false;        end
if ~isfield(param,'trial2plot'),    param.trial2plot     = 'correct';    end
if ~isfield(param,'binSize'),       param.binSize        = 500;          end
if ~isfield(param,'stepSz'),        param.stepSz         = 50;           end
if ~isfield(param,'interval'),      param.interval       = [-2000 10000];end
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end
if ~isfield(param,'baseline'),      param.baseline       = 2000;         end
if ~isfield(param,'plot'),          param.plot           = false;         end
    param.xtick = (0:2000:diff([param.interval(1) param.interval(2)]))/param.stepSz; % time ticks
    param.xticklabels = {'-2', 'Ini','StimOn', '4', '6', '8', '10'};
    param.timelabel = 'time (s)';

%% Initialize
toalignto = opt.alignto;
param.trial_change{1}   = 1; % force cell integer

%% Prepare treatments
% add levels accordingly, e.g. basal/treatment_present/post (+2) or basal/post (+1)
param.levels{1} = numel(events.(opt.trEvents{1}).trial{1});
param.levels{2} = 2;

% For blocks, add the first and last trials to complete the 8 blocks
if events.(opt.trEvents{1}).trial{1}(1) == 0
    param.trial_change{1} = [1 events.(opt.trEvents{1}).trial{1}(2:end) length(neurons.itiOn{1})];
else
    param.trial_change{1} = [1 events.(opt.trEvents{1}).trial{1}(1:end) length(neurons.itiOn{1})];
end

param.trial_change{2} = events.(opt.trEvents{2}).trial{1}(1:end);

%% fireRate full trial
for a = 1 % align to ini
    neuronSet = neurons.(toalignto{a});
    trialrange = 1:length(conditions.correct);

    % For each cluster
    for c = 1:length(neuronSet)
        % Do this for each requested situation (#levels)
        for p = 1:length(param.levels)
            %% Cluster's spike set selection
            if param.levels{p} > 2
                % Empy array for this cluster
                toCalculate = cell(1, param.nBlocks);

                for b = 1:param.nBlocks
                    % Take full set for the block
                    toCalculate{b} = neuronSet{c};
                
                    % Here, we find trials for the current block range
                    cndidx = ismember(trialrange, param.trial_change{1}(b):param.trial_change{1}(b+1));
                    
                    % And we keep ONLY those for Novel Stimuli (tr2)
                    if ~param.IncludeFS, FSidx = ismember(trialrange, events.tr2.trial{1}); 
                    else,                FSidx = ones(1,length(trialrange)); end
                    cndidx = trialrange(FSidx & cndidx);
                    
                    % We find the complementary set of trials (non-block (and FS, if so))
                    notcndidx = trialrange(setdiff(1:end,cndidx));
                    
                    % Remove the complementary set
                    toCalculate{b}(notcndidx) = {[]};
                    
                    % Within interest block set, empty trials of conditions of no interest
                    if strcmp(param.trial2plot, 'correct')
                        toCalculate{b}(conditions.aborted | conditions.omission | conditions.incorrect) = {[]};
                    elseif strcmp(param.trial2plot, 'incorrect')
                        toCalculate{b}(conditions.aborted | conditions.omission | conditions.correct) = {[]};
                    elseif strcmp(param.trial2plot, 'omission')
                        toCalculate{b}(conditions.aborted | conditions.correct | conditions.incorrect) = {[]};
                    elseif strcmp(param.trial2plot, 'allInitiated')
                        toCalculate{b}(logical(conditions.aborted)) = {[]};
                    end
                end

            elseif param.levels{p} == 2
                param.cond = 'tr2';
                % Empy array for this cluster
                toCalculate = cell(1,2);
    
                % Initial, full set for both levels 
                toCalculate{1} = neuronSet{c};
                toCalculate{2} = neuronSet{c};
    
                % Here, we index the tr2 code (13) for NS
                cndidx = events.(param.cond).trial{1};
                
                % find the complementary set of trials, the FS
                notcndidx = trialrange(setdiff(1:end,cndidx));
    
                if ~isempty(cndidx)
                    % Empty FS trials
                    toCalculate{1}(notcndidx) = {[]};
    
                    % Empty NS trials 
                    toCalculate{2}(cndidx) = {[]};
                    
                    % Also, empty Aborted and omited trials
                    if strcmp(param.trial2plot, 'correct')
                        toCalculate{1}(conditions.aborted | conditions.omission | conditions.incorrect) = {[]};
                        toCalculate{2}(conditions.aborted | conditions.omission | conditions.incorrect) = {[]};
                        % param.ext = '_Correct.png';
                    elseif strcmp(param.trial2plot, 'incorrect')
                        toCalculate{1}(conditions.aborted | conditions.omission | conditions.correct) = {[]};
                        toCalculate{2}(conditions.aborted | conditions.omission | conditions.correct) = {[]};
                        % param.ext = '_Incorrect.png';
                    elseif strcmp(param.trial2plot, 'omission')
                        toCalculate{1}(conditions.aborted | conditions.correct | conditions.incorrect) = {[]};
                        toCalculate{2}(conditions.aborted | conditions.correct | conditions.incorrect) = {[]};
                        % param.ext = '_Omission.png'; 
                    elseif strcmp(param.trial2plot, 'allInitiated')
                        toCalculate{1}(logical(conditions.aborted)) = {[]};
                        toCalculate{2}(logical(conditions.aborted)) = {[]};
                    end
                end
            end
        
            %% FireRate calculation
            % The fire rate is per (c)luster and per level (p) of analisys (blocks/treatments)
            for i = 1:length(toCalculate)
                spikes2use = ~cellfun(@isempty, toCalculate{i});
                [fireRate.sps{c,p}(i), ~, fireRate.NormMean{c,p}(i,:)] = calcFireRate(toCalculate{i}(spikes2use), opt, param);
            end
        end
    end
end

end % function end