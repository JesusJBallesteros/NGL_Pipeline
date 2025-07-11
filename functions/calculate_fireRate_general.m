function fireRate = calculate_fireRate_general(neurons, events, conditions, opt, param)
% Jesus. 06.03.2025.
% Description
% INPUTS
% OUTPUTS

%% Default options.
if ~isfield(param,'trial2plot'),    param.trial2plot     = 'allInitiated'; end % condition to plot
if ~isfield(param,'binSize'),       param.binSize        = 200;          end % ms
if ~isfield(param,'stepSz'),        param.stepSz         = 20;           end % ms
if ~isfield(param,'interval'),      param.interval       = [-2000 10000];end % ms range
if ~isfield(param,'smpRate'),       param.smpRate        = 1000;         end % smp/s
if ~isfield(param,'baseline'),      param.baseline       = -param.interval(1);end % ms
if ~isfield(param,'plot'),          param.plot           = true;         end % logic
if ~isfield(param,'blockchange'),   param.blockchange    = [];           end % trial number

%% Initialize
toalignto = opt.alignto;

%% Prepare treatments
param.levels = 1;
param.blockchange = (find(diff(conditions.block)>0)+1)';

% % E.G % add levels accordingly, e.g. basal/treatment_present/post (+2) or basal/post (+1)
% param.levels{1} = numel(events.(opt.trEvents{1}).trial{1});
% param.levels{2} = 2;
% 
% % For blocks, add the first and last trials to complete 8 blocks
% if events.(opt.trEvents{1}).trial{1}(1) == 0
%     param.trial_change{1} = [1 events.(opt.trEvents{1}).trial{1}(2:end) length(neurons.itiOn{1})];
% else
%     param.trial_change{1} = [1 events.(opt.trEvents{1}).trial{1}(1:end) length(neurons.itiOn{1})];
% end
% 
% param.trial_change{2} = events.(opt.trEvents{2}).trial{1}(1:end);

%% Trial indexing
for a = 1 % align only to ini (SUBJECT TO CHANGE AND LOOP several)
    neuronSet = neurons.(toalignto{a});
    % trialrange = 1:length(conditions.correct); % we could set the trial range (not INDEX)

    % For each cluster
    for c = 1:length(neuronSet)
       %% Cluster's spike set selection
       toCalculate = cell(1,param.levels); % preallocation
       for p = 1:param.levels 
            % Initial, full set for any level
            toCalculate{p} = neuronSet{c}; % actual set for this level

            % Find empty trials before any indexig (no activity at all)
            emptytrials = cellfun(@isempty, toCalculate{p});
            toCalculate{p}(emptytrials) = {NaN}; % make those Nan

            % Index any valid/non-valid conditions
            cndidx = ones(1,length(conditions.correct)); % all valid
            % E.G % if p == 2
            %       cndidx = events.(param.cond).trial{1}; % only a trial type 
            %       end
            
            % Empty non-valid
            toCalculate{p}(~cndidx) = {[]}; % empty complementary
                
            % Plus, empty specific conditions based on behavior
            if strcmp(param.trial2plot, 'allInitiated')
                toCalculate{p}(logical(conditions.aborted)) = {[]};
            else
                toCalculate{p}(logical(~conditions.(param.trial2plot))) = {[]};
            end
        end
    
        %% FireRate calculation
        % The fire rate is per (c)luster per level (p) of analisys (blocks/treatments), and per alignment
        for p = 1:length(toCalculate)
            param.cl = [a c p]; % pass cluster, and level (and aligment and...)
            spikes2use = ~cellfun(@isempty, toCalculate{p});
            % spikes2use = true(1,size(toCalculate{i},1));

            if ~isempty(param.blockchange) && c == 1
                for b = 1:size(param.blockchange,2)
                    substr = sum(spikes2use(1:param.blockchange(1,b))==0);
                    param.blockchange(1,b) = param.blockchange(1,b)-substr;
                end
            end

            [fireRate.sps{c,p}, fireRate.Norm{c,p}, fireRate.meanNorm{c,p}] = calcFireRate(toCalculate{p}(spikes2use), opt, param);

            if iscell(fireRate.sps{c,p}), fireRate.sps{c,p} = cell2mat(fireRate.sps{c,p}); end
            if iscell(fireRate.Norm{c,p}), fireRate.Norm{c,p} = cell2mat(fireRate.Norm{c,p}); end
            if iscell(fireRate.meanNorm{c,p}), fireRate.meanNorm{c,p} = cell2mat(fireRate.meanNorm{c,p}); end

        end
    end
end

%% Plot all clusters for session
if param.plot
    plot_multi_fireRate(fireRate.meanNorm, param, opt)
end

end % function end