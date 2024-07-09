%% NGL03_plotting
% To run after all standard data have been sorted, curated and saved. In
% principle, different tyoes of plots could be selected to be done or not,
% and new plots could be added for personalization. Parameters are taken
% from main script or from an additional one.
%
% Jesus 04.07.2024

%% 00. Check current inputs.
% Check if input variable exist already. Parse values.
if ~exist("input","var")
    input = struct( 'datadrive' , datadrive , ...   % force char array
                    'studyName' , studyname , ...   % force char array
                    'toolbox'   , toolbox   , ...   % force char array
                    'subjects'  , [], ...           % do NOT force char array
                    'dates'     , []        );      % do NOT force char array
    input.dates     = dates;    % place as it comes
    input.subjects  = subjects; % place as it comes
else
    disp('Using INPUTS from NGL01_MAIN.')
end

% Set default inputs and dependencies. In case NGL01 did not before.
input = set_default(input, opt);

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        input.run = [x y]; % Current run, to pass to functions.
        
        %% 02. Prepare to proceed with a single session.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);           
        mkdir(fullfile(opt.analysis,'genplots'))

        %% 03. Recover event data if not in workspace yet
        % Recover trial definitions created after event extraction and processing. 
        % Will have as many variations as requested at that time. Needs to be ran 
        % again to create new alignments.

        if ~exist('events','var'),    load(fullfile(opt.analysis, "events.mat")),      end
        if ~exist('neurons','var'),   load(fullfile(opt.analysis, "neurons.mat")),     end
%         if ~exist('trialdef','var'),  load(fullfile(opt.trialSorted, "trialdef.mat")), end
%         if ~exist('condition','var'), load(fullfile(opt.analysis, "condition.mat")),   end 
        
        %% 04.1 All clusters piled, ...
        % aligned to requested events, for all clusters
        if opt.pooledstats
            plot_pooledstats(neurons, opt, param)
        end
                   
        %% 04.2 Raster, trace, PSH, ISI
        % whole trial
        if opt.plot_trial
            if ~exist('spike','var'),     load(fullfile(opt.spikeSorted, "spike.mat")),    end
            plot_trialclusters(neurons, events, spike, opt, param)
        end

        % aligned to requested events, per cluster
        if opt.plot_align && numel(opt.alignto) > 1
            if ~exist('spike','var'),     load(fullfile(opt.spikeSorted, "spike.mat")),    end
            plot_alignedclusters(neurons, events, spike, opt, param)
        end

        clear events neurons spike
    end
end