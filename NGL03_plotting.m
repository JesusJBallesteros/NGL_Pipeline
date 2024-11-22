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
dosave = 1;

for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        input.run = [x y]; % Current run, to pass to functions.
        
        %% 02. Proceed only once
        if exist(fullfile(input.analysis, "SocialLearning_all.mat"),"file")
            load(fullfile(input.analysis, "SocialLearning_all.mat"));
            dosave = 0;
            break
        end
        
        %% 03. Prepare to proceed with a single session.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);           
        if ~exist(fullfile(opt.analysis,'newplots'),"dir")
            mkdir(fullfile(opt.analysis,'newplots'))
        end

        %% 04. Recover all Project data if not collected yet
        % Recover trial definitions created after event extraction and processing. 
        % Will have as many variations as requested at that time. Needs to be ran 
        % again to create new alignments.
        if ~exist('events','var'),    load(fullfile(opt.analysis, "events.mat")),      end
        if ~exist('neurons','var'),   load(fullfile(opt.analysis, "neurons.mat")),     end
        if ~exist('condition','var'), load(fullfile(opt.analysis, "condition.mat")),   end 
        if ~exist('spike','var'),     load(fullfile(opt.spikeSorted, "spike.mat")),    end
        if ~exist('blob','var') && opt.plot_SocLear
            load(fullfile(opt.analysis, "blob.mat")),
            blob = rmfield(blob,{'IDs' 'Pos' 'Par' 'Tracking'});
        end
        
        % Collect all data into single all variables: events, neurons, conditions, blob
        allneurons{x,y}     = neurons;
        allevents{x,y}      = events;
        allconditions{x,y}  = conditions;
        allspike{x,y}       = spike;
        allblobs{x,y}       = blob;
    
        clear events neurons spike blob spike
    end
end

%% 05. Save it for next iterations
if dosave
    save(fullfile(input.analysis, "SocialLearning_all.mat"), "allblobs","allspike","allconditions","allevents","allneurons")
end

%% Plot session by session
for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        input.run = [x y]; % Current run, to pass to functions.

        %% 04.1 All clusters piled, ...
        % aligned to requested events, for all clusters
        if opt.pooledstats
            plot_pooledstats(allneurons{x,y}, opt, param)
        end
                   
        %% 04.2 Raster, trace, PSH, ISI
        % whole trial
        if opt.plot_trial
            if ~exist('spike','var'),     load(fullfile(opt.spikeSorted, "spike.mat")),    end
            plot_trialclusters(allneurons{x,y}, allevents{x,y}, spike, allconditions{x,y}, opt, param)
        end
        
        % aligned to requested events, per cluster
        if opt.plot_align && numel(opt.alignto) > 1
            plot_alignedclusters(allneurons{x,y}, allevents{x,y}, spike, [], opt, param)
        end
        
        %% 04.03 SocialLearning-specific plots
        if opt.plot_SocLear
            plot_SocLear(allneurons{x,y}, allevents{x,y}, allspike{x,y}, allconditions{x,y}, allblobs{x,y}, opt, input, param)
        end

    end
end