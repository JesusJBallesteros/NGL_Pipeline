%% NGL03_NEWplotting
% Jesus 05.03.2025

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
if ~isfield(opt,'plot_SocLear'), opt.plot_SocLear = false; end
collect = 0;

%% 01. Recover all Project data if not collected yet
input.sessions = findSessions(input);
try load(fullfile(input.analysis, "data_all.mat"));
catch, collect = 1;
end

if collect
   for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        input.run = [x y]; % Current run, to pass to functions.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt); 

        % Recover trial definitions created after event extraction and processing. 
        % Will have as many variations as requested at that time. Needs to be ran 
        % again to create new alignments.
        if ~exist('events','var'),    load(fullfile(opt.analysis, "events.mat")),      end
        if ~exist('neurons','var'),   load(fullfile(opt.analysis, "neurons.mat")),     end
        if ~exist('condition','var'), load(fullfile(opt.analysis, "condition.mat")),   end 
        if ~exist('fireRate','var'),  load(fullfile(opt.analysis, "fireRate.mat")),    end
        if ~exist('fireRateNorm','var'), load(fullfile(opt.analysis, "fireRateNorm.mat")), end

        if ~exist('spike','var'),     load(fullfile(opt.spikeSorted, "spike.mat")),    end

        if ~exist('blob','var')     && opt.plot_SocLear
            load(fullfile(opt.analysis, "blob.mat")),
        end
        
        % Collect all data into single all variables: events, neurons, conditions, blob
        allneurons{x,y}     = neurons;
        allevents{x,y}      = events;
        allconditions{x,y}  = conditions;
        allfireRate{x,y}    = fireRate;
        allfireRateNorm{x,y}= fireRateNorm;

        allspike{x,y}       = spike;

        if opt.plot_SocLear
           allblobs{x,y}    = blob;
        end
        clear events neurons spike blob spike fireRate fireRateNorm conditions
    end
   end

% Save collected data
   if opt.plot_SocLear
    save(fullfile(input.analysis, "data_all.mat"),"allspike","allconditions","allevents","allneurons","allfireRate","allfireRateNorm","allblobs");
   else
    save(fullfile(input.analysis, "data_all.mat"),"allspike","allconditions","allevents","allfireRate","allfireRateNorm","allneurons");
   end
end

%% 02. Plot session by session
for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        input.run = [x y]; % Current run, to pass to functions.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt); 
        
        if ~exist(fullfile(opt.analysis,'newplots'),"dir")
            mkdir(fullfile(opt.analysis,'newplots'))
        end

        % 2.1 Raster, trace, PSH, ISI whole trace
        if opt.plot_trial
           plot_extinction(allneurons{x,y}, allevents{x,y}, allspike{x,y}, allconditions{x,y}, opt, param)
        end

        % 2.2 Raster, trace, PSH, aligned
        if opt.plot_align
            plot_extinction_aligned(allneurons{x,y}, allevents{x,y}, allspike{x,y}, allconditions{x,y}, opt, param)
        end

        % 2.2 Firing rate (or Norm Firing rate) color plots
        if opt.plot_fireRate
            param.save = input.analysis;
            plot_fireRateColMap(allfireRate{x,y}, allconditions{x,y}, param, opt)
        end

    end
end