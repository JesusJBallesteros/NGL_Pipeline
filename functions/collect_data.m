function [allneurons, allevents, allconditions, allspike, allblobs] = collect_data(input, opt)
%
%
%
%

if exist(fullfile(input.analysis, "data_all.mat"),"file")
    % If data has been already collected, eazy
    load(fullfile(input.analysis, "data_all.mat"), ...
        'allneurons', 'allevents', 'allconditions', 'allspike', 'allblobs');
else
    % if not, try to collect data from in a subject/session basis
    for x = 1:input.nsubjects % Go over subjects
        for y = 1:input.sessions(x).nsessions % Go over sessions
            input.run = [x y];
            [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);

            if ~exist('events','var'),    load(fullfile(opt.analysis, "events.mat")),              end
            if ~exist('neurons','var'),   load(fullfile(opt.analysis, "neurons.mat")),             end
            if ~exist('condition','var'), load(fullfile(opt.analysis, "condition.mat")),           end
            if ~exist('spike','var'),     load(fullfile(opt.spikeSorted, "spike.mat")),            end
            if ~exist('blob','var') && opt.plot_SocLear, load(fullfile(opt.analysis, "blob.mat")), end

            allneurons{x,y}     = neurons;      % Collect Subjects/session data
            allevents{x,y}      = events;       % idem
            allconditions{x,y}  = conditions;   % idem
            allspike{x,y}       = spike;        % idem
            allblobs{x,y}       = [];
            if opt.plot_SocLear, allblobs{x,y} = blob; end % Include blob tracking if Social

            clear neurons events conditions spike blob % Clean up Subjects/Session results after collection
        end % Sessions loop
    end % Subjects loop
    
    % Save for posterior uses
    save(fullfile(input.analysis, "data_all.mat"),'allspike','allconditions','allevents','allneurons');

    % Add blob data if required
    if exist('allblobs','var'), save(fullfile(input.analysis, "data_all.mat"), 'allblobs','-append');end  
end
end