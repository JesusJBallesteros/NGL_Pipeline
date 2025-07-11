function collect_allData(input, opt, type)
% Will loop throuhg all sessions and folders to get all data created until
% now. It will create project-wide variables and saved under \analysisCode
% folder. 
% INPUT: general pipeline 'input' and 'opt' variables, plus an explicit
%        string to determine the sort of data tha we need to collect and
%        save: 'spike', 'TFRtrial' or 'TFRcont'. Possibility to expand.
%
% OUTPUT: the data will be collected and saved as a mat file inside the 'analysis' folder 
%
% Jesus 22.04.2025

%% Determine data type and variables to collect
if strcmpi(type, 'spike')
    % Spike data collects the trialparsed spike variable and all other
    % relevant info
    filename = 'all_SpikeData.mat';
    varlist = {"allspike", "allconditions", "allevents", "allneurons"};
else
    % TFR data collects the Fieldtrip generated TFR data, ineither trial-parsed 
    % or continuous variants, and all other relevant info.
    filename = ['all_Data_', type, '.mat'];
    varlist = {"allspike", "allconditions", "allevents", "allneurons"};

    if contains(type,'cont'),       varlist{end+1} = "allTFR_continuous"; % DEVELOPMENT
    elseif contains(type,'trial'),  varlist{end+1} = "allTFR_trialparsed";
    end
end
if opt.plot_SocLear, varlist{end+1} = "allblobs"; end

%% Proceed
for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        % Prepare to proceed with a single session.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);           
        
        % Recover all files created on a session-by-session basis, to
        % agrupate on a single variable that can start to be treated at a
        % project level.
        allevents{x,y}      = load(fullfile(opt.analysis, "events.mat"));
        allneurons{x,y}     = load(fullfile(opt.analysis, "neurons.mat"));
        allconditions{x,y}  = load(fullfile(opt.analysis, "condition.mat"));
        if opt.plot_SocLear
            allblobs{x,y}   = load(fullfile(opt.analysis, "blob.mat"));
        end

        if strcmpi(type, 'spike')
            allspike{x,y} = load(fullfile(opt.spikeSorted, "spike.mat"));
        else
            foundfile = ls(fullfile(opt.spikeSorted,'*TFR*'));

            if contains(type,'cont')
                % DEVELOPMENT
                allTFR_continuous{x,y} = load(fullfile(opt.spikeSorted, foundfile));
            elseif contains(type,'trial')
                allTFR_trialparsed{x,y} = load(fullfile(opt.spikeSorted, foundfile));
            end
        end
    end
end

save(fullfile(input.analysis, filename), varlist{:})
disp('Spike data for the whole project colleted and saved as ''all_SpikeData.mat''.')

end