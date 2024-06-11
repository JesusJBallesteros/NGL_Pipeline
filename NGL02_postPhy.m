%% NGL02_postPhy (in progress)
% To run after manual curation of desired sessions is completed. Will read
% the resulting KS results after manual curation.
%
% Jesus 05.06.2024

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
            mkdir(opt.spikeSorted)
            mkdir(opt.trialSorted)

            %% 03. Extract preprocessed spikes and recover event data
            % Spike clusters after sorting and curation
            spike = loadSpikes(opt); % Also saves the result to \spikesorted
            if isfield(spike,"spike")
                spike = spike.spike; % Simplify loaded structure if needed
            end

            % Recover trial definitions created after event extraction and processing. 
            % Will have as many variations as requested at that time. Needs to be ran 
            % again to create new alignments.
            if ~exist('trialdef','var'),  load(fullfile(opt.trialSorted, "trialdef.mat")),    end
            
            %% 04. Iterate trough all units and sorts them into their trials
            % TODO fix Fieldtrip extraction
            [neurons, neurons_FT] = sort2trials(spike, trialdef, opt);
        
            %% ...

    end
end