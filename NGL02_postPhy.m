%% NGL02_postPhy (in progress)
% To run after manual curation of desired sessions is completed. Will read
% the resulting KS results after manual curation.
%
% Jesus 12.06.2024

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

%% 03. Proceed with data per session
for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
            input.run = [x y]; % Current run, to pass to functions.
            
            %% 04. Prepare to proceed with a single session.
            [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);           

            %% 03. Extract preprocessed spikes and recover event data.
            % Spike clusters after sorting and curation.
            spike = loadSpikes(opt); % Also saves the result to \spikesorted
            if isfield(spike,"spike"), spike = spike.spike; end % Simplify loaded structure if needed
            
            %% 04. Calculate the ISI histogram for all clusters
            [spike.isihist] = calc_isihist(spike);
            
            % Recover trial definitions created after event extraction and processing. 
            % Will have as many variations as requested at that time. Needs to be ran 
            % again to create new alignments.
            if ~exist('trialdef','var'), load(fullfile(opt.trialSorted, "trialdef.mat")); end
            
            %% 04. Iterate trough all units and sorts them into their trials.
            % Outputs are saved to data\analysis.
            % TODO fix Fieldtrip extraction
            [neurons, ~] = sort2trials(spike, trialdef, opt);
        
            %% ...

    end

%% 04. Proceed with data as whole
% TODO    

end