function [events, trialdef, EventRecord, opt] = EventProcess(input, opt)
% EventProcess  Extract and structure behavioural events from INTAN or Deuteron.
%
% PURPOSE:
%   Central extraction of events. Checks for already cached EventRecord,
%   trialdef, and events files; only re-extracts what is missing or newly
%   requested. After extraction, calls trialdefGen to build trial boundaries
%   and runs conditions_script.m to group trials into experimental conditions.
%   Results are saved to disk.
%
% USAGE:
%   [events, trialdef, EventRecord, opt] = EventProcess(input, opt)
%
% INPUTS:
%   input  - struct with input.run, input.sessions(x).info.fileformat,
%              input.exefile, and path fields
%   opt    - complete options struct; relevant:
%              .RetrieveEvents  (logical) whether to extract events at all
%              .alignto         (cell of char) event names for trial zero
%              .trEvents        (cell of char) ITI-period special events
%              .FolderProcDataMat, .trialSorted  output paths
%
% OUTPUTS:
%   events      - struct; one field per alignment event, each containing
%                   trial-aligned timestamps in ms (NGL convention)
%   trialdef    - (2 × n_alignments) cell array:
%                   row 1: event names (char)
%                   row 2: (nTrials × 3) arrays [start end t0] in ms
%   EventRecord - struct; raw event list with timestamps
%   conditions  - NOT OUTPUT, but SAVED to disk directly.
%   opt         - updated with opt.eventdef and opt.newEvent
%
% CACHING:
%   If EventRecord.mat, trialdef.mat, and events.mat all exist and all
%   requested alignto events are present in events, returns immediately
%   without re-extraction or re-computation. Only re-runs the parts that
%   are missing or newly requested.
%
% SAVED FILES to opt.FolderProcDataMat and opt.trialSorted:
%   EventRecord.mat, trialdef.mat, events.mat, condition.mat
%
% REQUIRES:
%   eventDefinitions.m and conditions_script.m in analysisCode\
%
% CALLS:
%   INTAN_ExtractEvents, Deuteron_ExtractEvents, trialdefGen
%
% Jesus 07.05.2026

%% Prepare.
% Read eventcode list
if ~isfield(opt,'eventdef')       
    opt.eventdef  = eventDefinitions(input.sessions(input.run(1)).info.fileformat);
end 

% Empty petition for new events
opt.newEvent = {}; % see if check==3

%% Create empty outputs
events      = []; % If remains empty, data shall be treated as continuous.
trialdef    = [];
EventRecord = [];
conditions  = [];
condition   = [];

%% Check for alredy collected events
check = 0;
if isfile(fullfile(opt.FolderProcDataMat, strcat('EventRecord.mat')))
    disp('EventRecord found. Loading.')
    load(fullfile(opt.FolderProcDataMat, strcat('EventRecord.mat')), 'EventRecord');
    check = check+1;
    
    if isfile(fullfile(opt.trialSorted, strcat('trialdef.mat')))
        disp('Trial definition file found. Loading')
        load(fullfile(opt.trialSorted, strcat('trialdef.mat')), 'trialdef');
        check = check+1;
    end

    if isfile(fullfile(opt.trialSorted, strcat('events.mat')))
        disp('Events have been collected. Loading')
        load(fullfile(opt.trialSorted, strcat('events.mat')), 'events');
        check = check+1;
    end

    if check == 3
        % New event aligment requested?
        if ~all(ismember(opt.alignto, fieldnames(events)))
            opt.newEvent = opt.alignto(~ismember(opt.alignto', fieldnames(events)));

            for n = 1:length(opt.newEvent)
                if regexp(opt.newEvent{n}, 'bhv')
                    bhvtype = '(\w+)(\d+)';
                    bhvtype = regexp(opt.newEvent{n},bhvtype,'tokens');
                    opt.newEvent{n,1} = bhvtype{1,1}{1,1};
                    opt.newEvent{n,2} = bhvtype{1,1}{1,2};

                    if ismember(opt.newEvent{n,1},fieldnames(events))
                        opt.newEvent{n,1} = [];
                        opt.newEvent{n,2} = [];
                    end
               end
            end
        else
            disp('Events exist and no new events were requested. Loaded.')
            return
        end
    end % Done here 
end

%% Retrieve Events and generate needed variables
if opt.RetrieveEvents
    if check < 1
        % Proceed to extract all events captured by DEUT/INTAN acquisition system,
        % stored along with the data and synchronized with it (proper timestamped).
        switch input.sessions(input.run(1)).info.fileformat
            case {'DT2', 'DF1'}
                disp('Retrieving events from Deuteron Event files using EXE.')
                [EventRecord, opt] = Deuteron_ExtractEvents(input, opt);
            
            case {'fileperch', 'filepertype'}
                % INTAN
                disp('Retrieving events from INTAN Dig-IN channels.')
                EventRecord = INTAN_ExtractEvents(input, opt);

            otherwise
                % It is FT, keep going.
        end
    end

    if check < 3
        % Then, based on the trial definitions (defaulted or given) create
        % an 'events' struct fitting the NGL convention
        disp('Creating trial definitions based on extracted EventRecord and Eventcodes descriptions.')
        [events, trialdef, opt.eventdef, EventRecord] = trialdefGen(EventRecord, opt);

    elseif check==3 && ~isempty(opt.newEvent)
        % Re create an additional field in events as requested
        disp('Re-Creating Events variable with additional alignments.')
        [events, trialdef, ~, EventRecord] = trialdefGen(EventRecord, opt, events, trialdef);
    end

%% Run the personalized script for the conditions to be extracted
run('conditions_script.m');
if isempty(conditions) || ~isempty(condition)
    conditions=condition; clear condition
end

%% Save this session events, trialdef and conditions variables.
save(fullfile(opt.FolderProcDataMat, strcat('EventRecord.mat')), 'EventRecord', '-v7.3');
save(fullfile(opt.trialSorted, strcat('trialdef.mat')), 'trialdef', '-v7.3');
save(fullfile(opt.trialSorted, strcat('events.mat')), 'events', '-v7.3');
save(fullfile(opt.trialSorted, strcat('condition.mat')), 'conditions', '-v7.3');

else
    disp('Events not requested. Skipped.')
    return
end

end