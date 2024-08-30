function [events, trialdef, EventRecord, conditions] = EventProcess(opt)
% Function meant to put together all possible ways to extract events from
% Deuteron and INTAN systems.
%
% Jesus 23.08.2024

%% Defaults.
if ~isfield(opt,'useexe'),          opt.useexe              = true;                 end
if ~isfield(opt,'ext'),             opt.ext                 = 'fileperch';          end
if ~isfield(opt,'eventdef'),        opt.eventdef            = eventDefinitions(opt.ext);   end
if ~isfield(opt,'trEvents'),        opt.trEvents            = [];                   end

opt.exefile = 'C:\Code\ephys-data-pipeline\toolboxes\Deuteron\software\Event_File_Reader_9_0.exe';

%% Prep. Make sure we create empty outputs
events      = []; % If remains empty, data shall be treated as continuous.
trialdef    = [];
EventRecord = [];

%% Check for alredy collected events
if isfile(fullfile(opt.FolderProcDataMat, strcat('EventRecord.mat')))
    disp('EventRecord found. Loading.')
    load(fullfile(opt.FolderProcDataMat, strcat('EventRecord.mat')), 'EventRecord');

    if isfile(fullfile(opt.trialSorted, strcat('trialdef.mat')))
        disp('Trial definition file found. Assuming events have been collected. Loading both.')
        load(fullfile(opt.analysis, strcat('events.mat')), 'events');
        load(fullfile(opt.trialSorted, strcat('trialdef.mat')), 'trialdef');
        load(fullfile(opt.FolderProcDataMat, strcat('EventRecord.mat')), 'EventRecord');
        load(fullfile(opt.analysis, strcat('condition.mat')), 'conditions');

        % Done here
        return
    end
end

%% Retrieve Events and generate needed variables
if opt.RetrieveEvents
    switch opt.ext
        case {'DT2', 'DF1'} 
            % DEFAULT. When extracting events from SD using EXE.
                % Proceed to extract all events captured by Deuteron acquisition system,
                % stored along with the data and synchronized with it (proper timestamped).
            disp('Retrieving events from Deuteron Event files using EXE.')
            [EventRecord, opt] = Deuteron_ExtractEvents(opt);
        
        case {'fileperch', 'filepertype'}
            % INTAN
            disp('Retrieving events from INTAN Dig-IN channels.')
            EventRecord = INTAN_ExtractEvents(opt);           
    end

    %% Then, based on the trial definitions (defaulted or given) create
    % an 'events' struct fitting the NGL convention
    disp('Creating trial definitions based on extracted EventRecord and Eventcodes descriptions.')
    [events, trialdef, opt.eventdef] = trialdefGen(EventRecord, opt);

    %% Create the conditions variable for trial indexing
    conditions = get_trialConditions(opt, events);

else
    disp('Events not requested. Skipped.')
    return
end

%% Save this session events, trialdef and conditions variables.
save(fullfile(opt.FolderProcDataMat, strcat('EventRecord.mat')), 'EventRecord', '-v7.3');
save(fullfile(opt.trialSorted, strcat('trialdef.mat')), 'trialdef', '-v7.3');
save(fullfile(opt.analysis, strcat('events.mat')), 'events', '-v7.3');
save(fullfile(opt.analysis, strcat('condition.mat')), 'conditions', '-v7.3');

end