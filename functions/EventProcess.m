function [events, trialdef, EventRecord] = EventProcess(opt)
% Function meant to put together all possible ways to extract events from
% Deuteron and INTAN systems.
%
% Jesus 17.10.2024

%% Defaults.
if ~isfield(opt,'useexe'),          opt.useexe              = true;                end
if ~isfield(opt,'ext'),             opt.ext                 = 'fileperch';          end
if ~isfield(opt,'eventdef'),        opt.eventdef            = eventDefinitions(opt.ext); end
if ~isfield(opt,'trEvents'),        opt.trEvents            = [];                   end
if ~isfield(opt,'addtime'),         opt.addtime             = 0;                    end

opt.exefile = 'C:\Code\ephys-data-pipeline\toolboxes\Deuteron\software\Event_File_Reader_9_0.exe';

%% Create empty outputs
events      = []; % If remains empty, data shall be treated as continuous.
trialdef    = [];
EventRecord = [];
conditions  = [];

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

    if isfile(fullfile(opt.analysis, strcat('events.mat')))
        disp('Events have been collected. Loading')
        load(fullfile(opt.analysis, strcat('events.mat')), 'events');
        check = check+1;
    end

    if check == 3, return; end % Done here 
end

%% Retrieve Events and generate needed variables
if opt.RetrieveEvents
    if check < 1
        % Proceed to extract all events captured by DEUT/INTAN acquisition system,
        % stored along with the data and synchronized with it (proper timestamped).
        switch opt.ext
            case {'DT2', 'DF1'}
                disp('Retrieving events from Deuteron Event files using EXE.')
                [EventRecord, opt] = Deuteron_ExtractEvents(opt);
            
            case {'fileperch', 'filepertype'}
                % INTAN
                disp('Retrieving events from INTAN Dig-IN channels.')
                EventRecord = INTAN_ExtractEvents(opt);
            otherwise
                % It is FT, keep going.
        end
    end

    % Then, based on the trial definitions (defaulted or given) create
    % an 'events' struct fitting the NGL convention
    disp('Creating trial definitions based on extracted EventRecord and Eventcodes descriptions.')
    [events, trialdef, opt.eventdef] = trialdefGen(EventRecord, opt);

    %% Run the personalized script for the conditions to be extracted
    run('conditions_script.m');

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