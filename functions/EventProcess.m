function [events, trialdef, EventRecord] = EventProcess(input, opt)
% Function meant to put together all possible ways to extract events from
% Deuteron and INTAN systems.
%
% Jesus 17.10.2024

%% Defaults.
if ~isfield(opt,'useexe'),          opt.useexe              = true;                 end
if ~isfield(opt,'ext'),             opt.ext                 = 'fileperch';          end
if ~isfield(opt,'eventdef'),        opt.eventdef            = eventDefinitions(opt.ext);   end % The script 'eventDefinitions.mat' must be inside your project file system, under 'analisysCode', and a template exists in the folder 'configfiles' of the toolbox
if ~isfield(opt,'trEvents'),        opt.trEvents            = [];                   end
if ~isfield(opt,'addtime'),         opt.addtime             = 0;                    end
opt.newEvent    = {};
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
        switch opt.ext
            case {'DT2', 'DF1'}
                disp('Retrieving events from Deuteron Event files using EXE.')
                [EventRecord, opt] = Deuteron_ExtractEvents(opt);
            
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