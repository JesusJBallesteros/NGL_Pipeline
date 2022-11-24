%% This Script is intended to serve as index of .NWB functionalities.
% TO BE redacted
% Most likely, funtions here should live as independent, to be called in
% any given pipeline. That, assuming they are not easier to work with in
% their native form.
% 
% MOST IMPORTANT, first git-clone 'nwbmat' toolbox.
% https://github.com/NeurodataWithoutBorders/matnwb
% 
% Then, it is worth visit: https://neurodatawithoutborders.github.io/matnwb/ 
%
% % I found useful the use of 'HDFView' software. With it, one can open a .nwb
% %  file and navigate it as a folder system. I don't think it allows to
% %  edit or add new info, tho.
%
% Last modified Jesus 10.2022

%% 00. Needed input
% Although this script makes sense more as a compilation of commands and 
% routines, I guess it could become a Script for a bulk of sessions from a
% single animal. Not sure. Anyways... locating things never hurts
input.mainfolder = 'C:\Code\Scripts\ephys-data-pipeline'; % string. Gets current folder by default
input.datafolder = "D:\Experiments\";   % char.
input.animal     = 'FRN';               % string. TES, FAT, FRN
input.dates      = {'20181029'};        % cell array. Probably only one session will be feeded at a time

%% 00. Dependencies
cd(input.mainfolder) % Code folder
addpath functions\
cd toolboxes\matnwb; % Get in toolbox folder
addpath(genpath(pwd)); % recursively add all subfolders
% generateCore()  % Install the API. Only first time, I think. We need to be at '\nwbmat'

%% Locate and get session
% Locate and get in session folder
if iscell(input.dates)
    input.dates = input.datafolder + input.animal + '\' + input.animal + '_' + input.dates(:);
    [sessions.folder,sessions.name,~] = fileparts(input.dates);
    sessions.folder = unique(sessions.folder);
    sessions.nSessions = length(sessions.name);
end

% To the session folder
cd(sessions.folder + '\' + sessions.name)

% Find the .nwb file.
sessions.nwbfile = dir('*.nwb'); 

% Read file schema. Not data.
nwb = nwbRead(sessions.nwbfile.name);

%% INTERIM
% From here, we are ready to read, edit or create information at data and
%  metadata levels. For now we basically have the RAW data from INTAN as
%  extracted from 'amplifier.dat' files, if using the typical 'fileperch' format. 
% I think, (TO TEST) if we have used 'lowpass.dat' files, the LFP should 
%  exist as well?
% We could take the raw voltage (30KHz) and filter and downsample it to LFP,
%  then save the LFP in the .nwb file separatedly from the raw. We could 
%  add electrode information (shanks, contact location, etc), also trial
%  information, video tracking (x,y coordinates), spike and cluster information... . 

%% Adding info
% After we have created the .nwb file directly from INTAN format, we find
% out that the metadata is not very exhaustive, and some data may be wrong. 
%  (Because the default transformation, not modified to not alter it's 
%  behavior, is very basic).
% However, it is possible, and not very complicated, to add metadata and data.
% Here are few relevant examples from the online tutorial, described.

%% General ID info:
    subject = types.core.Subject( ... % Generate a 'Subject' OBJECT with properties
        'subject_id', input.animal, ... % ID
        'age', 'NaN', ...               % age
        'description', 'Pigeon01', ...  % etc.
        'species', 'Columba sp', ...    %
        'sex', 'M' ...                  %
        );
    nwb.general_subject = subject;  % Then, add the object to the specific nwb field
    clear subject                   % We can get rid of this object now.
    
    % Now calling: 
     nwb.general_subject % should print the metadata, 
    % and so could be indexed into a single variable.

% Adding behavior data.

    %% Adding RAW EVENT CODES and TIMESTAMPS.
    %  'BehavioralEvents' is an interface for storing behavioral events. 
    %  We can use it for storing the timing of stimuli, trial outputs, 
    %  or peacking times.
        
        % For a generic set of arrays of events and timestamps
        event_codes      = [1, 3, 6, 8]; % Mix, uninterpreted, event codes in a continuous fashion
        event_timestamps = [1.0, 1.5, 1.8, 2.5]; % corresponding timestamps
         
        % Create a TimeSeries with appropiate properties
        time_series = types.core.TimeSeries( ...
            'data', event_codes, ...
            'timestamps', event_timestamps, ...
            'description', 'Event codes sent by software into acquisition system', ...
            'data_unit', 'NaN' ...
        );
         
        % Create an object for behavioral events
        behavioral_events = types.core.BehavioralEvents();

        % And include the TimeSeries containing the event codes
        behavioral_events.timeseries.set('trial_output', time_series);
         
        % behavior_processing_module = types.core.ProcessingModule("stores behavioral data.");  % if you have not already created it
        
        % In a processing module created for behavior, set the behavioral events
        behavior_processing_module.nwbdatainterface.set('BehavioralEvents', behavioral_events);
        
        % And add the behavior processing module to the main nwb data. (if not done yet)
        nwb.processing.set('behavior', behavior_processing_module); % if you have not already added it

        % get rid of objects
        clear event_codes event_timestamps time_series behavioral_events

    %% Adding SPATIAL coordinates.
    %  This could be useful for arena position, eye movements, peaking ...
    %  In this example I call it here 'spatialdata':
        spatialdata = types.core.SpatialSeries( ...          % create a SpatialSeries OBJECT
            'data', [linspace(0,10,100); linspace(0,8,100)], ... % [X1,X2,... ; Y1,Y2,...] positions
            'reference_frame', '(0,0) is bottom left corner', ... % Framework description
            'timestamps', linspace(0, 100)/200 ...              % time for each position
            );
        
        % Create a 'Top_view' OBJECT adding the 'spatialdata' to it. This would
        % assume it is 'possition in a plane' (i.e. arena). In this example 
        % imagine it is data of position recorded from the top of the arena.
        Top_view = types.core.Position('Top_view', spatialdata);
    
        % In parallel, create a 'processing module'. This is a "folder" apart from 
        %  'acquisition' that stores somehow processed data, to separate it from
        %  pure raw data. (Note that extracted position from a video is somehow
        %  'processed', versus the video frames themselves, whose would be the
        %  raw data). Let's name it 'behavior_module':
        behavior_module = types.core.ProcessingModule('description', 'contains behavioral data');
    
        % Insert the Top_view object into the module, as 'Position'. We annidate this because 
        %  we could have multiple positions from different angles, for other 
        %  porpouses (e.g. a screen-level camera). 
        % Or we could have other behavioral measurements (not positions).
        % To add ('set') it as 'Position' data, tagged as Top_view:
        behavior_module.nwbdatainterface.set('Position', Top_view);
    
        % Finally, we bring the behaviour module, with the position data on it,
        % to the .nwb file. We set it as behavior data under Processing.
        nwb.processing.set('behavior', behavior_module);
        clear Top_view behavior_module % We can get rid of the intermediates
    
        % Now we could obtain the data back from the .nwb file calling:
        read_bhv_series = nwb.processing.get('behavior'). ...     % Inside the behavior module
                            nwbdatainterface.get('Position'). ... % We have Position data
                            spatialseries.get('Top_view');        % coming from the Top_view
        % Which gives us an Object with all metadata plus '.data' and '.timestamps' 
        %  fields (data itself). Index them , for example, as:
        position_data = read_bhv_series.data; % for the x,y coordinates
        % Note that this data can be indexed to get only portions of it:
        % position_data = read_bhv_series.data(1:2, 1:10); both x,y for only 10 frames
    
        % To get used to this structures, I understand it as asking to get
        %  the 'behavior' that we have stored as processed data.
        % From it, we want data tagged as 'Position', and no ther (which could live here).
        % From it, we want the Positions tagged 'Top_view'. (Which are 'SpatialSeries')

    %% Adding TRIAL information.
    % This is a variation to adding EVENTS, In this case we would have already 
    %  extracted and interpreted the EVENTCODES to its actual meaning in a 
    %  trial context different approach. Having timestamp arrays for each event code,
    %  they can be added as the timing for each category.
    % Here is a long call, which creates a 'trials' object with desired fields 
        trials = types.core.TimeIntervals( ...  % Is a 'TimeInterval' table
                'colnames', {'start_time', 'stop_time', 'correct'}, ... % Headers
                'description', 'trial data and properties', ...         % module description
                'id', types.hdmf_common.ElementIdentifiers('data', 0:2), ... % three trials (0, 1, 2)
                'start_time', types.hdmf_common.VectorData( ...         % with vectors:
                    'data', [0.1, 1.5, 2.5], ...                        % timestamps for start_time
   	                'description','start time of trial in seconds' ... 
                    ), ...
                'stop_time', types.hdmf_common.VectorData( ...          
                    'data', [1.0, 2.0, 3.0], ...                        % for end_time
   	                'description','end of each trial in seconds' ...
                    ), ...
                'correct', types.hdmf_common.VectorData( ...
                    'data', [false, true, false], ...                   % and for output
   	                'description', 'whether the trial was correct') ...
                );
    
        % Add it to the .nwb file in the corresponding folder 'intervals_trials'
        nwb.intervals_trials = trials;
        clear trials % And get rid of the variable

%% TODO Adding Processed VOLTAGE data
% This could be useful to add the processed LFP, but potentially to add any 
%   spectrographic, coherence, etc. :

%% Searching for general acquisition info
read_general = nwb.general_devices.get('Intan USB Interface Board'); 

%% Searching and reading VOLTAGE data
% We can check if the ephys data is readable. The voltage data are stored
% in ElectricalSeries (a subclass of TimeSeries). These data are referenced
% to a set of rows in the electrodes table. It is data recorded directly by
% those electrodes, so it gets into '.nwb.acquisition' as RAW data.

% We can read the data schema with:
read_ephys_series = nwb.acquisition.get('ElectricalSeries'); 

% It gets the electrical series from acquisition. It is still an schema and,
%  if we want to get the actual voltage data, we can do it all at once or
%  index a smaller portion by:
    chIndex   = 1:5; % ch 1-5
    timeIndex = 10001:40000; % timestamps for 1 second

    % to call:
    data_chunk = read_ephys_series.data(:, timeIndex);
    % Which already gives a nCh*samples array of data.


%% Read events from INTAN's DIG IN to bring into the NWB file
[ev2, ~, ~] = readEvents('path',pwd);
%TODO

%% Searching and reading EVENTS data (DIG IN)
% TODO
% It uses the main code from 'readEvents', applied to data extracted from
% the .nwb files and not from INTAN files. It yields the same results.
% To read this data, we parallelize the reading of relatively big chunks of
% time from all possible DIG IN channels. Then we concatenate the pieces
% and proceed to index events based on single bit changes.
read_events_series = nwb.acquisition.get('TimeSeries_digital_input'); 

max_nBits = read_events_series.data.internal.maxSize(1);
numSmp    = read_events_series.data.internal.maxSize(2);

stp = 1000000;
t_chunks = 0:stp:numSmp-mod(numSmp,stp);
t_chunks = [t_chunks numSmp];

parfor tt = 1:length(t_chunks)-1
    % Read 4 first pins for timechunks
    % Flips them (first pin least significant)
    % Converts vector to string
    % Converts string to decimal as uint16
    bin{tt} = uint16(read_events_series.data(1:max_nBits,t_chunks(tt)+1:t_chunks(tt+1)))';
end
clear t_chunks stp max_nBits numSmp

dIN = cell2mat(bin');
clear bin

smpRate = 30000;
skipDur = .5;       % look 0.5 ms after the first change to assume all pins are there.    
smpSkip = smpRate*skipDur/1000;
smpSum  = 2;        % an event code is sums over n samples (to catch instabilities)

% Find those samples where any pin is 1
idx = double(find(any(dIN,2)));

% the start of event codes (1ms after the first pin up)
ev.smpInd = idx(diff([-2; idx])>2)+smpSkip;

% Get the events and their corresponding timestamps
ev.eventStandard = zeros(length(ev.smpInd),1);
for i=1:length(ev.smpInd)
    % convert to decimal, sum over n samples to catch inconsitencies
    ev.eventStandard(i) = binvec2dec(sum(dIN(ev.smpInd(i):ev.smpInd(i)+smpSum,1:4)));
end

ev.eventExtra = zeros(length(ev.smpInd),1);
dIN(:,1:4) = 0; % set first 4 pins to 0
for i=1:length(ev.smpInd)
    % convert to decimal, sum over n samples to catch inconsitencies
    ev.eventExtra(i) = binvec2dec(sum(dIN(ev.smpInd(i):ev.smpInd(i)+smpSum,1:end)));
end
clear dIN

% remove the delay caused by skipping samples
ev.smpInd = ev.smpInd-smpSkip;
clear smpRate skipDur smpSkip smpSum

% Find standard events 0, and substitute by associated extra event
idx = ev.eventStandard==0;
events.rawEvents = [ev.smpInd ev.eventStandard];
events.rawEvents(idx,2) = ev.eventExtra(idx);
clear idx

% Remove repetitions (consecutive equal eventcodes not expected)
idx = find(diff(events.rawEvents(:,2))==0)+1;
events.rawEvents(idx,:) = [];

% Find each start of trial and grab all events until the next one.
% If the resulting vector contains more than 3 events, consider it valid
% and include it in events.trialEvents. If not, discard it.
idx = find(events.rawEvents(:,2)==1);
c = 0;
trial = [];
for i = 1:length(idx)-1
    tr = events.rawEvents(idx(i):idx(i+1)-1,2);
    if length(tr)>3
        c = c+1;
        trial{c} = tr';
    end
end
events.trialEvents = trial';
clear tr trial idx c i