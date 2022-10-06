%% This Script is intended to serve as index of .NWB functionalities.
% TO BE redacted
% Most likely, funtions here should live as independent, to be called in
% any given pipeline. That, assuming they are not easier to work with in
% their native form.
% 
% MOST IMPORTANT, first git-clone or download the 'nwbmat' toolbox.
% It is in my Code folder, in server TMP.
% https://github.com/NeurodataWithoutBorders/matnwb
% 
% Then, it is worth visit: https://neurodatawithoutborders.github.io/matnwb/ 
%
% % I found useful the use of 'HDFView' software. With it, one can open a .nwb
% %  file and navigate it as a folder system. I don't think it allows to
% %  edit or add new info, tho.
%
% Last modified Jesus 06.10.2022

%% 00. Needed input
% Although this script makes sense more as a compilation of commands and 
% routines, I guess it could become a Script for a bulk of sessions from a
% single animal. Not sure. Anyways... locating things never hurts
input.mainfolder = pwd;                 % string. Gets current folder by default
input.datafolder = "D:\Experiments\";   % char.
input.animal     = 'FAT';               % string. TES, FAT, FRN
input.dates      = {'20220909'};        % cell array. Probably only one session will be feeded at a time

%% 00. Dependencies
cd(input.mainfolder) % Code folder
cd functions\toolboxes\matnwb; % Get in toolbox folder
addpath(genpath(pwd)); % recursively add all subfolders
% generateCore()  % Install the API. Only first time, I think. Needs to be in 'nwbmat' folder

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

% General ID info:
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

% Adding BEHAVIOR data:
% This could be useful for arena position, eye movements, peaking ...
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

% Adding Trials information:
    % Here is a long call, which creates a 'trials' object with many fields 
    trials = types.core.TimeIntervals( ...  % Is a 'TimeInterval' table
            'colnames', {'start_time', 'stop_time', 'correct'}, ... % Headers
            'description', 'trial data and properties', ...         % module description
            'id', types.hdmf_common.ElementIdentifiers('data', 0:2), ... % three trials (0, 1, 2)
            'start_time', types.hdmf_common.VectorData( ...         % with vectors:
                'data', [0.1, 1.5, 2.5], ...                        % for start_time
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

% TODO Adding Processed VOLTAGE data:
% This could be useful to add the processed LFP, but potentially to add any 
%   spectrographic, coherence, etc. :
% TODO

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

