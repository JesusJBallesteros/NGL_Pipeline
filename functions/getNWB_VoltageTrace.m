function volt_data = getNWB_VoltageTrace(input,chIndex,tIndex,dwnsmpl)
cd(input.mainfolder) % Code folder
cd toolboxes\matnwb; % Get in toolbox folder
addpath(genpath(pwd)); % recursively add all subfolders

if iscell(input.dates)
    input.dates = input.datafolder + input.animal + '\' + input.animal + '_' + input.dates(:);
    [sessions.folder,sessions.name,~] = fileparts(input.dates);
    sessions.folder = unique(sessions.folder);
    sessions.nSessions = length(sessions.name);

    % To the session folder
    cd(sessions.folder + '\' + sessions.name)
else
    try
        % To the session folder
        cd(input.dates)
    catch
        disp('Input dates are not valid.')
        return
    end
end

% Find the .nwb file.
sessions.nwbfile = dir('*.nwb'); 

% Read file schema. Not data.
nwb = nwbRead(sessions.nwbfile.name);

% read_general = nwb.general_devices.get('Intan USB Interface Board'); 

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
%     chIndex   = 1:5:32; % ch 1-5
if strcmp(tIndex,'all')
    tstamps = size(read_ephys_series.timestamps);
    time = [1:dwnsmpl:tstamps];
else    
    time = 1:dwnsmpl:(3000*60)+1; % timestamps for 1st minute
end

% to call:
volt_data = read_ephys_series.data(chIndex, time);
% Which already gives a nCh*samples array of data.

end