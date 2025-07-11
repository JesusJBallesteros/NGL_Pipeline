%% Add a few details about your probe. 
% Add new one if needed, keep old ones as archive.
clear all

% E32+R-50-S2-L10-200 NT
% mapName = 'chanMapE32-S2.mat';
% poly    = 1;  % integer. Linear (1) / politrode (>1) configuration.
% nChan   = 32; % integer, number of contacts.
% nShank  = 2;  % integer, number of shanks.
% interCh = [200, 50]; % um between contacts [x, y])

% % E32Tri
mapName = 'chanMapE32-Tri.mat';
poly    = 3;    % integer. Linear (1) / politrode (>1) configuration.
nChan   = 32;   % integer, number of contacts.
nShank  = 1;    % integer, number of shanks.
interCh = [50, 50]; % integer, um between contacts [x, y])
rows    = 12;   % integer, specify total rows (distribution tends to not be homogeneous
offset  = 25;    % integer, um offset between central and lateral colums in politrodes

% Sampling frequency here, because of reasons.
fs = 32000; 

%% Generate map. Few details might need customization, mostly for non-linear configurations.
% List of channel indices (and give an index to dead channels too). 
% chanMap(1) is the row in the raw binary file for the first channel. 
chanMap = 1:nChan;

% Declare which channels are "connected" in this normal ordering, meaning not 
% dead or non-ephys data
connected = true(nChan, 1);

% Now we define the horizontal (x) and vertical (y) coordinates of these
% channels. In um here but the scaling doesn't really matter. 
if poly == 1
    xcoords = repmat(interCh(1):interCh(1):interCh(1)*nShank, [1, nChan/nShank])';
else
    % customized for non-existing channels in Atlas E32-Tri
    xcoords = repmat(interCh(2):interCh(2):interCh(2)*poly, [1, rows])';  
    xcoords([31,33,34,36],1) = 0; % removes non-existing contacts in these lateral columns
    xcoords(xcoords==0) = [];
end

if poly == 1
    ycoords = sort(repmat(interCh(2)*(nChan/nShank):-interCh(2):interCh(2), [1, nShank])'); % sort is helpful
else
    % customized for non-existing channels in Atlas E32-Tri
    ycoords = sort(repmat(interCh(2)*rows:-interCh(2):interCh(2), [1, poly])'); % sort is helpful

    % Offset central column
    ycoords(2:poly:poly*(round(nChan/poly)+1),1) = ycoords(2:poly:poly*(round(nChan/poly)+1),1) - offset;

    ycoords([31,33,34,36],1) = 0; % removes non-existing contacts in these lateral columns
    ycoords(ycoords==0) = [];
end

% Often, multi-shank probes can be organized into groups of channels that cannot 
% possibly share spikes with the rest of the probe. This helps the algorithm 
% discard noisy templates shared across groups. 
% All == 1, If all channels are on the same shank or there is not obvious grouping.
kcoords =repmat(1:1:nShank, [1, nChan/nShank])';

% Plot resulting map.
scatter(xcoords,ycoords,[],"black","filled");
 title(mapName)

%% Save.
save(mapName, 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs')