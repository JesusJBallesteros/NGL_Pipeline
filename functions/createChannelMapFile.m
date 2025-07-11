%% Creates a channel map file
% 'kcoords' is used to forcefully restrict templates to channels in the same
% channel group. Creates labeled clusters to better differentiate grouped
% tetrodes, or shanks.
%
% An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is:
%    ops.criterionNoiseChannels = 0.2;
%     if < 1, it will be treated as a fraction of the total number of clusters
%     if > 1, it will be treated as the "effective number" of channel groups 
%     at which to set the threshold. So if a template occupies more than this many 
%     channel groups, it will not be restricted to a single channel group. 
%
% Created from Kilosort original script, with additions from Sara's 'Howto
% createachannelmap' document. Generalized and simplified. TESTING.
%
% Jesus 10.05.2023 

%% Probe info 
chanMapName = 'chanMapDefault'; % Keep 'chanMap*' nomenclature.
Nchannels   = 32;
connected   = true(Nchannels, 1); % Keep all
% fs          = input.sessions(input.run(1)).info.amplifier_sample_rate; % sampling frequency (why?)
fs          = ops.fs;%input.sessions(input.run(1)).info.amplifier_sample_rate; % sampling frequency (why?) %above did not work

%% Probe physical layout.
% Write it as a vector. Write the HEADSTAGE Pins in the order that
% represents the probe layout, i.e a pure linear one:
chanMap     = 1:Nchannels;

% % Or (Masahiro's Neuronexus tetrode)
% chanMap     = [18 19 17 21  7 22  6 23 ...
%                 0  1  4  2 16  3 20  5 ...
%                14 15 13 11 12 31 10 27 ...
%                28 29 26 30 25  8 24  9];

% Zero index it for KS:
chanMap0ind = chanMap - 1; 

%% Simple linear
xcoords   = ones(Nchannels,1);
ycoords   = [1:Nchannels]';
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

% % Tetrode map example
% xcoords   = repmat([1 2 3 4]', 1, Nchannels/4);
%     xcoords   = xcoords(:);
% ycoords   = repmat(1:Nchannels/4, 4, 1);
%     ycoords   = ycoords(:);
% kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

%% Save map
% Change name 
save(fullfile(input.analysisCode, [chanMapName '.mat']), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs');
