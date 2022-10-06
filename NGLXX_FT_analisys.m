%% Jesus' Pipeline to analyze FieldTrip formatted files
% TO BE redacted
% TO BE DONE
% One of the first functions here will split the continous dataset into
% trials, according to a pre-fixed length or a trial scheme given by the
% user.
% 
% Last modified Jesus 30.09.2022

% Dependencies
cd(input.mainfolder)
addpath functions\
addpath functions\toolboxes\fieldtrip_light

% Set defaults
input.mainfolder = 'C:\Code\Scripts\NGL_ephys_data_pipeline';
input.datafolder = "D:\Experiments";
ft_defaults

%% 00. Needed input
input.animal   = 'TES';

% Chunk continous data into trials? For now, fixed length. 
% TODO: Give Trial info as Fieldtrip cfg to directly create it.
input.parsing.cfg     = [];
 input.parsing.cfg.length  = 5; % in sec.

%% 01. Find sessions
% Get FT files in animal folder
input.folder = "D:\Experiments\" + input.animal;
cd(input.folder)
input.files = dir('*_FT.mat');
