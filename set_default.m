function set_default(input)
% 'set_default' adds the dependencies, included under the main folder.
%
% It reads the inputs, if any, and validates them.
%
% If not valid values are provided, it ask for the necessary ones and puts 
% them in correct format. If any necessary input continues to be invalid or
% empty, throws error.
% 
% For missing optionals, it uses defaults.
%
% Jesus. 24.01.2023

%% Check necessary inputs
% Toolbox Main folder
if ~isfield(input,'mainfolder') || isempty(input.mainfolder)
    input.mainfolder = inputdlg('Toolbox absolute folder:',...
                        'Toolbox foldder', [1 50], {'C:\Code\Scripts\ephys-data-pipeline'});
    
    if isempty(input.mainfolder)
        error('input.mainfolder is empty. Please provide a valid one and try again.')
    else
        input.mainfolder = char(input.mainfolder);
    end

elseif ~ischar(input.mainfolder)
    input.mainfolder = char(input.mainfolder);
end
    
% Data Main folder
if ~isfield(input,'datafolder') || isempty(input.datafolder)
    input.datafolder = inputdlg('Data absolute folder:',...
                        'Data folder', [1 50], {'D:\Experiments\'});
        
    if isempty(input.datafolder)
        error('input.datafolder is empty. Please provide a valid one and try again.')
    else
        input.datafolder = string(input.datafolder);
    end

elseif ~isstring(input.datafolder)
    input.datafolder = string(input.datafolder);
end
    
% Animal code
if ~isfield(input,'animal') || isempty(input.animal)
    input.animal = inputdlg('Code:',...
                    'Animal code', [1 50], {'DOE'});
        
    if isempty(input.datafolder)
        error('input.animal is empty. Please provide a valid one and try again.')
    else
      input.animal = char(input.animal);
    end

elseif ~ischar(input.animal)
     input.animal = char(input.animal);
end
    
%% Optional Inputs
% Dates
if ~isfield(input,'dates') || isempty(input.dates) || ~iscell(input.dates)
    input.dates  = 'all';
end

% Pipelines
if ~isfield(input,'useNWB') || isempty(input.useNWB),           input.useNWB      = true; end
if ~isfield(input,'ExtractData') || isempty(input.ExtractData), input.ExtractData = true; end

% Plots
if ~isfield(input,'plots'),   input.plots     = []; end
if ~isfield(input,'test_ch'), input.test_ch   = []; end

%% Set Dependencies
cd(input.mainfolder)
addpath functions\
addpath toolboxes\fieldtrip_light
addpath toolboxes\Deuteron
addpath toolboxes\Viewer
addpath(genpath('toolboxes\multitaper_prerau'))
ft_defaults

%% Send input to base workspace
assignin('base','input', input);
end