%% Jesus' Pipeline to analyze FieldTrip formatted files
% TO BE redacted
% TO BE DONE
% One of the first functions here will split the continous dataset into
% trials, according to a pre-fixed length or a trial scheme given by the
% user.
% 
% Last modified Jesus 07.03.2023

% If not provided, a promp will ask for them or 'set_default' will use the 
% defaults. It will also put them in the correct format if an incorrect one 
% was given.
input.mainfolder = 'C:\Code\Scripts\ephys-data-pipeline';   % Default: 'C:\Code\Scripts\ephys-data-pipeline';
input.datafolder = 'D:\Experiments\';                       % Default: 'D:\Experiments\';
input.animal     = '420';                                   % i.e '420' or FAT;
input.processed  = 'processed';                             % Default: 'processed'. Subfolder to be created in session folder
input.dates       = {'20230220_Deut'} ;                     % can be left empty, 'all', or a list like:
    % {'20230217_01'   '20230217_02'  ...
    %  '20230220_Deut' '20230221_Deut'...
    %  '20230222_Deut' '20230220_Int' ...
    %  '20230221_Int'  '20230222_Int'    }; 

% These apply to FieldTrip-ready .mat files, only.
input.plots      = []; % An logic array of 0/1s, to ask for specific plots. See details.
input.test_ch    = []; % An array of numerals for channels to plot.

%% Options for the different wrappers
opt = struct();
    opt.parsing  = true;
    opt.length   = 5; % Can be a fix value, or can use 'Events'

%% 00. Check inputs, set defaults and dependencies.
set_default(input);

%% 01. Find and list sessions. 
% Read requested sessions from specified animal folder.
% This 'sessions' variable can be used for summary, book keeping and
% debugging at the end of the pipeline. But it will not be saved
% automatically. % smt TODO?
sessions = findSessions(input);

%% 02. Loop sessions to process
for ss = 1:sessions.nSessions
    % Progress report
    txt = sprintf('\n --> Session %d out of %d. Session name: %s \n', ss, sessions.nSessions, sessions.list(ss).name);
    fprintf(txt);

    %% 03. Check file type and versions
    % Navigate to session raw data folder.
    cd(fullfile(sessions.folder,sessions.list(ss).name,input.processed));
    












end
