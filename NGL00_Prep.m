%% NGL00_Prep
% Running this script will create and set a complete folder structure in the
% PC, to store and process new projects. It is necessary to have the right 
% folder system, so the following scripts can find the data and store it
% according to the lab standards.
%
% Last updated 
% JESUS 27.10.2023

%% All data for a project is stored in a main HD or SSD unit. 
% Define it here:
datadrive = 'F:'; % Drive unit in local PC.

%% All projects should have a name. 
% Define it here:
studyName = 'SPPtest'; % A descriptive and unique name to the project.

%% The project folder contains a 'readme.txt' file.
% It contains details about the project. 
% Info can be added later. To add it now, just write it in 'txtcontent'.
% The default content is:
readmecontent = [   "Study name: "                                   , ...
                    "Readme date: "                                  , ...
                    "Person (1) responsible for data repository: "   , ...
                    "Person(s) responsible for study: "              , ...
                    "Hardware used: "                                , ...
                    "Related Publication(s): "                       , ...
                    "Short description of study: "                           ];

%% The folder system will be created under 'datadrive:\studyName\'
% Define project folder.
projectFolder = fullfile(datadrive, studyName);

if ~exist(projectFolder,"dir")
    txt = sprintf('Folder system for project "%s" will be created. \n', studyName);
    fprintf(txt);

    % Create project folder and change current directory to it.
    mkdir(projectFolder);
    cd(projectFolder);
    
    % Create 'readme.txt' file.
    fileID = fopen('readme.txt','w');
    
    % Fill content into txt file and close it.
    fprintf(fileID,'%s \r\n', readmecontent);
    fclose(fileID);
    
    % Create first level subfolders.
    mkdir(projectFolder, 'analysisCode');
    mkdir(projectFolder, 'data');
    mkdir(projectFolder, 'manuscript');
    mkdir(projectFolder, 'paradigmCode');
    mkdir(projectFolder, 'training');
    
    % Create data secondary subfolders.
    cd(fullfile(projectFolder, 'data'))
    mkdir('analysis');
    mkdir('preprocessing');
    mkdir('raw');
    mkdir('spikesorted');
    mkdir('trialsorted');
    
    % Update result
    txt = sprintf('Folder system for project "%s" created. Done. \n', studyName);    
    fprintf(txt);
else
    % Update result
    txt = sprintf('Folder system for project "%s" already exists. Nothing changed. \n', studyName);
    fprintf(txt);
end

% Clear workspace.
clear all