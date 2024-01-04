%% NGL00_Prep
% Running this script will create and set a complete folder structure in the
% PC, to store and process new projects. It is necessary to have the right 
% folder system, so the following scripts can find the data and store it
% according to the lab standards.
% 
% JESUS 04.01.2023

%% The folder system will be created under 'datadrive:\studyName\'
% Data is stored in a main HD or SSD unit. 
input.datadrive = [input.datadrive ':\'];
    
% Project folder.
projectFolder = fullfile(input.datadrive, input.studyName);

if ~exist(projectFolder,"dir")
    txt = sprintf('Folder system for project "%s" will be created. \n', input.studyName);
    fprintf(txt);

    % Create project folder and change current directory to it.
    mkdir(projectFolder);
    cd(projectFolder);
    
    % Create 'readme.txt' file.
    fileID = fopen('readme.txt','w');
    
    % Fill content into txt file and close it.
    fprintf(fileID,'%s \r\n', input.readme);
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
    txt = sprintf('Folder system for project "%s" created. Done. \n', input.studyName);    
    fprintf(txt);
    
    clear fileID txt
else
    % Update result
    txt = sprintf('Folder system for project "%s" already exists. Nothing changed. \n', input.studyName);
    fprintf(txt);

    clear fileID
end