%% NGL00_Prep
% Running this script will create and set a complete folder structure in the
% PC, to store and process new projects. It is necessary to have the right 
% folder system, so the following scripts can find the data and store it
% according to the lab standards.
% 
% JESUS 04.01.2023

%% The folder system will be created under 'datadrive:\studyname\'
% Data is stored in a main HD or SSD unit. 
if ~contains(datadrive,':\')
    datadrive = [datadrive ':\'];
end

% Project folder.
projectfolder = fullfile(datadrive, studyname);

if ~exist(projectfolder,"dir")
    txt = sprintf('Folder system for project "%s" will be created. \n', studyname);
    fprintf(txt);

    % Create project folder and change current directory to it.
    mkdir(projectfolder);
    cd(projectfolder);
    
    % Create 'readme.txt' file.
    fileID = fopen('readme.txt','w');
    
    % Fill content into txt file and close it.
    fprintf(fileID,'%s \r\n', readmecontent);
    fclose(fileID);
    
    % Create first level subfolders.
    mkdir(projectfolder, 'analysisCode');
    mkdir(projectfolder, 'data');
    mkdir(projectfolder, 'manuscript');
    mkdir(projectfolder, 'paradigmCode');
    mkdir(projectfolder, 'training');
    
    % Create data secondary subfolders.
    cd(fullfile(projectfolder, 'data'))
    mkdir('analysis');
    mkdir('preprocessing');
    mkdir('raw');
    mkdir('spikesorted');
    mkdir('trialsorted');
    mkdir('behaviour');
    
    % Update result
    txt = sprintf('Folder system for project "%s" created. Done. \n', studyname);    
    fprintf(txt);
    
    clear fileID txt
else
    % Update result
    txt = sprintf('Folder system for project "%s" already exists (or not). But nothing changed. \n', studyname);
    fprintf(txt);

    clear txt
end

clear projectfolder