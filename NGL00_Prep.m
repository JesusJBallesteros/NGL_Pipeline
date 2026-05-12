%% NGL00_Prep
% Running this script at the time of setting up a project will create and 
% set a complete folder structure in the PC, to store and process data within.
% It is necessary to have the right folder system, so the following scripts 
% can find the data and store it according to the lab standards.
% A call during other NGLXX scripts will collect and sort the user inputs into
% proper structures expected by the pipeline.
% 
% Last Version. 12.05.2026

%% The folder system will be created under 'datadrive:\studyname\'
% Data is stored in a main HD or SSD unit. 
if exist('datadrive','var')
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
        warning('NOW is a good time to check your CONFIG files. They should go into your analysisCode folder')
        
        clear fileID txt
    else
        % Update result
        txt = sprintf('Folder system for project "%s" located. \n', studyname);
        fprintf(txt);
        warning('ALWAYS check your CONFIG files. They should be inside your analysisCode folder')
    
        clear txt
    end
end

% At the beggining of any other NGLXX script, the existence of opt and input
% structs will be checked.
if exist('opt','var')
    if ~exist("input","var")
        input = struct( 'datadrive' , datadrive , ...
                        'studyName' , studyname , ...
                        'subjects'  , [], ...
                        'dates'     , [], ...
                        'Areas'     , []        );
        input.dates    = dates;
        input.subjects = subjects;
        input.Areas    = areas;
    end

    clear areas subjects dates datadrive studyname
end
 clear projectfolder readmecontent