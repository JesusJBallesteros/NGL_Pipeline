function intan2NWB_wrapper(sessions)
% Make sure there is no '.nwb' files in directory. Then, run the wrapper
% for the INTANtoNWB tool.
% INPUT:    sessions = struct with folder, name and number of sessions
% OUTPUT:   none explicit.
%           It generates a new file with extension .nwb in the folder of origin
% By Jesus J. Ballesteros 09.2022

%% 00. Tell where the python folder with 'IntanToNWB' scripts is
IntanToNWB_folder = 'C:\Code\Python39\IntanToNWB';

%% 01. Check for files ending in .nwb in directory
files = dir('*.nwb'); 

    %% If there is none, proceed
    if isempty(files)
        % Warn about file being process.
        disp('- Will convert session to NWB format. This may take a moment.');

        % List all files in origin.
        files = dir('*.*'); 
        files(1:2) = []; % Remove '.' and '..' outputs

        % Copy one by one.
        for i=1:length(files)
            fprintf('- Copying file %d of %d.\n', i , length(files));
            [copy.status, copy.msg] = copyfile(files(i).name,IntanToNWB_folder);
        end

        % Navigate to python folder
        cd(IntanToNWB_folder);

        %% Run IntanToNWB python routine using the bypasser.
        % Bring the new file back and delete the copies
        disp('- Check for Python enviroment...');
        convert2nwb()

        %% Find and move the new .nwb file to original data folder
        nwbfile = dir('*.nwb'); 
        if ~isempty(nwbfile)
            disp('- Done! Moving NWB file back to original folder...');
            movefile(nwbfile.name, strcat(sessions.folder,'\',sessions.name));
        else
            disp('- Something went wrong. Cannot find NWB files.');
            return
        end
        
        % Delete data files not necessary anymore
        delete(files(:).name);
        
        % Navigate back to original data folder
        cd(strcat(sessions.folder, '\', sessions.name));
        movefile(nwbfile.name,string([sessions.name + '.nwb']))

    else
    %% If there is any, exit the function and continue
        % Warn about existing .nwb files
        disp('- A file in the NWB format has been found for this session. Skipping.');
        return
    end    
end