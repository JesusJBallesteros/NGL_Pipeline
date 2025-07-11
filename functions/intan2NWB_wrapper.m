function intan2NWB_wrapper(input, opt)
% Makes sure there is no '.nwb' files in directory. Then, runs the wrapper
% for the INTANtoNWB tool.
% INPUT:    
%   input    struct, general inputs to script. Needs the path to python folder.
%   sessions struct, info for sessions: folder, name and number of sessions.
%   ss       int, numeral of processing session
% OUTPUT:   none explicit.
%           It generates a new file with extension .nwb in the /processed folder
%
% Version 07.03.2023 Jesus

% 00. Tell where the python folder with 'IntanToNWB' scripts is
% Added to 'set_default' now. If not found, check that it is working.
IntanToNWB_folder = input.pyfolder;

% Check for .nwb files in directory
files = dir('*.nwb'); 

    %% 01. If there is none, proceed
    if isempty(files)
        % Warn about file being process.
        disp('- Will convert session to NWB format. This may take a moment.');

        % List all files in origin.
        files = dir('*.dat'); 
        files = [files; dir('*.rhd')];

        % Copy one by one to Python folder.
        for i=1:length(files)
            fprintf('- Copying file %d of %d.\n', i , length(files));
            [copy.status, copy.msg] = copyfile(files(i).name, IntanToNWB_folder);
        end

        % Navigate to python folder
        cd(IntanToNWB_folder);

        %% Run IntanToNWB python routine using the bypasser.
        % WHAT IT IS: Substitute input layer to avoid using Jupyter Notebook GUI.
        % It bypasses 'IntanToNWB.ipynb' and 'ConverterUI.py' scripts, so they are
        % not longer necessary.
        %
        % WHAT IT DOES: It navigates to the Python folder, runs the routine, copies the .nwb file
        % back to the original data folder and navigates back to it.
        %
        % TODO: Add inputs from Matlab and pass them to the script for easyness
        %
        % WHAT IT CONTAINS: 'my_IntanToNWB.py' script with default settings as:
        %   intan_filename="info.rhd"       % Default header filename when multiple files are created with INTAN
        %   nwb_filename="info.nwb"         % Default new file name. Suggested to be changed afterwards.
        %   session_description=None        % Do not add descriptions, since the INTAN python wrapper seems to not really like it.
        %   blocks_per_chunk=2000           % This can be modified to change how much data is proccesed per bulk.
        %   use_compression=False           % No compression
        %   compression_level=4             % Left default, not used if above False
        %   lowpass_description=None        % Not defined
        %   highpass_description=None       % Not defined
        %   merge_files=None                % No merging header files
        %   subject=None                    % No metadata info at subject level
        %   manual_start_time=datetime.now()% Use current time
        %
        % The Python Script is reproduced here. Copy it in a new file and add it to
        % the folder where the IntanToNWB scripts are as '.py' file. Modify the
        % name in pyrunfile("new_name.py"). 
        % 'my_IntanToNWB.py' file starts below (remove % % symbols when pasting):
        % % from ConvertIntanToNWB import *
        % % from datetime import datetime
        % % """ Callback function for when the user begins conversion.
        % % 
        % % Parameters
        % % ----------
        % % None
        % % 
        % % Returns
        % % -------
        % % None
        % % """
        % %         
        % % # Get NWB settings and pass these to NWB conversion
        % % intan_filename="info.rhd"
        % % nwb_filename="info.nwb"
        % % session_description=None
        % % blocks_per_chunk=2000
        % % use_compression=False
        % % compression_level=4
        % % lowpass_description=None
        % % highpass_description=None
        % % merge_files=None
        % % subject=None
        % % manual_start_time=datetime.now()
        % % 
        % % convert_to_nwb(intan_filename=intan_filename,
        % %                nwb_filename=nwb_filename,
        % %                session_description=session_description,
        % %                blocks_per_chunk=blocks_per_chunk,
        % %                use_compression=use_compression,
        % %                compression_level=compression_level,
        % %                lowpass_description=lowpass_description,
        % %                highpass_description=highpass_description,
        % %                merge_files=merge_files,
        % %                subject=subject,
        % %                manual_start_time=manual_start_time)
        % File ends above.
        %
        % By Jesus J. Ballesteros 09.2022
        disp('- Checking for Python enviroment...');
        % Verify Python Configuration
        pe = pyenv; disp(pe);
        
        % Run the routine
        disp('- Conversion in progress...');
        pyrunfile("NGL_IntanToNWB.py");

        %% Find and move the new .nwb file to processed data folder
        nwbfile = dir('*.nwb'); 
        if ~isempty(nwbfile)
            disp('- Done! Moving NWB file back to original folder...');
            movefile(nwbfile.name, opt.FolderProcDataMat);
        else
            disp('- Something went wrong. Cannot find NWB files.');
            return
        end
        
        % Delete data files not necessary anymore
        delete(files(:).name);
        
        % Navigate to saving folder
        cd(opt.FolderProcDataMat);
        % Edif file name.
        movefile(nwbfile.name,strcat([opt.SavFileName, '.nwb']))

    else
    %% If there is any, exit the function and continue
        % Warn about existing .nwb files
        disp('- A file in the NWB format has been found for this session. Skipping.');
        return
    end    
end