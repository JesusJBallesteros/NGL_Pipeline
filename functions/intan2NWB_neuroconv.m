function intan2NWB_neuroconv(input, opt)
% Makes sure there is no '.nwb' files in directory. Then, runs the wrapper
% for the INTANtoNWB tool.
% INPUT:    
%   input        struct, general inputs to script. Needs the path to python folder.
%   opt.sessions struct, info for sessions: folder, name and number of sessions.
%   opt.ss       int, numeral of processing session
% OUTPUT:   none explicit.
%           It generates a new file with extension .nwb in the /processed folder
%
% Version 16.01.2025 Jesus

% 00. Tell where the python folder with 'IntanToNWB' scripts is
% Added to 'set_default' now. If not found, check that it is working.
% IntanToNC_folder = input.NCfolder;

% Check for .nwb files in output directory and come back
currdir = pwd;
cd(opt.FolderProcDataMat)
    files = dir('*.nwb'); 

    %% 01. If there is none, proceed
    if isempty(files)
        % Warn about file being process.
        disp('- Will convert session to NWB format. This may take a moment.');

        % Navigate to python folder
        cd(input.NCfolder);
        
        %% Check NeuroConv enviroment
        % Call enviroment status
        pe = pyenv('Version', input.NCfolder);
        
        % Check if pyenv is set, or kill any residual process running
        if pe.ExecutionMode && pe.Status > 0
            terminate(pyenv) % Terminate process
            pe = pyenv; % Recall enviroment status
        
            % Proceed to start enviroment
            if pe.Status == "Terminated"
                % And only if properl'y terminated, reset it
                pe = pyenv('ExecutionMode', 'OutOfProcess');
                py.list; % a call to restart the Interpreter
                pe = pyenv; % Recall enviroment status
            else
                error('Something went wrong while reloading Python Interpreter. Restart Matlab.')
            end
        end
        
        % Display current status if Loaded
        if pe.Status == "Loaded"
            disp(append('Python enviroment set as version: ', pe.Version))
        else
            error('Something went wrong with the Python enviroment setup.')    
        end
        
        %% Prepare argument to send to the python script
        % The argument is given as a single string, so we can prepare the pieces to
        % put them all together at the end.
        command.script = "master_neuroconv.py"; % Our script that wraps the call to neuroconv
        command.s1 = " '"; % To introduce the necessary 's before the argument.
        command.s2 = "'"; % To introduce the necessary 's after the argument.
        command.var1 = "C:\Code\miniconda3\envs\neuroconv\Lib\site-packages\neuroconv"; % var1 is the absolute path to the kilosort library in the python enviroment
        command.var2 = string(fullfile(opt.PathRaw, 'info.rhd')); %,  % var2 is the absolute path to the INTAN header file
        command.var3 = string(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '.nwb'])); %,  % var3 is the absolute path to the .bin file has been created
        % command.var4 = string();
        
        % Put all strings together for a full argument
        command.full = append(command.script, ...
            command.s1, command.var1, command.s2, ...
            command.s1, command.var2, command.s2, ...
            command.s1, command.var3, command.s2 ... % If you add more varX, this one needs a comma at the end, before '...'
            ); % Add more 'command.s1, command.varX, command.s2 ...' for new variables, and make sure you collect them in the python script

%         cd(input.NCfolder)
        
        % Clear cache
        if isfolder("__pycache__")
            rmdir __pycache__ s
        end

        % Run the routine
        disp('- Conversion in progress...');
        % Run the wrapper script with the given arguments
        pyrunfile(command.full)
        
        % Terminate python process
        terminate(pyenv)

        cd(currdir)
    else
        % If there is already one, exit the function and continue. Warn about it
        disp('- A file in the NWB format has been found for this session. Skipping.');
        return
    end    
end