function intan2NWB_neuroconv(input, opt)
% Makes sure there is no '.nwb' files in directory. Then, runs the wrapper
% for the INTANtoNWB tool.
%
% INPUT:    
%   input        struct, general inputs to script. Needs the path to python folder.
%   opt.sessions struct, info for sessions: folder, name and number of sessions.
%   opt.ss       int, numeral of processing session
%
% OUTPUT:   no explicit output.
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
    if isempty(files) || files.bytes < 1e6
        % Warn about file being process.
        disp('- Will convert session to NWB format. This may take a moment.');

        % Navigate to python folder
        cd(input.NCfolder);
        
        %% Check NeuroConv enviroment
        terminate(pyenv)

        % Call enviroment status
        pe = pyenv(Version=fullfile(input.NCfolder,'python.exe'), ExecutionMode="OutOfProcess");
        
        % Check if pyenv is set, or kill any residual process running
        if pe.ExecutionMode && pe.Status > 0
            terminate(pyenv)
            pe = pyenv;
        
            if pe.Status == "Terminated"
                pe = pyenv('ExecutionMode', 'OutOfProcess');
                py.list;
                pe = pyenv;
            else
                py.list;
                pe = pyenv;
            end
        end
        
        % Display current status if Loaded
        if pe.Status == "Loaded"
            disp(append('Python enviroment set as version: ', pe.Version))
        else
            error('Something went wrong with the Python enviroment setup.')
        end
        
        %% Prepare argument to send to the python script
        % Arguments:
        %   var1  absolute path to the INTAN header file (info.rhd)
        %   var2  absolute path for the output .nwb file
        %   var3  absolute path to the project nwb_metadata.yaml in analysisCode\
        command.script = "master_neuroconv.py";
        command.s1     = " '";
        command.s2     = "'";
        command.var1   = string(fullfile(opt.PathRaw, 'info.rhd'));
        command.var2   = string(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.nwb']));
        command.var3   = string(fullfile(input.analysisCode, 'nwb_metadata.yaml'));

        command.full = append(command.script, ...
            command.s1, command.var1, command.s2, ...
            command.s1, command.var2, command.s2, ...
            command.s1, command.var3, command.s2  ...
            );

        %% Copy wrapper script to the NC environment folder, then run from there.
        % pyrunfile resolves scripts relative to the current directory, so we
        % copy master_neuroconv.py next to python.exe (same pattern as
        % master_kilosort4.m copies its script to input.KSpyfolder).
        copyfile(fullfile(input.analysisCode, 'master_neuroconv.py'), input.NCfolder, 'f');

        % Clear any stale bytecode cache
        if isfolder("__pycache__")
            rmdir __pycache__ s
        end

        % Run the conversion
        disp('- Conversion in progress...');
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