function cfg = NGL_machineConfig()
    % These are the only values that should differ between lab machines.
    % Each user maintains a local copy. 
    %
    % The values in the file stored within the pipeline code should work in
    % any of the Workstations at NBL lab.

    % Where the pipeline code is stored:
    cfg.toolbox      = 'C:\Code\ephys-data-pipeline';

    % Where the 'kilosort' environment has been created 
    cfg.KSpythonExe  = 'C:\Users\ACN\miniconda3\envs\kilosort\';

    % Where the 'phy' environment has been created
    cfg.PHYpythonExe = 'C:\Users\ACN\miniconda3\envs\phy2\';

    % The exe file used within the 'neuroconv' environment
    cfg.NCpythonExe  = 'C:\Code\miniconda3\envs\neuroconv\python.exe';

end