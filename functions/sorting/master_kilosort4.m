function master_kilosort4(input, varargin)
% master_kilosort4  MATLAB -> Python wrapper for Kilosort 4 spike sorting.
%
% PURPOSE:
%   Calls Kilosort 4 (which runs entirely in Python) from MATLAB via pyrunfile.
%   Activates the KS4 Python environment, copies the wrapper script and
%   parameters file into it, then executes Kilosort on the session's .bin file.
%   KS4 output (templates, spike times, cluster assignments) is written to
%   opt.KSfolder.
%
% USAGE:
%   master_kilosort4(input, opt)
%   master_kilosort4(input)    % opt defaults used
%   Called from NGL01_Main stage 04; do not call directly.
%   For multi-area runs, NGL01_Main calls this once per unique area,
%   overriding opt.KSfolder and opt.KSchanMapFile before each call.
%
% INPUTS:
%   input  - struct from set_default + prepforsession; relevant fields:
%              .KSpythonExe     path to Kilosort 4 Conda environment
%              .KSpyfolder      path to kilosort package inside that env
%              .analysisCode    path to analysis code folder (chanMap location)
%   opt    - (optional) struct; relevant fields:
%              .FolderProcDataMat  folder containing the .bin file
%              .KSfolder           output folder for KS4 results
%              .SavFileName        session name
%              .numChannels        number of electrode channels (total in .bin)
%              .KSchanMapFile      chanMap filename ('' = linear) or full path
%
% OUTPUT:
%   Kilosort 4 output directory at opt.KSfolder containing:
%     spike_times.npy, spike_templates.npy, templates.npy,
%     cluster_group.tsv, params.py, and associated files.
%
% REQUIREMENTS:
%   - Kilosort 4 installed in the KS Python environment:
%       conda activate kilosort && pip install kilosort[gui]
%   - GPU (CUDA) strongly recommended for performance.
%   - Only one pyenv can be active per MATLAB session; restart MATLAB between
%     runs if switching Python environments.
%
% CITE:
%   Pachitariu et al. (2024). Kilosort4: https://github.com/MouseLand/Kilosort
%
% MATLAB wrapper: Jesus J. Ballesteros, 08.2024
% Last modified 08.05.2026 (Jesus)
%
%% INSTALL Python requirements and kilosort4
%  1. To be able to use Kilosort4 at all. This will be setup once per
%  computer and, in principle, not anymore.
%   - Install a Anaconda distribution, if non-existing. For example, miniconda:
%       https://docs.anaconda.com/free/miniconda/miniconda-install/
%   - Create a python enviroment.
%       In conda prompt type: 'conda create --name kilosort python=3.9'
%   - Activate and install Kilosort:
%       'conda activate kilosort'
%       'python -m pip install .[gui]'
%   - Install GPU pytorch:
%       'pip uninstall torch'
%       'conda install pytorch pytorch-cuda=11.7 -c pytorch -c nvidia'
%
%% USE
% Python script wrapper calls run_kilosort with:
%   settings['n_chan_bin'] = numChannels
%   filename = /path/to/session.bin
%   probe_name = /path/to/chanMap.mat
%   results_dir = opt.KSfolder
%
%% DESCRIPTION
% 'run_kilosort' key parameters (see kilosort/parameters.py for full list):
%   'n_chan_bin'    - MUST be set; total channels in binary file (incl. disconnected)
%   'fs'           - sampling rate (default 30000)
%   'Th_universal' - spike detection threshold for universal templates (default 9)
%   'Th_learned'   - spike detection threshold for learned templates (default 8)
%   results_dir    - output directory; defaults to data_dir/kilosort4 if not set

%% Defaults
if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

% Channelmap files are found under analysisCode unless a full path is given.
if ~isfield(opt,'KSchanMapFile') || isempty(opt.KSchanMapFile)
    opt.KSchanMapFile = ls(fullfile(input.analysisCode, 'chanMap*.mat'));
end

%% Check existence of previous KS4 results
% Use opt.KSfolder so per-area checks work correctly in multi-area mode.
if isfolder(opt.KSfolder)
    content = dir(opt.KSfolder);
    if length(content) > 10
        disp(['KS4 results already present in: ', opt.KSfolder, ' -- skipping. Delete folder to re-run.'])
        return
    end
end
clear content

%% Set up Kilosort enviroment
pe = pyenv(Version=[input.KSpythonExe,'python.exe'], ExecutionMode="OutOfProcess");

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

if pe.Status == "Loaded"
    disp(append('Python enviroment set as version: ', pe.Version))
else
    error('Something went wrong with the Python enviroment setup.')
end

%% Prepare argument to send to the python script
command.script = "master_kilosort4.py";
command.s1 = " '";
command.s2 = "'";

% var1: absolute path to kilosort library in python env
command.var1 = string(input.KSpyfolder);

% var2: absolute path to the .bin file
command.var2 = string(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '.bin']));

% var3: total channel count in the .bin (n_chan_bin for KS4)
command.var3 = string(opt.numChannels);

% var4: absolute path to chanMap/probe file.
% If opt.KSchanMapFile is already a full path (multi-area per-area file), use directly.
% Otherwise join with analysisCode (single-area default).
if isfile(opt.KSchanMapFile)
    command.var4 = string(opt.KSchanMapFile);
else
    command.var4 = string(fullfile(input.analysisCode, opt.KSchanMapFile));
end

% var5: results directory for this KS run.
% Explicit opt.KSfolder enables per-area output folders in multi-area mode.
command.var5 = string(opt.KSfolder);

command.full = append(command.script, ...
    command.s1, command.var1, command.s2, ...
    command.s1, command.var2, command.s2, ...
    command.s1, command.var3, command.s2, ...
    command.s1, command.var4, command.s2, ...
    command.s1, command.var5, command.s2  ...
    );

%% RUN
% Select and copy parameters file into the KS4 environment folder.
% Priority order:
%   1. parameters_<Area>.py  — area-specific settings (e.g. parameters_NCL.py)
%   2. parameters.py         — single shared parameters file
%   3. def_parameters.py     — project default, fallback
% The area label is the last path component of opt.KSfolder (e.g. 'NCL', 'STR').
% In single-area mode opt.KSfolders is absent and areaLabel stays empty,
% so the area-specific lookup is skipped.
cd(input.analysisCode)

% Determine area label
if isfield(opt, 'KSfolders')
    [~, areaLabel] = fileparts(opt.KSfolder);
else
    areaLabel = '';
end

% Find the best-matching parameters file
paramFile = '';
if ~isempty(areaLabel)
    candidate = fullfile(input.analysisCode, sprintf('parameters_%s.py', areaLabel));
    if isfile(candidate), paramFile = candidate; end
end
if isempty(paramFile)
    candidate = fullfile(input.analysisCode, 'parameters.py');
    if isfile(candidate), paramFile = candidate; end
end
if isempty(paramFile)
    candidate = fullfile(input.analysisCode, 'def_parameters.py');
    if isfile(candidate), paramFile = candidate; end
end
assert(~isempty(paramFile), 'NGL:noParameters', ...
    'No parameters .py file found in %s. Expected parameters_<Area>.py, parameters.py, or def_parameters.py.', ...
    input.analysisCode);

% Copy parameters file into KS env (always as "parameters.py")
copyfile(paramFile, fullfile(input.KSpyfolder, 'parameters.py'), 'f');
fprintf('KS4 parameters: %s\n', paramFile);

% Copy the MATLAB->Python run script
copyfile(fullfile(input.analysisCode, 'master_kilosort4.py'), input.KSpyfolder, 'f');

cd(input.KSpyfolder)
if isfolder("__pycache__")
    rmdir __pycache__ s
end

pyrunfile(command.full)

terminate(pyenv)

%% Post-run: relocate KS4 output to the intended area folder if needed.
% Some KS4 versions ignore results_dir and always write to
% <data_dir>\kilosort4. When that happens and this is a multi-area run
% (detectable by opt.KSfolders being present), move the output folder.
% Conditions for the move:
%   - the default kilosort4\ folder exists and has actual files
%   - opt.KSfolder differs from that default location (area run)
%   - the intended area folder exists but is empty (pre-created by prepforsession)
defaultKSout = fullfile(opt.FolderProcDataMat, 'kilosort4');
if isfield(opt, 'KSfolders') && isfolder(defaultKSout) && ~strcmp(defaultKSout, opt.KSfolder)
    defaultContent = dir(defaultKSout);
    areaContent    = dir(opt.KSfolder);
    if length(defaultContent) > 3 && length(areaContent) <= 4
        rmdir(opt.KSfolder);                  % remove the pre-created empty area folder
        movefile(defaultKSout, opt.KSfolder); % rename kilosort4\ -> <Area>\ ("move")
        fprintf('KS4 output relocated: kilosort4/ -> %s\n', opt.KSfolder);
    end
end

end
