function master_kilosort4(input, varargin)
% Implementation of a MATLAB wrapper to the newly developed Kilosort4, which runs
% completely under python. To be used with the API version, programatically, 
% during a regular session processing.
%
% Fundamental algorithm from MOUSELAND KILOSORT GITHUB. Cite the toolbox and paper:
% https://github.com/MouseLand/Kilosort
% General documentation: https://kilosort.readthedocs.io/en/latest/
%
% This function and the python wrapper by Jesus J. Ballesteros 08.2024

%% INSTALL Python requirements and kilosort4
%  1. To be able to use Kilosort4 at all. This will be setup once per
%  computer and, in principle, not anymore.
%   - Install a Anaconda distribution, if non-existing. For example, miniconda:
%       https://docs.anaconda.com/free/miniconda/miniconda-install/
%   - Create a python enviroment. 
%       Open a conda prompt which has 'conda' for python 3 in the path:
%       For that, search in windows start for 'conda' and open the CONDA prompt
%       In there, type (without '): 'conda create --name kilosort python=3.9'
%   - Activate the python enviroment. To get into the enviroment, so all
%       packages and libraries can only be used within itself. In conda prompt
%       type: 'conda activate kilosort'
%   - Install Kilosort from source. Type one of the following:
%       'python -m pip install .' 
%       'python -m pip install .[gui]' (recommended)
%   - Install GPU version for pytorch. Type, one after the other:
%       'pip uninstall torch'
%       'conda install pytorch pytorch-cuda=11.7 -c pytorch -c nvidia'

%% USE
% We have created a python script, containing all necessary instructions to run 'kilosort_run'. Basically:
% >>
%   import sys, os
%   move to working directory (enviroment)
%   from kilosort import run_kilosort, DEFAULT_SETTINGS
%   settings = DEFAULT_SETTINGS
%   settings['data_dir'] = 'project/data/preprocessed/animal/session/binfile.bin'
%   settings['n_chan_bin'] = numChannels
%   ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate = run_kilosort(settings=settings, probe_name='chanMap*.mat')
% <<

%% DESCRIPTION 
% 'run_kilosort' is the main function. From KS documentation:
%  kilosort.run_kilosort.run_kilosort(
%   settings, (dic) Specifies a number of configurable parameters used throughout the spike sorting pipeline. See kilosort/parameters.py for a full list of available parameters. 
%       'n_chan_bin': def 385. MUST be specified here (all other settings are optional).
%       'fs': def 30000. Sampling frequency.
%       'batch_size': def 60000. Number of samples included in each batch of data.
%       'nblocks'; def 1. Number of non-overlapping blocks for drift correction (additional nblocks-1 blocks are created in the overlaps).
%       'Th_universal': def 9. Spike detection threshold for universal templates. Th(1) in previous versions of Kilosort.
%       'Th_learned': def 8. Spike detection threshold for learned templates. Th(2) in previous versions of Kilosort.
%       'tmin': def 0. Time in seconds when data used for sorting should begin. By default, begins at 0 seconds.
%       'tmax': def 'inf'. Time in seconds when data used for sorting should end. By default, ends at the end of the recording.
%       'nt':  def 61. Number of samples per waveform. Also size of symmetric padding for filtering.
%       'artifact_threshold': def 'inf'. If a batch contains absolute values above this number, it will be zeroed out under the assumption that a recording artifact is present. By default, the threshold is infinite (so that no zeroing occurs).
%       'nskip': def 25. Batch stride for computing whitening matrix.
%       'whitening_range': def 32. Number of nearby channels used to estimate the whitening matrix.
%       'binning_depth': def 5. For drift correction, vertical bin size in microns used for 2D histogram.
%       'sig_interp': def 20. For drift correction, sigma for interpolation (spatial standard deviation). Approximate smoothness scale in units of microns.
%       'nt0min': def 'None'. Sample index for aligning waveforms, so that their minimum or maximum value happens here. Default of 20.
%       'dmin': def 'None'. Vertical spacing of template centers used for spike detection, in microns. Determined automatically by default.
%       'dminx': def 'None'. Horizontal spacing of template centers used for spike detection, in microns. Determined automatically by default.
%       'min_template_size': def 10. Standard deviation of the smallest, spatial envelope Gaussian used for universal templates.
%       'template_sizes': def 5. Number of sizes for universal spike templates (multiples of the min_template_size).
%       'nearest_chans': def 10. Number of nearest channels to consider when finding local maxima during spike detection.
%       'nearest_templates': def 100. Number of nearest spike template locations to consider when finding local maxima during spike detection.
%       'templates_from_data': def 'true'.  Indicates whether spike shapes used in universal templates should be estimated from the data or loaded from the predefined templates.
%       'n_templates':  def 6. Number of single-channel templates to use for the universal templates (only used if templates_from_data is True).
%       'n_pcs': def 6. Number of single-channel PCs to use for extracting spike features (only used if templates_from_data is True).
%       'Th_single_ch': def 6. For single channel threshold crossings to compute universal-templates. In units of whitened data standard deviations. 
%       'acg_threshold': def 0.20. Fraction of refractory period violations that are allowed in the ACG compared to baseline; used to assign "good" units.
%       'ccg_threshold': def 0.25. Fraction of refractory period violations that are allowed in the CCG compared to baseline; used to perform splits and merges.
%       'cluster_downsampling': def. 20. Inverse fraction of nodes used as landmarks during clustering (can be 1, but that slows down the optimization). 
%       'cluster_pcs': def 64.  Maximum number of spatiotemporal PC features used for clustering.
%       'duplicate_spike_bins': def 15.  Number of bins for which subsequent spikes from the same cluster are assumed to be artifacts. A value of 0 disables this step.        
%   probe = None, (dict; optional.) – A Kilosort4 probe dictionary, as returned by kilosort.io.load_probe.
%   probe_name = None, (str; optional.) – Filename of probe to use, within the default PROBE_DIR. Only include the filename without any preceeding directories. Will ony be used if probe is None. Alternatively, the full filepath to a probe stored in any directory can be specified with settings = { probe_path : …}.
%   filename = None, (str or Path; optional.) – Full path to binary data file. If specified, will also set data_dir = filename.parent.
%   data_dir = None, (str or Path; optional.) – Specifies directory where binary data file is stored. Kilosort will attempt to find the binary file. This works best if there is exactly one file in the directory with a .bin, .bat, .dat, or .raw extension. Only used if filename is None. Also see kilosort.io.find_binary.
%   file_object = None,  (array-like file object; optional.) – Must have  shape  and ‘dtype’ attributes and support array-like indexing (e.g. [:100,:], [5, 7:10], etc). For example, a numpy array or memmap. Must specify a valid filename as well, even though data will not be directly loaded from that file.
%   results_dir = None, (str or Path; optional.) – Directory where results will be stored. By default, will be set to data_dir /  kilosort4 .
%   data_dtype = None, (str or type; optional.) – dtype of data in binary file, like  int32  or np.uint16. By default, dtype is assumed to be ‘int16’.
%   do_CAR = True, (bool; default=True.) – If True, apply common average reference during preprocessing (recommended).
%   invert_sign = False, (bool; default=False.) – If True, flip positive/negative values in data to conform to standard expected by Kilosort4.
%   device = None, (torch.device; optional.) – CPU or GPU device to use for PyTorch calculations. By default, PyTorch will use the first detected GPU. If no GPUs are detected, CPU will be used. To set this manually, specify device = torch.device(<device_name>).
%   progress_bar = None, (tqdm.std.tqdm or QtWidgets.QProgressBar; optional.) – Used by sorting steps and GUI to track sorting progress. Users should not need to specify this.
%   save_extra_vars = False, (bool; default=False.) – If True, save tF and Wall to disk after sorting.
% )

%% Defaults
if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

% Channelmap files are to be found under '\analysisCode' !!!
if ~isfield(opt,'KSchanMapFile') || isempty(opt.KSchanMapFile)
    opt.KSchanMapFile   = ls(fullfile(input.analysisCode, 'chanMap*.mat')); % load the map in this folder if not defined
end 
    
%% Check existence of previous KS4 results
if isfolder(fullfile(opt.FolderProcDataMat, 'kilosort4'))
    content = dir(fullfile(opt.FolderProcDataMat, 'kilosort4'));
    if length(content)>10
        disp('It seems like KS4 was already ran in this session. To re-run, delete the folder under preprocessing.')
        return
    end
end
clear content

%% Set up Kilosort enviroment
% Call enviroment status
%pe = pyenv;
pe = pyenv(Version=[input.KSpython,'python.exe'], ExecutionMode="OutOfProcess");

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
        % error('Something went wrong while reloading Python Interpreter. Restart Matlab.')        
        py.list; % a call to restart the Interpreter
        pe = pyenv; % Recall enviroment status
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
command.script = "master_kilosort4.py"; % Our script that wraps the call to kilosort_run
command.s1 = " '"; % To introduce the necessary 's before the argument.
command.s2 = "'"; % To introduce the necessary 's after the argument.
command.var1 = string(input.KSpyfolder); % var1 is the absolute path to the kilosort library in the python enviroment
command.var2 = string(fullfile(opt.FolderProcDataMat, [opt.SavFileName, '.bin'])); %,  % var2 is the absolute path to the .bin file has been created
command.var3 = string(opt.numChannels); % Give number of channels as string Will convert to int within python script)
command.var4 = append(input.analysisCode, opt.KSchanMapFile); % Absolute path to the probe map.
% command.var5 = string();

% Put all strings together for a full argument
command.full = append(command.script, ...
    command.s1, command.var1, command.s2, ...
    command.s1, command.var2, command.s2, ...
    command.s1, command.var3, command.s2, ...
    command.s1, command.var4, command.s2 ... % If you add more varX, this one needs a comma at the end, before '...'
    ); % Add more 'command.s1, command.varX, command.s2 ...' for new variables, and make sure you collect them in the python script

%% RUN
% Make sure we use the project's parameters
cd(input.analysisCode)
projfiles = string(ls("*.py"));

% Some projects might use more than one probe. CAREFUL!
if length(projfiles)>3 % if project uses only one probe, there should be no more than 3 .py files
    if contains(opt.KSchanMapFile, 'S2') % This would depend on the specific probes used
        copyfile(string(fullfile(input.analysisCode,projfiles{3})), string(fullfile(input.KSpyfolder, 'parameters.py')),'f');
    elseif contains(opt.KSchanMapFile, 'Poly3') % This would depend on the specific probes used
        copyfile(projfiles{2}, [input.KSpyfolder '\parameters.py'],'f');
    else
        error('Your specific configuration for Kilosort4 does not seem to be listed.')
    end
else
    % One single probe would mean there is 3 files, being the parameters' the 2nd one.
    copyfile(projfiles{3}, input.KSpyfolder,'f');
end

copyfile(projfiles{1},input.KSpyfolder,'f'); 

% Move to the kilosort enviroment working directory
cd(input.KSpyfolder)

% Clear cache
if isfolder("__pycache__")
    rmdir __pycache__ s
end

% Run the wrapper script with the given arguments
pyrunfile(command.full)

% Terminate python process
terminate(pyenv)
