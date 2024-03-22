function master_kilosort4(input, varargin)
% General documentation:
% https://kilosort.readthedocs.io/en/latest/

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
% We have created a python script, containing all necessary instructions
% to be ran. It looks basically as:
%
%   from kilosort import run_kilosort, DEFAULT_SETTINGS
%   settings = DEFAULT_SETTINGS
%   settings['data_dir'] = 'project/data/preprocessed/animal/session/binfile.bin'
%   settings['n_chan_bin'] = 32
%   ...
%   ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate = run_kilosort(settings=settings, probe_name='chanMap.mat')
%
% Where 'run_kilosort' is the main function. From documentation:
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

%% Input arguments check
if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Defaults, if not given.
% Config and Channelmap files are to be found under '\analysisCode'
if ~isfield(opt,'KSchanMapFile') || isempty(opt.KSchanMapFile),         opt.KSchanMapFile   = ls(fullfile(input.analysisCode, 'chanMap*.mat')); end 
if ~isfield(opt,'spkTh') || isempty(opt.spkTh),                         opt.spkTh           = 6; end 
    if opt.spkTh < 0; opt.spkTh = abs(opt.spkTh); end
if ~isfield(opt,'CAR') || isempty(opt.CAR),                             opt.CAR             = 1;  end 
    
%% Set up Kilosort enviroment
% Call enviroment status
pe = pyenv;

% Kill any running process
if pe.ExecutionMode && pe.Status > 0
    terminate(pyenv) % Terminate process
    pe = pyenv; % Recall enviroment status

    % Check for correct termination
    if pe.Status == "Terminated"
        pe = pyenv('Version', [input.KSpyfolder,'\python.exe'], 'ExecutionMode', 'OutOfProcess');
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
    disp('Something went wrong with the Python enviroment setup.')    
end

%% Prepare command
command.script = "master_kilosort4.py";
command.s1 = " '";
command.s2 = "'";
command.var1 = string(input.KSpyenv_NGL);
command.var2 = string(opt.FolderProcDataMat);
command.var3 = string(opt.numChannels);
command.var4 = append(input.analysisCode,opt.KSchanMapFile);
% command.var5 = string();

command.full = append(command.script, ...
    command.s1, command.var1, command.s2, ...
    command.s1, command.var2, command.s2, ...
    command.s1, command.var3, command.s2, ...
    command.s1, command.var4, command.s2 ...
    );

%% 02 RUN
cd(input.KSpyenv_NGL)
pyrunfile(command.full)

terminate(pyenv) % Terminate process
