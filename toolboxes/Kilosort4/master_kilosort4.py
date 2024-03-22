# This file, together with a customized 'parameters.py' are copied to the kilosort enviroment 
# folder before the call happes in MATLAB.

# MATLAB call lives in 'master_kilosort4.m'. This call is of the shape:
# pyrunfile("master_kilosort4.py \                                          % master script (this file)
#           'C:\code\miniconda3\envs\kilosort\Lib\site-packages\kilosort' \ % enviroment place
#           '...\data\preprocessing\913\20240219' \                         % bin file absolute location
#           '32' \                                                          % number of channels
#           '...\analysisCode\chanMapE32-S2_linearized_DeutSN11.mat' \      % absolute path to probe map
#           "                                                               % end of string
#
# TODO: optimize all this.

import sys
import os

# Take matlab argument to change working directory
os.chdir(sys.argv[1])

# Here we can start getting kilosort packages
from kilosort import run_kilosort, DEFAULT_SETTINGS

# Set the settings. DEFAULT_SETTINGS can be customized in 'parameters.py' 
settings = DEFAULT_SETTINGS

# Some common settings can be overrided by passing them as arguments in MATLAB
settings['data_dir'] = sys.argv[2]
settings['n_chan_bin'] = int(sys.argv[3])

# Run the thing!
ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate = \
    run_kilosort( \
    settings=settings, \
    probe_name=sys.argv[4] \
    )
