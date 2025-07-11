import sys
import os

# Take matlab arguments (ksdir, datadir, nchan, ...)
ksdir = sys.argv[1]

# Go to ksdir
os.chdir(ksdir)
from kilosort import run_kilosort, DEFAULT_SETTINGS

# Set the settings
settings = DEFAULT_SETTINGS
#settings['data_dir'] = sys.argv[2]
settings['n_chan_bin'] = int(sys.argv[3])

# Run
ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = \
    run_kilosort(settings=settings, \
                 probe_name=sys.argv[4], \
                 filename=sys.argv[2] \
                )
