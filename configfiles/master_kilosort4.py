import sys
import os

# Take matlab arguments: ksdir, binfile, nchan, chanmap, results_dir
# sys.argv[1] = ksdir       - path to the kilosort package inside the python env
# sys.argv[2] = binfile     - full path to the .bin recording file
# sys.argv[3] = nchan       - total number of channels in the .bin file (n_chan_bin)
# sys.argv[4] = chanmap     - full path to the probe/chanMap .mat file
# sys.argv[5] = results_dir - output folder for this KS run (enables per-area output)
ksdir       = sys.argv[1]
binfile     = sys.argv[2]
nchan       = int(sys.argv[3])
chanmap     = sys.argv[4]
results_dir = sys.argv[5]

# Go to ksdir
os.chdir(ksdir)
from kilosort import run_kilosort, DEFAULT_SETTINGS

# Set the settings
settings = DEFAULT_SETTINGS
settings['n_chan_bin'] = nchan

# Run
ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = \
    run_kilosort(settings=settings,  \
                 probe_name=chanmap, \
                 filename=binfile,   \
                 results_dir=results_dir \
                )
