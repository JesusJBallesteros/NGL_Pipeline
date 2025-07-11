Files in this folder are to be saved into your project '\analisysCode'. 

Maybe not all are necessary, that would depend on your specific needs.
(i.e. you use Kilosort2 instead of Kilosort4, or you have an specific channel map file).

In general the files are:

- scripts with parameters for different processing steps in the pipeline.
  (e.g. kilosortConfig.m, parameters_default.py)
  These would need to be modified, in most cases, to work properly on your dataset.

- channel map files for specific silicon probes or combinations of probe-setup
  (e.g. chanMapE32-S2.mat and chanMapE32-S2_DeutSN11.mat)

- Information files like this one or HowToData.m