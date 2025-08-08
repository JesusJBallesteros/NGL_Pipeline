# CAD*opti*

This repository contains MATLAB code for a machine learning algorithm for the detection of functional cell assemblies at their optimal time scale (CAD*opti*) from simultaneously recorded spike trains.
##

The algorithm extends from the cell assembly detection algorithm (CAD) presented in
[*Russo, Durstewitz, Cell assemblies at multiple time scales with arbitrary lag constellations. Elife. 2017 6*](https://elifesciences.org/articles/19428) and available at [here](https://github.com/DurstewitzLab/Cell-Assembly-Detection), and returns, for each set of coordinated units, only the configuration (in terms of activity pattern and temporal resolution) of maximal temporal coordination (lowest p-value). 

Similarly to CAD, CAD*opti* detects functional cell assemblies of arbitrary size through an agglomerative scheme, which recursively tests and adds new units to pre-tested assemblies. The new agglomeration algorithm proceeds as follows: Each unit pair is tested at all time scales specified in `BinSizes`. For each significant unit-pair, the time scale at which the pairwise interaction has the lowest p-value is selected. This step fixes the assembly characteristic time scale, any further agglomeration step will then be performed at that temporal resolution. Significance of all tests performed at each agglomeration step is corrected for multiple comparisons with the Bonferroni-Holm method, only significant activity patterns pass to the next step. At the end of each step, if the complete set of units composing an assembly is present also in a different assembly of same size, only the assembly with lowest p-value is kept for further agglomeration. This latter case applies when different unit-pairs, part of an overarching assembly, converge to the same assembly upon agglomeration of other assembly-units.

CAD/CAD*opti* comparison: CAD*opti* returns a sub-selection of the assemblies detected by CAD, automatically selecting for the time scale at which assembly units are best coordinated. For this reason, CAD*opti* is faster than CAD. Note that if interested in the specific question of whether the same set of units is coordinated at two time scales, CAD*opti* has to be run separately on the two time scales of interest. 


Code tested for MATLAB R2018a. In case of questions or comments, please contact eleonora.russo@zi-mannheim.de. 
Please credit the source and cite *Russo and Durstewitz 2017* and *Oettl et al. 2020* when using the code in any form of publication.


*Copyright*: 2020 Russo Eleonora, Dept. of Theoretical Neuroscience, Central Institute of Mental Health, Medical Faculty Mannheim, Heidelberg University.

### Usage

For a short tutorial on how to use CAD*opti* and visualize the result, please, run `Tutorial_CADopti.m`

