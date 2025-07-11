# BasicDataProcessingAndAnalysisTools



## What this is

This directory contains the lab's shared data processing and analysis functions:

- Matlab scripts
- Matlab functions
- Documentation on how to do analyses (reference links to places where more details can be obtained from (**also check the [Wiki](https://gitlab.ruhr-uni-bochum.de/ngl/basicdataprocessingandanalysistools/-/wikis/home)!**)

## Add your files!

Lab members are obliged to add new tools to this directory. This means if anyone has developed something new (e.g., new analysis code) it has to be added in a way that makes it accessible to everyone in the lab.


# How to treat this directory

Please add a note below into the table of contents whenever you add something.

***

## Table of contents

calcFireRate - takes spike times and calculates binned firing rate (spks/s), called by plotRaster

calculateMutualInformation - calculates an information theory metric about mutual dependence of two variables, can be used for example in classifier analysis

calculatePEV - takes results of an ANOVA table and calculates effect size

generateSpikeTrain - creates a spike train for an arbitrary amount of bins and trials

nanMeanSterrHistogram - takes binned input (e.g., firing rate per bin) and plots mean line and SEM as shaded area around the mean, 
called by plotPSTH

nanste - basic helper function, calculates standard error of the mean, called by nanMeanSterrHistogram

pcaTest - script running a PCA on neuronal data, example data provided in additionalInformation

plotPSTH - plots a PSTH

plotRaster - plots a dot raster diagram

spikingNeuronSimulator - simulates neuronal spiking, can be used to estimate power of statistical analysis of fire rate comparisons

## Dependencies

- nanMeanSterrHistogram -> nanste
- plotPSTH -> calcFireRate
- plotPSTH -> nanMeanSterrHistogram
- spikingNeuronSimulator -> calcFireRate
- spikingNeuronSimulator -> generateSpikeTrain
- spikingNeuronSimulator -> nanMeanSterrHistogram
- spikingNeuronSimulator -> plotRaster
- spikingNeuronSimulator -> plotPSTH
