function plot_genstuff(sp)

%% basic quantification of spiking plot
depthBins = 0:50:750;
ampBins = 0:20:min(max(spikeAmps),500);
recordingDur = sp.st(end);

[pdfs, cdfs] = computeWFampsOverDepth(ap.amps, sp.depths, ampBins, depthBins, recordingDur);

plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);

%% Computing some useful details about spikes/neurons (like depths)
[~, ~, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

end