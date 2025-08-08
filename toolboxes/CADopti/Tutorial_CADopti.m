%% %%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%
clear all

load('Data.mat');  % Data.mat contains a set of simulated parallel spike trains 
% containing 5 cell assemblies of 5 units each embedded in a set of 50
% units. Assembly time scales and activity patterns are described in Russo
% and Durstewitz 2017.

% spM := matrix containing the parallelly recorded spike time series
nneu=size(spM,1);  % nneu is the number of recorded units

%% %%%%%%%%%%%%%%%%%%%%%%%% DETECTION %%%%%%%%%%%%%%%%%%%%%%%%
% Range of temporal resolution to test
BinSizes=[0.015 0.025 0.04 0.06 0.085 0.15 0.25 0.4 0.6 0.85 1.5];
% Range of lags to test
MaxLags=[10 10 10 10 10 10 10 10 10 10 10];

% assembly detection:
[As_across_bins,As_across_bins_index,assembly]=CADopti(spM,MaxLags,BinSizes);


% other detection options:
% to specity the reference lag 
% ref_lag=2; % more strict correction on the non stationarity, you find less assemblies
% ref_lag=4;% less strict correction on the non stationarity, you find more assemblies
% [As_across_bins,As_across_bins_index,assembly]=CADopti(spM,MaxLags,BinSizes, ref_lag);

% to detect only unit pairs
% [As_across_bins,As_across_bins_index,assembly]=CADopti(spM,MaxLags,BinSizes, [],[],[],2);   

% to change the significance level (this algorithm, with respect to CAD, it has a harsher multiple 
% comparison correction on the significance level, you might therefore also relax a bit the 0.05 alpha level)
% [As_across_bins,As_across_bins_index,assembly]=CADopti(spM,MaxLags,BinSizes,[],0.1);   

%% %%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%
% Visualization of detected assemblies in the assembly assignment matrix (Russo and Durstewitz, 2017)
% see comments inside "assembly_assignment_matrix.m" for details on how to fix the function paramethers

display='raw';
% display='clustered';
clf
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins,nneu,BinSizes,display);

%% %%%%%%%%%%%%%%%%%%%%%%%% ASSEMBLY ACTIVATION %%%%%%%%%%%%%%%%%%%%%%%%
% visualization of the activity of the assembly in time
% see comments inside "assembly_activity_function.m" for details on how to fix the function parameters

lagChoice='duration';
act_count = 'full';
[assembly_activity]=assembly_activity_function(As_across_bins,assembly,spM,BinSizes,lagChoice,act_count);

clf
for i=1:length(assembly_activity)
    subplot(size(assembly_activity,1),1,i)
    plot(assembly_activity{i}(:,1),assembly_activity{i}(:,2));
    hold on
end

