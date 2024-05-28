%Example script PCA
%
% This script runs a PCA on (z-scored) firing rate data (spikes per bin), 
% and plots the result for the first three components. Time bins indicate
% temporal progression across trial. Based on data recorded by Erica
% Fongaro ('chgDtctPCue'). Processing of data based on analysis in Ott &
% Nieder, 2016 (https://doi.org/10.1093/cercor/bhw244).
%
%INPUTS
%  * 'path of data'   : String to location where data is stored
%  * 'anl'            : subject index (1 or 2)
%  * 'rgn'            : region index (1)
%  * 'zFR1' - 'zFR3'  : matrix of z-scored firing rates of neurons 
%                      (trials x timeBins x neurons) for trial types 1 - 3
%  * 'selRange'       : vector defining time bins for analysis
%  * 'locCols'        : color for different trials types
%
%OUTPUTS
%   * One Figure with two subplots (PCA space across time points and
%   Euclidean distance between time points of different trial types)

% VERSION HISTORY:
% Author:         Lukas Hahn
% Version:        1.0.0
% Last Change:    19.12.2023
%
% 19.12.2023, Lukas: v1.0.0 release version
%% INPUT
load('C:\Users\hahnlu1m\Desktop\storedValuesForPCA')
anl = 1;
rgn = 1;
zFR1 = storeVals.zFR1{anl,rgn};
zFR2 = storeVals.zFR2{anl,rgn};
zFR3 = storeVals.zFR3{anl,rgn};

selRange = 1:1:151; %selected time bins

locCols =  [0.9290 0.6940 0.1250; 0 0.3470 0.6410; 0.8 0 0];
locCols1 = [0.929 0.694 0.125].*linspace(1,0.1,151)'; %color gradient 1
locCols2 = [0 0.347 0.641].*linspace(1,0.1,151)'; %color gradient 2
locCols3 = [0.8 0 0].*linspace(1,0.1,151)'; %color gradient 3

% events during trial (bin numbers indexing zFR)
pCON = 11;
pCOFF = 21;
smpON = 46;
smpOFF = 66;
comON = 121;
chcON = 141;
%% setup and run PCA
zHlp1 = squeeze(mean(zFR1(:,selRange,:),1,'omitnan'));
zHlp2 = squeeze(mean(zFR2(:,selRange,:),1,'omitnan'));
zHlp3 = squeeze(mean(zFR3(:,selRange,:),1,'omitnan'));

PHlpComb = [zHlp1; zHlp2; zHlp3];

clear score
[coeff,score,latent,tsquared,explained] = pca(PHlpComb);

subplot(1,2,1)
for t=1:size(zHlp1,1)
    if t==1
        markerT = '*'; % maker of trial start
        sz = 100;
    elseif t==size(zHlp1,1)
        markerT = 'x'; % marker of trial end
        sz = 100;
    else
        markerT = 'o'; % marker of bins between start and end
        sz = 20;
    end
    %plot PC scores (first three dimensions) per time point
    %color gradients go from light (start) to dark (end) in trial
    scatter3(score(t,1),score(t,2),score(t,3),...
        'MarkerFaceColor',locCols1(t,:),'MarkerEdgeColor',locCols1(t,:),...
        'Marker',markerT,'SizeData',sz);
    hold on
    scatter3(score(t+size(zHlp1,1),1),score(t+size(zHlp1,1),2),...
        score(t+size(zHlp1,1),3),'MarkerFaceColor',locCols2(t,:),...
        'MarkerEdgeColor',locCols2(t,:),...
        'Marker',markerT,'SizeData',sz);
    scatter3(score(t+2*size(zHlp1,1),1),score(t+2*size(zHlp1,1),2),...
        score(t+2*size(zHlp1,1),3),'MarkerFaceColor',locCols3(t,:),...
        'MarkerEdgeColor',locCols3(t,:),...
        'Marker',markerT,'SizeData',sz);
end

t = 1:size(zHlp1,1);
plot3(score(t,1),score(t,2),score(t,3),'Color',locCols(1,:),'LineWidth',3)
hold on
plot3(score(t+size(zHlp1,1),1),score(t+size(zHlp1,1),2),...
    score(t+size(zHlp1,1),3),'Color',locCols(2,:),'LineWidth',3)
plot3(score(t+2*size(zHlp1,1),1),score(t+2*size(zHlp1,1),2),...
    score(t+2*size(zHlp1,1),3),'Color',locCols(3,:),'LineWidth',3)

axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
box off
set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',3,'TickDir','out')
legend({'cue top' 'cue mid' 'cue bot'})
%% euclidean distance between points in different trial types
eucDist = nan(size(zHlp1,1),3);
for t=1:size(zHlp1,1)
    eucDist(t,1:3) = pdist([score(t,1) score(t,2) score(t,3); ...
        score(t+size(zHlp1,1),1) score(t+size(zHlp1,1),2) ...
        score(t+size(zHlp1,1),3); score(t+2*size(zHlp1,1),1) ...
        score(t+2*size(zHlp1,1),2) score(t+2*size(zHlp1,1),3)]);
end
subplot(1,2,2)
plot(eucDist(:,1),'Color',mean([locCols(1,:); locCols(2,:)],1),...
    'LineWidth',3) %top to mid
hold on
plot(eucDist(:,2),'Color',mean([locCols(1,:); locCols(3,:)],1),...
    'LineWidth',3) %top to bot
plot(eucDist(:,3),'Color',mean([locCols(2,:); locCols(3,:)],1),...
    'LineWidth',3) %mid to bot

legend({'\Delta top mid' '\Delta top bot' '\Delta mid bot'})
xlim([1 size(zHlp1,1)])
xlabel('time (bin)')
ylabel('Euclidean dist. between trial types (\Delta)')
box off
set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',3,'TickDir','out')

f = gcf;
f.Children(4).Position = [-0.015 0.1 0.8 0.75];
f.Children(3).Position = [0.1 0.9 0.1 0.05];
f.Children(2).Position = [0.77 0.1 0.2 0.8];
f.Position = [10 10 1000 900];
f.Color = 'white';