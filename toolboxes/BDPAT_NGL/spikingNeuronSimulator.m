%% spike train for different firing rates
%%
% Use this script to simulate spike trains. Spike events are based on a
% Poisson distribution (user input determines lambda). Automatically plots
% a spike dot raster and PSTH. Automatically calculates a ttest on an
% chosen bin to estimate power.
%
% Dependencies (requires the following custom functions):
%   calcFireRate
%   generateSpikeTrain
%   nanMeanSterrHistogram
%   nanste
%   plotRaster
%   plotPSTH


% VERSION HISTORY:
% Author:         Lukas Hahn
% Version:        1.0.0
% Last Change:    25.01.2024
%
% 25.01.2024, Lukas: v1.0.0 release version

%%
maxReps = 300; %number of repetitions to test power

producePlot = true; %plot of spike raster and PSTH, if true will produce
% one plot per outer cycle, if false will produce no plot

fullInterval = 2000; %length of the full interval for simulation, should be
% a multiple of binSize
binSize = 200; %in ms

%trials per condition (1 set of conditions per row)
trls = [20 20; 30 30; 20 30; 30 20];

FR = [3 6]; %defined firing rates (1 set of conditions per row)

%Note: When both variables (trls and FR) have the same number of rows they
% will be tested as corresponding row-pairs. If one of the two has only 1
% row and the other more than 1 then the higher row variable will cycle
% through rows and keep the lower variable stable. If both have more than 1
% row, but not the same amount, the code will not execute!

tstBin = 4; %bin for fire rate testing

% stimulation options, requires input per FR
stimCond = [0 0]; %true or false index for stimulation of respective FR

%start of stimulation (bin no.), use 0 to have no effect
stimBinS = [0 0];
%end of stimulation (bin no.), use 0 to have no effect
stimBinE = [0 0];

%factor of stimulation, >1 (e.g. 2 doubles the firing rate),
% <1 (e.g. 0.5 halves firing rate)
stimulationFactor = [1 1];
%interval start within a bin where stimulation
% (increase or decrease of spiking) occurs
stimMin = [0 0];
%interval end within a bin where stimulation
% (increase or decrease of spiking) occurs
stimMax = [200 200];

testPower = nan(size(FR,1),1);
rasterLineSize = 3;
cols = [0 0 0; 1 0 0; 0 0.8 0; 0 0 1; 0.5 0.8 0; 0.2 0.5 0.8];
upperY = max(FR)*max(stimulationFactor)+5; %y-limit for plotting in PSTH
%% check input
if size(FR,1)>1 && size(trls,1)>1 && size(FR,1)~=size(trls,1)
    error('nuber of rows of FR and trls are incompatible')
end

for i=1:size(stimCond,2)
    if stimCond(i)~=0
        if any([stimBinS(i)==0 stimBinE(i)==0])
            error(['Invalid stimulation bin start and/or end, check ' ...
                'stimBinS and stimBinE'])
        end
        if stimMax(i)<=0 || stimMax(i)>binSize
            stimMax(i) = binSize;
            warning(['Invalid maximum for stimMax. stimMax was set to ' ...
                num2str(binSize)])
        end
        if stimMin(i)<0 || stimMin(i)>=stimMax(i)
            stimMin(i) = 0;
            warning('Invalid value for stimMin. stimMin was set to 0')
        end
    end
end
%%
%cycle through this variable to test how power changes, automatically
%selects larger of the two variables, defaults to FR
% (e.g., stable trial numbers at changing FR or the other way around)
if size(trls,1)>size(FR,1)
    cycleVar = trls;
elseif size(trls,1)<size(FR,1)
    cycleVar = FR;
else
    cycleVar = FR;
end
for ex=1:size(cycleVar,1)
    %check variables and adjust as necessary
    if size(trls,1)<size(FR,1)
        trls = ones(size(FR,1),2).*trls;
    elseif size(trls,1)>size(FR,1)
        FR = ones(size(trls,1),2).*FR;
    end
    %%
    pd = cell(1,length(FR));
    pdStim = cell(1,length(stimCond));
    for FRInd=1:size(FR,2)
        %(e.g., for 10 spks per second expect 2 spikes per 200 ms)
        pd{1,FRInd} = makedist('Poisson','lambda',...
            FR(ex,FRInd)/(1000/binSize));
        for s=1:length(stimCond)
            if stimCond(s)==1
                %adjusted spiking distribution
                pdStim{1,s} = makedist('Poisson','lambda',...
                    (FR(ex,FRInd)/(1000/binSize))*stimulationFactor(s));
            end
        end
    end

    h = nan(maxReps,1);
    p = nan(maxReps,1);
    for reps=1:maxReps
        fireRate = cell(2,1);
        trialCounter = 1;
        f = 0;
        for FRInd=1:size(FR,2)
            f = f+1;
            %e.g. 2 seconds in bins of binSize, for 100 trials
            spksPerSec = random(pd{1,FRInd},...
                fullInterval/binSize,trls(ex,FRInd));

            %stimulation/inhibition
            if stimCond(FRInd)==1 && stimBinS(FRInd)~=0 && ...
                    stimBinE(FRInd)~=0
                spksPerSec(stimBinS(FRInd):stimBinE(FRInd),:) = ...
                    random(pdStim{1,FRInd},...
                    stimBinE(FRInd)-stimBinS(FRInd)+1,trls(ex,FRInd));
            end

            spkTrain = generateSpikeTrain(spksPerSec,binSize,...
                'conditionIndex',FRInd,...
                'stimCondition',stimCond,...
                'stimFactor',stimulationFactor,...
                'stimStart',stimBinS,...
                'stimEnd',stimBinE,...
                'stimMin',stimMin,...
                'stimMax',stimMax);

            if producePlot && reps==1
                figure(ex)
                subplot(2,1,1)
                trialCounter = plotRaster(spkTrain,trialCounter,...
                    cols(f,:),rasterLineSize,'lines');
                xlim([0 fullInterval])
                set(gca,'YTick',0:20:trialCounter,'YTickLabel',...
                    0:10:trialCounter/2)
                xlabel('time (ms)')
                ylabel('repetitions (trials)')
                box off
                set(gcf,'color','white')
                set(gca,'LineWidth',3,'XColor',[0 0 0],'YColor',[0 0 0],...
                    'TickDir','out','FontSize',15)
                subplot(2,1,2)
                plotPSTH(spkTrain,binSize,binSize,[0 fullInterval],1000,...
                    'plotCol',cols(f,:),'meanline','-','smoothPlot',0);
                xlim([0.5 fullInterval/binSize+0.5])%center bins in  middle
                xlabel('bin no. (middle of bin)')
                ylabel('smoothed firing rate (spks/s)')

                for s=1:length(stimCond)
                    if stimCond(s)==1 %draw lines of stimulation interval
                        line([stimBinS(s) stimBinS(s)],[0 upperY],...
                            'LineWidth',3,'Color',[0.5 0.5 0.5])
                        line([stimBinE(s) stimBinE(s)],[0 upperY],...
                            'LineWidth',3,'Color',[0.5 0.5 0.5])
                    end
                end
                box off
                set(gca,'LineWidth',3,'XColor',[0 0 0],'YColor',[0 0 0],...
                    'TickDir','out','FontSize',15)
                set(gcf,'color','white','Position',[200 70 700 800])
            end

            fireRate(f,1) = calcFireRate(spkTrain,binSize,binSize,...
                [0 fullInterval],1000);
        end
        [h(reps,1), p(reps,1)] = ttest2(fireRate{1,1}(:,tstBin),...
            fireRate{2,1}(:,tstBin));
    end
    testPower(ex) = sum(h,'omitnan')/reps; %achieved power
end
%% power simulation results
figure(ex+1)
plot(testPower,'Marker','o','LineWidth',3)
hold on
xlim([0.9 length(testPower)+0.1])
ylim([0 1])
xlabel('condition')
ylabel('simulated estimate of power')
box off
set(gca,'LineWidth',3,'XColor',[0 0 0],'YColor',[0 0 0],'TickDir','out',...
    'FontSize',15)
set(gcf,'color','white','Position',[40 70 700 800])