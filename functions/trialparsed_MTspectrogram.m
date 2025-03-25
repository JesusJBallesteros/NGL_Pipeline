function [TFR] = trialparsed_MTspectrogram(FT_data, conditions, param, opt)
%
%
%

%% Defaults
param.multi = false;
t_baseline = 4; % in sec
nblocks = 1;

% Figure size
if strcmpi('adaptive', param.size)
    param.screen.size   = get(0, 'ScreenSize');  
    param.screen.width  = param.screen.size(3);
    param.screen.height = param.screen.size(4);
else
    param.screen.width  = param.size(1); 
    param.screen.height = param.size(2);
end

%% Wavelet default configuration
% Absolute calculation
cfg = [];
 cfg.method     = 'wavelet';
 cfg.output     = 'pow';      % power output
 cfg.foi        = 1:1:50;  % frequencies from 1 to 80 Hz
 cfg.width      = 7;          % wavelet width (can be tuned based on your analysis)
 cfg.keeptrials = 'yes';    
 cfg.toi        = -4:.1:8;  % time vector from -4 s to 8 s (adjust time resolution as needed)
 cfg.pad        = 'nextpow2';

% Relative calculation
cfgb              = [];
 cfgb.baseline     = [-t_baseline 0];
 cfgb.baselinetype = 'zscore';
 cfgb.stimType     = 0; % adapt to whatever (EG. 1=FS, 0=NS) 

%% Check existence of multiple levels
if max(unique(conditions.block)) > 1
    nblocks = max(unique(conditions.block));
    param.multi = true; 
end

%% 1 Proceed
for b = 1:nblocks % per block
    % define valid trials
    cfgt = [];         
     cfgt.trials = find(conditions.block == b & conditions.stimulus == cfgb.stimType & ~isnan(FT_data.cfg.trl(:,3)));
        blockData = ft_redefinetrial(cfgt, FT_data);
    
    % Skip if there are no trials in this block
    if isempty(blockData.trial), continue; end

    % Perform Time-Frequency Representation (TFR) Analyses 
    TFR.abs{b} = ft_freqanalysis(cfg, blockData);

    % Relativize to baseline.
    TFR.rel{b} = ft_freqbaseline(cfgb, TFR.abs{b});

    % % Average over trials for each case
    % cfg.keeptrials = 'no';    
    %     TFR.absAV{b} = ft_freqanalysis(cfg, blockData);
    %     TFR.relAV{b} = ft_freqbaseline(cfgt, TFR.absAV{b});
end
clear blockData b

%% 2 Plot dB scaled TFRs.
% Plot Configuration
cfgp           = [];
 cfgp.figure    = 'gca';
 cfgp.colormap  = hot;
 cfgp.xlim      = 'maxmin';
 cfgp.ylim      = [0 max(cfg.foi)];
 cfgp.clim      = 'maxmin';
 cfgp.colorbartext = cfgb.baselinetype;
 cfgp.interactive = 'no';
 cfgp.fontsize  = 12;
figtitl = ['TFR_StimOn2_Blocks_NS_',cfgb.baselinetype,'_1-50'];

% Plot itself
fig = figure('Visible', param.visible, 'Position', [0 0 700 400]);
set(fig, 'Position', [0, 0, round(param.screen.height), round(param.screen.width)]); % Set fig size as screen
for b = 1:nblocks
 cfgp.title = sprintf('Block %d',b); 
    subplot(nblocks/2, 2, b)
    try ft_singleplotTFR(cfgp, TFR.rel{b}); hold on
    catch, end
    xlim([-3.1 6.9]);
    clim([0 3]);
    xticks(TFR.rel{b}.time(1):2:TFR.rel{b}.time(end));
    xticklabels(round(TFR.rel{b}.time(1):2:TFR.rel{b}.time(end)));
    xlabel('Time (s)'); ylabel('Frequency (Hz)');    
    hold off
    xline(0, '--w', 'LineWidth', 2);
    xline(-2, ':w', 'LineWidth', 2);
end
sgtitle(fig,['Session Type: ', conditions.type, '. AllCh']); 

%% Save the figure as .png
if ~isfolder(fullfile(opt.analysis,'plots', 'TFR')), mkdir(fullfile(opt.analysis,'plots', 'TFR')); end
exportgraphics(fig, fullfile(opt.analysis,'plots', 'TFR', [figtitl '.png']),'Resolution', param.Resolution);
close all hidden

%% Save the data
save(fullfile(opt.analysis,[figtitl '.mat']), 'TFR','-mat');
end