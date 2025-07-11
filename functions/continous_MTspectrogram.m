function [TFR] = continous_MTspectrogram(FT_data, conditions, param, opt)
%{Intended to calculate and plot the Multitaper Spectrogram of continuous
%  data from anesthesia sessions.
% Requires a FieldTrip ready file, normally generated from an ETALO- 
% preprocessed continuous data. It's format would be:
% 'ANIMAL_Experiment_YYMMDD_FT.mat'
%
%
%
%
%}

% MTM
cfg              = [];
 cfg.channel     = FT_data.label(1);
 cfg.method      = 'mtmconvol'; % uses MT convolution 
 cfg.output      = 'pow';       % Outputs power and cross-spectra
 cfg.taper       = 'dpss';      % Slepian sequence
 cfg.foi         = 1:1:40;     % freqs of interest
 cfg.t_ftimwin   = ones(length(cfg.foi),1).*1; % 5s windows
 cfg.toi         = '50%';       % stimate every .5s 
 cfg.tapsmofrq   = 2;          % smooth over +/- 2 Hz
 cfg.keeptrials  = 'no';       % Do (not) keep indiv. trial data
 cfg.polyremoval = 0;
 cfg.pad         = 'nextpow2';

%% Figure size
if strcmpi('adaptive', param.size)
    param.screen.size   = get(0, 'ScreenSize');  
    param.screen.width  = param.screen.size(3);
    param.screen.height = param.screen.size(4);
else
    param.screen.width  = param.size(1); 
    param.screen.height = param.size(2);
end

% t_baseline = varin.baseline(1); % in sec. 36 minutes

% for a=1:length(varin.arrays) % per array
    %% 1 Perform Time-Frequency Representation (TFR) Analyses 
    % abs = ft_freqanalysis(cfg, FT_data.(varin.arrays{a}));
    abs = ft_freqanalysis(cfg, FT_data);
    clear FT_data

    % 2 Transform to dB.
    cfg2              = [];
     cfg2.baselinetype = 'db';
     cfg2.baseline     = 'no'; 
     % cfg2.baseline     = [0 t_baseline]; % Baseline from 0 to approx. 36 min
        TFR.dB = ft_freqbaseline(cfg2, abs);
        % TFR.dB.(varin.arrays{a}) = ft_freqbaseline(cfg2, abs);
    clear abs
    
    % 4 Plot dB scaled TFRs.
    cfg2           = [];
     cfg2.colormap  = hot;
     cfg2.xlim      = 'maxmin';
     cfg2.ylim      = [0 max(cfg.foi)];
     cfg2.zlim      = [0 500];
     cfg2.colorbartext = 'Norm. Pow. (dB)';
     cfg2.interactive = 'no';
     % cfg2.colorbartext = 'Norm. Pow. v. Baseline (dB)';
     cfg2.fontsize  = 12;
     cfg2.title = 'allCh';    
     % cfg2.title = varin.arrays{a};    
    
    fig = figure('Visible', param.visible, 'Position', [0 0 700 400]);
    set(fig, 'Position', [0, 0, round(param.screen.width), round(param.screen.height)]); % Set fig size as screen
        ft_singleplotTFR(cfg2, TFR.dB); hold on
        ylim([0 40]);
        xticks(TFR.dB.time(1):600:TFR.dB.time(end));
        xticklabels(TFR.dB.time(1)/60:10:TFR.dB.time(end)/60);
        xlabel('Time (min)'); ylabel('Frequency (Hz)');
        
        hold off
        % xline(t_baseline, '--k', 'LineWidth', 2);

    % figure;
    % ft_singleplotTFR(cfg2, TFR.dB.(varin.arrays{a})); hold on
    %   xticks(FT_data.(varin.arrays{a}).time{1}(1):1800:FT_data.(varin.arrays{a}).time{1}(end));
    %   xticklabels(FT_data.(varin.arrays{a}).time{1}(1)/60:30:FT_data.(varin.arrays{a}).time{1}(end)/60);
    %   xlabel('Time (min)'); ylabel('Frequency (Hz)');
    %   sdf('Spectrogram');
    %   xline(t_baseline, '--k', 'LineWidth', 2);
    %   hold off;

    % Save the figure as .png
    if ~isfolder(fullfile(opt.analysis,'plots', 'TFR')), mkdir(fullfile(opt.analysis,'plots', 'TFR')); end
    exportgraphics(fig, fullfile(opt.analysis,'plots', 'TFR', ...
                        'Cont_allCh.png'), ...
                        'Resolution', param.Resolution);
    close all hidden
    
    % % 3 Z-score to baseline to measure change.
    % cfg2              = [];
    %  cfg2.baselinetype = 'zscore';
    %  cfg2.baseline     = [0 t_baseline]; % Baseline from 0 to approx. 36 min
    % 
    % TFR.Zsc.(varin.arrays{a}) = ft_freqbaseline(cfg2, TFR.dB.(varin.arrays{a}));
    %
    % % 5 Plot ZScored TFRs.
    % cfg2           = [];
    %  cfg2.zlim      = [-3 3]; % Fix the zscore scale
    %  cfg2.colormap  = varin.map.redblue; %
    %  cfg2.colorbartext = 'Pow. change v. baseline (Z)'; %'ZScore'; %'dB change';
    %  cfg2.fontsize  = 12;
    % 
    %  cfg2.xlim  = [FT_data.(varin.arrays{a}).time{1}(1) FT_data.(varin.arrays{a}).time{1}(end)];
    %  cfg2.title = varin.arrays{a};    
    % 
    % figure;
    % ft_singleplotTFR(cfg2, TFR.Zsc.(varin.arrays{a})); hold on
    %   xticks(FT_data.(varin.arrays{a}).time{1}(1):1800:FT_data.(varin.arrays{a}).time{1}(end));
    %   xticklabels(FT_data.(varin.arrays{a}).time{1}(1)/60:30:FT_data.(varin.arrays{a}).time{1}(end)/60);
    %   xlabel('Time (min)'); ylabel('Frequency (Hz)');
    %   sdf('Spectrogram');
    %   xline(t_baseline, '--k', 'LineWidth', 2);
    %   hold off;
    % 
    % % Save the figure as .png
    % savname = strcat(varin.saveplot,'\TFR_',cfg2.title,'_Zsc.png');
    % saveas(gcf,savname);
    % close gcf
% end
       
end