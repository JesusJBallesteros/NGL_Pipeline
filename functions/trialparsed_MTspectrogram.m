function [TFR, cfg] = trialparsed_MTspectrogram(FT_data, conditions, param, opt)
% This function takes trial-parsed FieldTrip formatted data, the conditions
% file and parameters and options to calculate the Time-frequency representation 
% (TFR) with the method of choice. For wide-ranges 'superlets' is
% recommended. By default, it will keep all trials information.
% 
% So far, it has a block-by-block approach, with early and late sub-block
% iterations. For each of these, it will perform for each type of stimuli
% present (e.g., Familiar vs Novel = 2x). Before each call to 'ft_freqanalysis'
% is preceded by an artifact rejection function. Once TFR absolute is 
% calculated, the value of one condition will be substracted from the
% other stimuli type trials. 
%
% Baseline calculations would be performed at the time of plotting. 
% 
% It will save the single session data and output it to save to a general
% file, outside the function, if requested.
%
% Plotting options can be added at the end of the function.
%
% Jesus. 14.04.2025

%% Defaults
if ~isfield(param,'multi'),     param.multi     = false;        end
if ~isfield(param,'testname'),  param.testname  = 'trial_TFR_'; end
if ~isfield(opt,'blocks'),      opt.blocks      = 'all';        end
toload = 0;

%% Check existence of multiple levels
if isstring(opt.blocks)
    opt.blocks = max(unique(conditions.block));
    param.multi = true; 
elseif ~isstring(opt.blocks)
    if opt.blocks > 1; param.multi = true; end
end

%% Check for processed data
fTFR = ls(fullfile(opt.analysis, [param.testname '*'])); % 'ASL_*
if ~isempty(fTFR)
    % If existing load
    if isfile(fullfile(opt.analysis, fullfile(fTFR)))
        toload = 1;
        disp('Loading existing TFR data ...')
        load(fullfile(opt.analysis, fullfile(fTFR)));
    end
end

%% TFR Absolute calculation
if ~toload
    global ft_default;
    ft_default.notification.warning = [];
    cfg = [];
    for f = 1:(size(opt.freqInterest,1))
        cfg{f}.method = opt.TFRmethod;
        cfg{f}.output = 'pow'; % Outputs power and cross-spectra
        cfg{f}.pad    = 'nextpow2';
        cfg{f}.keeptrials  = 'yes'; % Do (not) keep indiv. trial data
        cfg{f}.foi    = opt.freqInterest{f}; % frequencies
        
        if strcmpi(opt.TFRmethod, 'wavelet')
         % Wavelet default configuration
         cfg{f}.width      = 7; % wavelet width (can be tuned based on your analysis)
         cfg{f}.toi        = -4:.1:8; % time vector from -4 s to 8 s (adjust time resolution as needed)
         
        elseif strcmpi(opt.TFRmethod, 'mtmconvol')
         % Multitaper default configuration
         cfg{f}.taper       = 'hanning';      % sequence
         cfg{f}.t_ftimwin   = 3./cfg{f}.foi;  % 3 cycles per time window
         cfg{f}.toi         = linspace(-4,8,length(cfg{f}.foi));   % stimate
         cfg{f}.tapsmofrq   = 4;          % smooth over +/- Hz
    
        elseif strcmpi(opt.TFRmethod, 'superlet')
         cfg{f}.toi     = -4:opt.timeResol:8; % time vector from -4 s to 8 s (adjust time resolution as needed)
         cfg{f}.width   = opt.width{f};
         cfg{f}.combine = opt.combine;
         
         AF{f} = calculate_superlet_order(cfg{f}.foi, opt.superletOrder{f});
         cfg{f}.order   = AF{f};
        end
    
        %% 1 Proceed
        for b = 1:opt.blocks
            % For both stimuli NS/FS
            for st = 1:2 % st==1 -> NS, st==2 -> FS (bc 0=NS, 1=FS)
                % get valid trials and halven into early/last
                seltrials = find(conditions.correct == 1 & conditions.block == b & conditions.stimulus == (st-1) & ~isnan(FT_data.cfg.trl(:,3)));
                midblock = floor(size(seltrials,1)/2);
                blocktrials = {seltrials(1:midblock); ...
                               seltrials(midblock+1:end)};
    
                % Do first/last half of the block
                for i = 1:2
                    cfgt = [];
                     cfgt.trials = blocktrials{i};
                    blockData = ft_redefinetrial(cfgt, FT_data);
                
                    % Skip if there are no trials in this block
                    if isempty(blockData.trial), continue; end
                    
                    % Artifact detection and rejection
                    if opt.artifdet
                        blockData.cfg.trl = FT_data.cfg.trl(blockData.cfg.trials,:);
                        blockData = artifact_detRej_lfp(blockData, opt);
                    end

                    % Perform Time-Frequency Representation (TFR) Analyses 
                    TFR{f,b}.abs{i,st} = ft_freqanalysis(cfg{f}, blockData);
                end
            end
        
            % Again for Early/late, substract Abs pow, FS from NS (NS-FS)
            for i = 1:2
                % copy Abs NS data
                TFR{f,b}.rem{i} = TFR{f,b}.abs{i,1};
    
                % calculate Abs FS mean across trials
                meanFS = mean(TFR{f,b}.abs{i,2}.powspctrm, 1, 'omitnan');
    
                % Substract the Abs FS mean power from each Abs NS power trial
                for t = 1:size(TFR{f,b}.rem{i}.powspctrm, 1)
                    TFR{f,b}.rem{i}.powspctrm(t,:,:,:) = TFR{f,b}.rem{i}.powspctrm(t,:,:,:) - meanFS;
                end
            
                % empty 'cumtapcnt'
                TFR{f,b}.rem{i}.cumtapcnt = [];
            end
        end
    end
    clear blockData
    
    % Save TFR data
    filtitl = [param.testname, opt.SavFileName, '_',opt.TFRmethod];
    save(fullfile(opt.analysis,[filtitl '.mat']), 'TFR', 'cfg','-mat');
end

% % 2 Plot Baselined TFRs.
% plot_superletsTFR_extintion(TFR, param, opt)

if opt.trialbytrial
    % 3 Plot trial-by-trial TFRs.
    plot_superletsTFR_extintion_tbt(TFR, param, opt)
end

if opt.chbych
    % 4 Plot trial-by-trial TFRs.
    plot_superletsTFR_extintion_chbych(TFR, param, opt)
end

end