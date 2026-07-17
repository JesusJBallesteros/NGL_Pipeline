function [TFR, cfg] = computeTrialparsedTFR_ASL(FT_data, condition, param, opt, alignName)
% computeTrialparsedTFR_ASL  Social-learning-specific trial-parsed TFR.
%
% RELOCATION (26.06.2026, LFP Pass 2):
%   Moved from functions/analysis/trialparsed_MTspectrogram.m to this
%   projects/socialLearning/ folder because the block-early/late x NS/FS
%   x NS-FS-subtraction logic is specific to the ASL (extintion) task
%   and non-ASL projects would crash on condition.stimulus / .block.
%   The generic (project-agnostic) core lives at
%   functions/analysis/computeTrialparsedTFR.m; the top-level
%   trialparsed_MTspectrogram.m is now a thin shim that dispatches on
%   opt.proj_socialLearning.
%
% This function takes trial-parsed FieldTrip formatted data, the condition
% file and parameters and options to calculate the Time-frequency representation
% (TFR) with the method of choice. For wide-ranges 'superlets' is
% recommended. By default, it will keep all trials information.
%
% So far, it has a block-by-block approach, with early and late sub-block
% iterations. For each of these, it will perform for each type of stimuli
% present (e.g., Familiar vs Novel = 2x). This whole block-level layout is
% project-specific (social-learning task); Pass 2 of the LFP refactor
% splits this into a generic core plus a projects/socialLearning/ wrapper.
% Once TFR absolute is calculated, the value of one condition will be
% substracted from the other stimuli type trials.
%
% Baseline calculations would be performed at the time of plotting.
%
% It will save the single session data and output it to save to a general
% file, outside the function, if requested.
%
% Plotting options can be added at the end of the function.
%
% USAGE:
%   [TFR, cfg] = trialparsed_MTspectrogram(FT_data, condition, param, opt);
%   [TFR, cfg] = trialparsed_MTspectrogram(FT_data, condition, param, opt, alignName);
%
% INPUTS:
%   FT_data   - FieldTrip trial-parsed data for ONE alignment.
%   condition - condition struct (from NGL02_postPhy).
%   param     - plotting/tuning params. .testname optional (see below).
%   opt       - resolved options; see the LFP-side entries in optSchema.
%   alignName - (optional) alignment tag ('itiOn' | 'stim2' | ...). When
%               supplied it is BAKED INTO the save filename and cache
%               scan, preventing the pre-26.06.2026 bug where NGL02_LFP
%               iterated opt.alignto and each iteration overwrote the
%               previous alignment's TFR file (silent data loss). When
%               omitted, behaves like the legacy call (single-align
%               callers only; discouraged).
%
% CACHE:
%   Existing cache read uses an EXACT filename match keyed by
%   (testname, session, method, alignName). Pre-26.06.2026 first-match
%   ls(...) bug is fixed here.
%
% Last modified 26.06.2026 (Jesus) - Pass 1: alignName arg + baked into
%                                     save filename + tightened cache
%                                     scan; internal artifact-rejection
%                                     re-run removed (NGL02_LFP does it
%                                     upstream on the whole FT_data).

%% Defaults
if nargin < 5 || isempty(alignName), alignName = ''; end
if ~isfield(param,'multi'),     param.multi     = false;        end
if ~isfield(param,'testname'),  param.testname  = 'trial_TFR_'; end
if ~isfield(opt,'blocks'),      opt.blocks      = 'all';        end
toload = 0;

% Alignment tag threaded into every emitted filename. Empty -> legacy
% behaviour (no _<align> suffix); NGL02_LFP always passes a non-empty
% tag, so the empty branch is for one-off standalone use.
if isempty(alignName)
    alignSlug = '';
else
    alignSlug = ['_' char(alignName)];
end

%% Check existence of multiple levels
if isstring(opt.blocks)
    opt.blocks = max(unique(condition.block));
    param.multi = true;
elseif ~isstring(opt.blocks)
    if opt.blocks > 1; param.multi = true; end
end

%% Check for processed data
% Cache scan is an EXACT match against the filename produced by this
% call, not the pre-26.06.2026 ls(<testname>*) first-match sweep that
% happily loaded another alignment's TFR into memory.
cacheStem = [param.testname, opt.SavFileName, alignSlug, '_', opt.TFRmethod];
cachePath = fullfile(opt.analysis, [cacheStem '.mat']);
if isfile(cachePath)
    toload = 1;
    fprintf('Loading existing TFR data (%s)...\n', cachePath);
    load(cachePath, 'TFR', 'cfg');
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
                seltrials = find(condition.correct == 1 & condition.block == b & condition.stimulus == (st-1) & ~isnan(FT_data.cfg.trl(:,3)));
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

                    % NOTE (Pass 1, 26.06.2026): the pre-Pass-1 code re-ran
                    % artifact_detRej_lfp here even though NGL02_LFP already
                    % handles that upstream on the whole FT_data before
                    % calling this function. That was double-filtering /
                    % double-zeroing. Removed. If you need per-block
                    % artifact treatment, do it in NGL02_LFP before the
                    % call, not inside the compute function.

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

    % Save TFR data. Filename baked with alignSlug so a subsequent
    % opt.alignto iteration cannot silently overwrite THIS alignment's
    % TFR (pre-26.06.2026 bug). cachePath was computed above so this
    % path matches the cache-read path exactly.
    save(cachePath, 'TFR', 'cfg', '-mat');
    fprintf('Saved TFR data -> %s\n', cachePath);
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