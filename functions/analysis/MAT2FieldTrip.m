function MAT2FieldTrip(data, opt, varargin)
% MAT2FieldTrip  Convert pseudo-FieldTrip struct to proper FT format and save.
%
% PURPOSE:
%   Takes the pseudo-FieldTrip struct produced by intan2MAT_wrapper and
%   calls ft_checkdata + ft_redefinetrial to produce:
%   (1) A continuous FT file (*_FTcont.mat) — always produced unless skipped
%   (2) Trial-parsed FT files (*_<event>.mat) — one per alignment event in trialdef
%   Both outputs are saved to opt.trialSorted. Skips if the continuous file
%   already exists (re-run safety).
%
% USAGE:
%   MAT2FieldTrip(data, opt)                    % continuous only
%   MAT2FieldTrip(data, opt, trialdef)          % continuous + trial-parsed
%   MAT2FieldTrip(data, opt, trialdef, cont)    % explicit control
%     cont = true  → force continuous output even if trialdef is supplied
%
% INPUTS:
%   data      - pseudo-FieldTrip struct from intan2MAT_wrapper (.label,
%               .trial, .time, .sampleinfo)
%   opt       - options struct; must contain:
%                 .trialSorted    output folder for .mat files
%                 .SavFileName    session name for output filenames
%   trialdef  - (optional) (2 × nAlignments) cell array from trialdefGen:
%                 {1,i} alignment event name (char)
%                 {2,i} [nTrials × 3] double [start end t0] in MILLISECONDS
%   cont      - (optional, logical) if true, force continuous processing
%               even when trialdef is provided (produces both outputs)
%
% OUTPUTS (saved to opt.trialSorted):
%   <SavFileName>_FTcont.mat          continuous FieldTrip data (FT_data)
%   <SavFileName>_<event>.mat         trial-parsed FieldTrip data per alignment
%
% NOTES:
%   - trialdef times are in ms; conversion to samples uses fs_lfp derived
%     from data.time{1} directly (not from opt.dwnsmplRate).
%   - cfg.trl offset is set so that FieldTrip's t=0 aligns to the t0 column
%     of trialdef: cfg.trl(:,3) = cfg.trl(:,1) - round(t0/1000 * fs_lfp).
%
% CALLS:
%   ft_checkdata, ft_redefinetrial (FieldTrip)
%
% Jesus 12.06.2024

%% Check input emptyness as a purpously skipped step
if isempty(data)
    disp('Fieldtrip proper formatting was skipped too.')
    return
end

%% If data is loaded, but is the result of a failed run 
if isfield(data,"FT_data")
    data = data.FT_data;
end

%% Check inputs, decide subprocessing
if nargin < 3
   cont = true; % do not trial parse
   disp('No trial definition was given. Data treated as continuous.');
   trialdef = [];
elseif nargin == 3
    trialdef = varargin{1};
    if isempty(trialdef)
        cont = true; % do not trial parse
        disp('Trial definition was empty. Data treated as continuous.');
    else
        cont = false; % trial parse
        disp('Trial definition found. Data will be trial-parsed.');        
    end
elseif nargin == 4
    trialdef = varargin{1};
    cont = varargin{2};
    disp('Trial definition found but continuous treatment forced. Data will be process both ways.');        
end

%% Continuous treatment (trialdef missing or empty)
if cont
    if isfile(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_FTcont.mat')))
        disp('A Fieldtrip-formatted file found in this directory, skipping.')
    else
        % Check that Fieldtrip likes what we have (it should).
        FT_data = ft_checkdata(data);
        
        % Then give the FT_data a proper 'continous' state.
        cfg = [];
        cfg.continuous = 'yes';
        
        FT_data = ft_redefinetrial(cfg, FT_data);
        clear cfg
        
        % Save this session data.
        save(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_FTcont.mat')), 'FT_data', '-v7.3');
        clear data
    end
end

%% Trial-parsed treatment 
if ~isempty(trialdef) 
    % Check that Fieldtrip likes what we have (if not forced before).
    if ~exist('FT_data',"var")
        FT_data_cont = ft_checkdata(data);
    else
        FT_data_cont = FT_data;
        FT_data = [];
    end

    % We may have more than one event to align things to.
    for i=1:size(trialdef,2)
        if isfile(fullfile(opt.trialSorted, strcat(opt.SavFileName, '_', trialdef{1,i} ,'.mat')))
            disp('A Fieldtrip trialparsed file found, skipping...')
        else
            % get sample rate
            fs_lfp = 1 / (data.time{1}(2) - data.time{1}(1));

            % Then proceed to trial-parse the FT_data. Use 'ft_redefinetrial'
            cfg = [];
            % cfg.trl = trialdef{2,i}; trial boundaries passed in ms instead of samples
            cfg.trl = round(trialdef{2,i} / 1000 * fs_lfp);
    
            % Re-set the offset of the trial definition for FT to get it.
            cfg.trl(:,3) = cfg.trl(:,1) - cfg.trl(:,3);
    
            % Redefine trials
            FT_data = ft_redefinetrial(cfg, FT_data_cont);
        
            % Update FT header info manually
            FT_data.hdr.nTrials = length(FT_data.trial);
        
            % Save this session data.
            save(fullfile(opt.trialSorted, strcat(opt.SavFileName, '_', trialdef{1,i} ,'.mat')), 'FT_data', '-v7.3')
        end
    end

end

end