function MAT2FieldTrip(data, opt, varargin)
% Wraps the process to transform a simple .mat file into one with
% appropiate format for further processing with FieldTrip toolbox.
% Options are, to create a 'continuous' FT file, (one, large trial) or to
% create a trial-parsed FT file, for which we need the eventcodes.
%
% INPUT: data, data struct as pseudo-FT format (from previous step)
%        opt, parameters and paths (needed only to save files)
%        trialdef, a (trials,3) array with start/end/t0 times. Optional
%        cont, logic, to force continuous data treatment if true. Optional.
%
% OUTPUT: 

% Jesus 12.06.2024

%% Check existence of FT files
if isempty(data)
    disp('Fieldtrip proper formatting was skipped too.')
    return
end

%% Check inputs
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

%% If data is loaded as a result of a failed run (to deprecate)
if isfield(data,"FT_data")
    data = data.FT_data;
end

%% Continuous treatment (trialdef missing or empty)
if cont
    % Check that Fieldtrip likes what we have (it should).
    FT_data = ft_checkdata(data);
    
    % Then give the FT_data a proper 'continous' state.
    cfg = [];
    cfg.continuous = 'yes';
    
    FT_data = ft_redefinetrial(cfg, FT_data);
    clear cfg
    
    % Save this session data.
    save(fullfile(opt.trialSorted, strcat(opt.SavFileName,'_FTcont.mat')), 'FT_data', '-v7.3');
    clear data
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
        % Then proceed to trial-parse the FT_data. Use 'ft_redefinetrial'
        cfg = [];
        cfg.trl = trialdef{2,i};

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