function MAT2FieldTrip(data, opt, varargin)
% Wraps the process to transform a simple .mat file into one with
% appropiate format for further processing with FieldTrip toolbox.
% Options are, to create a 'continuous' FT file, (one, large trial) or to
% create a trial-parsed FT file, for which we need the eventcodes.

% Jesus 05.01.2024

if nargin < 3
   cont = true; % do not trial parse
   disp('No trial definition was given. Data treated as continuous.');
else           
    trialdef = varargin{1};
    if isempty(trialdef)
        cont = true; % do not trial parse
        disp('Trial definition was empty. Data treated as continuous.');
    else
        cont = false; % trial parse % TODO
        disp('Trial definition found. Data will be trial-parsed.');        
    end
end

if isfield(data,"FT_data")
    data = data.FT_data;
end

if cont 
    % Check that Fieldtrip likes what we have (it should).
    FT_data = ft_checkdata(data);
    clear data
    
    % Then give the FT_data a proper 'continous' state.
    cfg = [];
    cfg.continuous = 'yes';
    
    FT_data = ft_redefinetrial(cfg, FT_data);
    clear cfg
    
    % Save this session data. Generates a file with continous data for a
    % SINGLE session only into the session folder.
    save(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_continous_FT.mat')), 'FT_data', '-v7.3')

else
    % Check that Fieldtrip likes what we have (it should).
    FT_data = ft_checkdata(data);
    clear data

    % Let's put trialdef times into ms and round to the nearest integer
    trialdef = round(trialdef*1000);

    cfg = [];
    cfg.trl = trialdef;
    FT_data = ft_redefinetrial(cfg, FT_data);

    FT_data.hdr.nTrials = length(FT_data.trial);

    % Save this session data. Generates a FT file with trialparsed data.
    save(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_tparsed_FT.mat')), 'FT_data', '-v7.3')

end

end