function mat2FieldTrip(data, opt, varargin)
% Wraps the process to transform a simple .mat file into one with
% appropiate format for further processing with FieldTrip toolbox.
% Options are, to create a 'continuous' FT file, (one, large trial) or to
% create a trial-parsed FT file, for which we need the eventcodes.
%
% Version 01.03.2023 Jesus

if nargin < 3,  stream = 1; 
else,           stream = 2;
                EventRecord = varargin{1};
end

switch stream
    case 1
        % Check that Fieldtrip likes what we have (it should).
        FT_data = ft_checkdata(data);
        clear data
        
        % Then give the FT_data a proper 'continous' state.
        cfg = [];
        cfg.continuous = 'yes';
        
        FT_data = ft_redefinetrial(cfg, FT_data);
        clear cfg
        
        % Save this session data. Generates a file with continous data for a
        %   SINGLE session only into the session folder.
        save(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_continous_FT.mat')), 'FT_data', '-v7.3')

    case 2

        % TODO
end

end