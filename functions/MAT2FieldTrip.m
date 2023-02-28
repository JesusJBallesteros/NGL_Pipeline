function MAT2FieldTrip(input, data, sessions, ss, opt, varargin)
% Wraps the process to transform a simple .mat file into one with
% appropiate format for further processing with FieldTrip toolbox.
% Options are, to create a 'continuous' FT file, (one, large trial) or to
% create a trial-parsed FT file, for which we need the eventcodes.

if nargin < 5,  stream = 'cont'; 
else,           stream = 'parsed';
                EventRecord = varargin{1};
end


if ~isfield(opt,'PathRaw'),           opt.PathRaw           = pwd;                                 end
if ~isfield(opt,'FolderProcDataMat'), opt.FolderProcDataMat = sessions.info{ss}.savefolder;        end
if ~isfield(opt,'SavFileName'),       opt.SavFileName       = sessions.list(ss).name;              end

% Navigate to the processed folder
saveFolder = opt.FolderProcDataMat;
cd(saveFolder);

switch stream
    case strcmp(stream,'cont')
        % Check that Fieldtrip likes what we have (it should).
        FT_data = ft_checkdata(data, 'feedback' ,'yes');
        clear data
        
        % Then give the FT_data a proper 'continous' state.
        cfg = [];
        cfg.continuous = 'yes';
        
        FT_data = ft_redefinetrial(cfg, FT_data);
        clear cfg
        
        % Save this session data. Generates a file with continous data for a
        %   SINGLE session only into the session folder.
        save(strcat(sessions.list(ss).name,'_continous_FT.mat'), ...
            'sessions', 'input', 'FT_data', '-v7.3')

    case strcmp(stream,'parsed')

        % TODO
end

end