function MAT2FieldTrip(data, opt, varargin)
% Wraps the process to transform a simple .mat file into one with
% appropiate format for further processing with FieldTrip toolbox.
% Options are, to create a 'continuous' FT file, (one, large trial) or to
% create a trial-parsed FT file, for which we need the eventcodes.
%
% Version 01.03.2023 Jesus

if nargin < 3,  stream = 1; 
else           
    trialdef = varargin{1};
    if isempty(trialdef), stream = 1;
    else, stream = 2;
    end
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
        disp('Data had no trial definition, so it was saved as continuous.');
        save(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_continous_FT.mat')), 'FT_data', '-v7.3')

    case 2
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
        disp('Data had a trial definition, so it was trial parsed and saved.');
        save(fullfile(opt.FolderProcDataMat, strcat(opt.SavFileName,'_tparsed_FT.mat')), 'FT_data', '-v7.3')

end

end