function [FT_data_NoArtif] = artifact_detRej_lfp(FT_data, opt)
% A basic wrapper for the FieldTrip functions 'ft_artifact_zvalue' and
% 'ft_rejectartifact' which typically go together. Outputs the final
% FT structure free of artifacts where those have been replaced with zeros
% (apparently TFR calculations don't like NaNs)
%
% Jesus 15.04.2025

%% Default
if ~isfield(opt,'artZvalue'),         opt.artZvalue             = 10;                    end
if ~isfield(opt,'rejValue'),          opt.rejValue              = 'zero';                end

%% Proceed
cfg = [];
    cfg.trl                         = FT_data.cfg.trl;
    cfg.artfctdef.zvalue.channel    = 'all';
    cfg.artfctdef.zvalue.cutoff     = opt.artZvalue;
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;                
    
    % The optional configuration settings (see below) are:
      cfg.artfctdef.zvalue.artfctpeak       = 'yes';
      cfg.artfctdef.zvalue.interactive      = 'no';
      cfg.artfctdef.zvalue.zscore           = 'yes';
    
    [~, artifact] = ft_artifact_zvalue(cfg, FT_data);

    % The following configuration options are supported
    cfg = [];
      cfg.artfctdef.reject          = opt.rejValue;
      cfg.artfctdef.zvalue.artifact = artifact;

    [FT_data_NoArtif] = ft_rejectartifact(cfg, FT_data);
end