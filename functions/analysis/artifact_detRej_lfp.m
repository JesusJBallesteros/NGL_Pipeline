function [FT_data_NoArtif] = artifact_detRej_lfp(FT_data, opt)
% artifact_detRej_lfp  z-value artifact detection and rejection on
%                      FieldTrip LFP data.
%
% PURPOSE:
%   Thin wrapper around FieldTrip's ft_artifact_zvalue + ft_rejectartifact
%   pair. Detects samples whose z-score exceeds opt.artZvalue and
%   replaces them with opt.rejValue (default 'zero' because TFR
%   downstream dislikes NaN). Used by NGL02_LFP when opt.artifdet=true.
%
% USAGE:
%   FT_data = artifact_detRej_lfp(FT_data, opt)
%
% INPUTS:
%   FT_data - FieldTrip data struct. Required: FT_data.cfg.trl.
%   opt     - resolved options struct. Used fields (defaulted by set_default):
%               .artZvalue  z-value cutoff for ft_artifact_zvalue.
%               .rejValue   replacement for rejected samples
%                           ('zero'|'nan'|numeric scalar).
%
% OUTPUT:
%   FT_data_NoArtif - FT_data with artifacts replaced.
%
% Last modified 23.07.2026 (Jesus)

%% Proceed
cfg = [];
    cfg.trl                         = FT_data.cfg.trl;
    cfg.artfctdef.zvalue.channel    = 'all';
    cfg.artfctdef.zvalue.cutoff     = opt.artZvalue;
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;    
    cfg.artfctdef.zvalue.continuous = 'no';
    
    % The optional configuration settings (see below) are:
      cfg.artfctdef.zvalue.artfctpeak       = 'yes';
      cfg.artfctdef.zvalue.interactive      = 'no';
      cfg.artfctdef.zvalue.zscore           = 'yes';
    
    [~, artifact] = ft_artifact_zvalue(cfg, FT_data);

    % The following configuration options are supported
    cfg = [];
      cfg.artfctdef.reject          = 'partial';
      % cfg.artfctdef.value           = 'nan';
      cfg.artfctdef.zvalue.artifact = artifact;

    [FT_data_NoArtif] = ft_rejectartifact(cfg, FT_data);
end