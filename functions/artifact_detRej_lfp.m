function [FT_data_NoArtif] = artifact_detRej_lfp(FT_data, condition, opt, param)
    cfg = [];
    cfg.trl                         = FT_data.trial(~isnan(FT_data.cfg.trl(:,3))');
    cfg.artfctdef.zvalue.channel    = 'all';
    cfg.artfctdef.zvalue.cutoff     = 15;
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;                
    
    % The optional configuration settings (see below) are:
      cfg.artfctdef.zvalue.artfctpeak       = 'no';
      cfg.artfctdef.zvalue.interactive      = 'no';
      cfg.artfctdef.zvalue.zscore           = 'yes';
    
    [~, artifact] = ft_artifact_zvalue(cfg, FT_data);

    % The following configuration options are supported
    cfg = [];
      cfg.artfctdef.reject          = 'partial';
      cfg.artfctdef.zvalue.artifact = artifact;
      cfg.artfctdef.value           = 'nan';

    [FT_data_NoArtif] = ft_rejectartifact(cfg, FT_data);
end