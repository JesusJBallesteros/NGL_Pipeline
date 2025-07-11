function LFP_Fieldtrip(neurons, spike, trialdef, input, opt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if ~isfield('artifact_rejection',opt),  opt.artifact_rejection  = false; end
if ~isfield('FLIP',opt),                opt.FLIP                = false; end
if ~isfield('chgDtctPCue',opt),         opt.chgDtctPCue         = false; end

FT_data = [];

disp('Loading FT continuous file...')
load(fullfile(opt.analysis, input.sessions(input.run(1)).info.files.name));
if isfield(FT_data,"FT_data")
    FT_data = FT_data.FT_data;
end   % Simplify loaded structure if needed
    
% Obtain or create trial definition to pass to FT
if isfile('trialdef.mat')
    load(fullfile(opt.trialSorted, "trialdef.mat"));
else 
    if isfile(fullfile(opt.analysis, "events.mat"))
       load(fullfile(opt.analysis, "events.mat"));
       [~, trialdef, ~] = trialdefGen(events, opt, 1);
       save(fullfile(opt.trialSorted, "trialdef.mat"), 'trialdef');
    else
        warning('Neither trial definitions or events found for this session.')
    end
end

if opt.chgDtctPCue % Specific Step
  trialdef{2,1} = ceil(trialdef{2,1}/32);                
  MAT2FieldTrip(FT_data, opt, trialdef); 
  clear FT_data events trialdef
end

% If artifact retection
if op.artifact_rejection
    cfg = [];
     cfg.trl         = FT_data.cfg.trl;
     cfg.continuous  = 'no';
     cfg.artfctdef.zvalue.channel    = 'all';
     cfg.artfctdef.zvalue.cutoff     = 20;
     cfg.artfctdef.zvalue.trlpadding = 0;
     cfg.artfctdef.zvalue.fltpadding = 0;
     cfg.artfctdef.zvalue.artpadding = 0;                
     cfg.artfctdef.zvalue.artfctpeak       = 'no';
     cfg.artfctdef.zvalue.interactive      = 'no';
     cfg.artfctdef.zvalue.zscore           = 'yes';
        [~, artifact] = ft_artifact_zvalue(cfg, FT_data);

    % The following configuration options are supported
    cfg = [];
     cfg.artfctdef.reject          = 'partial';
     cfg.artfctdef.zvalue.artifact = artifact;
        [FT_data_art] = ft_rejectartifact(cfg, FT_data);
end

% vFLIP Analysis. Developed for chgDtctPCue
  if opt.FLIP
     laminaraxis = 0:0.05:1.55;
     freqaxis = 1:150;
     [FLIP, relpow, ~] = vFLIP_NGL(FT_data, laminaraxis, freqaxis, 0);
  end

end