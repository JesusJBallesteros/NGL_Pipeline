function LFP_Fieldtrip(neurons, spike, trialdef, input, opt)
% LFP_Fieldtrip  DEPRECATED (moved to functions/_deprecated/ on 26.06.2026).
%
% WHY DEPRECATED:
%   * `if ~isfield('artifact_rejection', opt, ...)` (lines below) — arg
%     order swapped, always false; three misleading defaults.
%   * `if op.artifact_rejection` — typo (`op` vs `opt`), runtime NGL:undefined
%     if reached.
%   * Loads `input.sessions(...).info.files.name` from `opt.analysis` —
%     but FTcont lives in `opt.FolderProcDataMat`, not `opt.analysis`.
%   * References `opt.chgDtctPCue`, `opt.FLIP`, `opt.artifact_rejection`
%     — none are in the schema.
%   * Not called from anywhere in the pipeline. Stale scaffolding from
%     the single-project era.
%
% WHAT TO USE INSTEAD:
%   Per-session quick-look LFP  -> NGL02_LFP (schema-driven, area-aware).
%   vFLIP spectrolaminar mapping -> functions/video/vFLIP_NGL.m (planned
%       to be wired into a caller as part of the NGL07_LFPanalysis rollout).
%   Research-grade LFP analyses  -> NGL07_LFPanalysis (planned; will
%       cover trial-parsed TFR, oscillation detection, phase, spike-field,
%       LFP-behavior regression).
%
% This file is kept only so a lingering `LFP_Fieldtrip(...)` call in an
% old project script fails with a discoverable location rather than a
% missing-function error.
%
% Last modified 26.06.2026 (Jesus) - deprecation header added.
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