function [TFR, cfg] = trialparsed_MTspectrogram(FT_data, condition, param, opt, alignName)
% trialparsed_MTspectrogram  Dispatch shim (project-specific vs generic).
%
% PURPOSE (as of LFP Pass 2, 26.06.2026):
%   This used to be the monolithic ASL-specific compute function. It is
%   now a THIN DISPATCHER that picks the right implementation based on
%   opt.proj_socialLearning:
%
%       opt.proj_socialLearning = true  -> computeTrialparsedTFR_ASL
%           (functions/analysis/projects/socialLearning/), which runs
%           the block-early/late x NS/FS x NS-FS-subtraction pipeline
%           that ONLY makes sense on the extintion task.
%
%       opt.proj_socialLearning = false -> computeTrialparsedTFR
%           (generic core), which returns per-band FT freq structs with
%           .powspctrm keeping per-trial detail; downstream analyses
%           (compareByBlock, planned NGL07_LFPanalysis regressions and
%           spike-field pipelines) work with that per-trial payload.
%
% USAGE:
%   [TFR, cfg] = trialparsed_MTspectrogram(FT_data, condition, param, opt);
%   [TFR, cfg] = trialparsed_MTspectrogram(FT_data, condition, param, opt, alignName);
%
% INPUTS / OUTPUTS:
%   See the docstrings of computeTrialparsedTFR / computeTrialparsedTFR_ASL
%   for the semantics of TFR / cfg on each branch.
%
% BACKWARD COMPATIBILITY:
%   NGL02_LFP calls this function directly; keeping the name preserves
%   that call site. Any project script that used to call
%   trialparsed_MTspectrogram will keep working - it just gets routed
%   to the correct implementation.
%
% SEE ALSO:
%   computeTrialparsedTFR (generic core),
%   computeTrialparsedTFR_ASL (project-specific ASL wrapper),
%   compareByBlock (generic block-contrast primitive that supersedes the
%       hardcoded early/late loops the ASL wrapper used to embed).
%
% Last modified 26.06.2026 (Jesus) 

    if nargin < 5, alignName = ''; end

    if isfield(opt,'proj_socialLearning') && opt.proj_socialLearning
        [TFR, cfg] = computeTrialparsedTFR_ASL(FT_data, condition, param, opt, alignName);
    else
        [TFR, cfg] = computeTrialparsedTFR(FT_data, condition, param, opt, alignName);
    end
end
