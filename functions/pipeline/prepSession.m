function [input, opt] = prepSession(input, opt)
% prepSession  Per-session scaffolding shared by NGL02_postPhy and
%              NGL02_LFP. Wraps the three steps that always run at the
%              top of a session iteration: prepforsession, the NGL01
%              output pre-flight, and the preprocInfo opt overlay.
%
% PURPOSE:
%   Both NGL02_postPhy (spike) and NGL02_LFP (LFP) need the same per-
%   session preamble before doing their work: set the session paths on
%   opt, confirm NGL01 outputs exist, and overlay the preprocessing-time
%   opt snapshot. Centralising it here keeps the two scripts thin and
%   guarantees they handle these steps identically.
%
% USAGE:
%   [input, opt] = prepSession(input, opt)
%   Called once per iteration of the session loop, after input.run has
%   been set to [subject_idx session_idx].
%
% INPUTS:
%   input  - struct from set_default, with input.run = [x y] already set.
%   opt    - resolved options struct (post-set_default).
%
% OUTPUTS:
%   input  - updated by prepforsession (info field, etc.).
%   opt    - updated by prepforsession (session paths) then overlaid with
%            the saved-at-NGL01-time opt where available.
%
% STEPS:
%   1. prepforsession(input, opt)
%      Sets opt.FolderProcDataMat, opt.spikeSorted, opt.trialSorted,
%      opt.analysis, opt.KSfolder / opt.KSfolders, etc.
%   2. checkNGL01Outputs(input, opt)
%      Errors up front if NGL01 outputs are missing for THIS session.
%      Gated by opt.doSpikething / opt.doLFPthing so each calling script
%      only requires what it actually needs.
%   3. applyPreprocInfo overlay
%      Tries the per-session preprocInfo.mat first; falls back to the
%      master snapshot in analysisCode\; falls back to a warning if
%      neither exists. NGL01-owned fields (numChannels, alignto, filters,
%      etc.) are overwritten; NGL02-owned flags are preserved.
%
% SEE ALSO:
%   prepforsession, checkNGL01Outputs, loadPreprocInfo, applyPreprocInfo.
%
% Last modified 29.05.2026 (Jesus) - extracted from NGL02_postPhy during
%                                    the NGL02 spike/LFP split (#21)

%% 1. Per-session path resolution.
[input, opt] = prepforsession(input, opt);

%% 2. Pre-flight: confirm NGL01 outputs exist for this session.
%   Errors up front rather than letting failure surface deep inside
%   loadSpikes / sort2trials / FT_data loads.
checkNGL01Outputs(input, opt);

%% 3. Overlay preprocessing-time opt for this session.
%   Pull the per-session preprocInfo.mat written by NGL01 and overlay
%   NGL01-owned fields (numChannels, alignto, kilosort settings, filters,
%   etc.) onto the current opt. NGL02-owned flags (doSpikething,
%   doLFPthing, popDyn, ...) are preserved. Falls back to the master file
%   in analysisCode\, then to a warning if neither exists.
x = input.run(1); y = input.run(2);
subjStr = input.subjects(x).name;
sessStr = input.sessions(x).list{y};

try
    sessionInfo = loadPreprocInfo(opt.FolderProcDataMat, 'session');
    fprintf('prepSession: applying per-session preprocInfo (NGL01 run %s, branch %s)\n', ...
            sessionInfo.toolboxVersion.hash, sessionInfo.toolboxVersion.branch);
    opt = applyPreprocInfo(opt, sessionInfo);
catch ME
    if strcmp(ME.identifier, 'NGL:loadPreprocInfo:notFound')
        try
            masterInfo = loadPreprocInfo(input.analysisCode, 'master');
            warning('NGL02:fallbackMaster', ...
                ['Per-session preprocInfo not found for %s/%s. Falling back ', ...
                 'to analysisCode\\preprocInfo_lastRun.mat (NGL01 run %s, branch %s).'], ...
                subjStr, sessStr, ...
                masterInfo.toolboxVersion.hash, masterInfo.toolboxVersion.branch);
            opt = applyPreprocInfo(opt, masterInfo);
        catch
            warning('NGL02:noPreprocInfo', ...
                ['No preprocInfo snapshot found for %s/%s. Proceeding with opt ', ...
                 'as set in NGL_SetAndRunMe; if NGL01 was run with different ', ...
                 'acquisition/sorting options, results may be inconsistent.'], ...
                subjStr, sessStr);
        end
    else
        rethrow(ME);
    end
end
end
