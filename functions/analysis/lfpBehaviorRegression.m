function reg = lfpBehaviorRegression(TFR, behavior, opt) %#ok<INUSD>
% lfpBehaviorRegression  Cluster-permutation regression of TFR power on behavior.
%
% STATUS: STUB (LFP Pass 3, 26.06.2026).
%
%   The wiring for this analysis is planned but not usable end-to-end
%   yet because the UPSTREAM CONTRACT is not fulfilled: NGL06_videoAnalysis
%   currently emits only the rendered mp4/png. It does NOT save the
%   cleaned per-frame position + head-direction time series as a MATLAB
%   sidecar (.mat). Without that sidecar, this function has nothing to
%   regress LFP power against on a per-trial basis.
%
% WHEN THIS FUNCTION WILL BECOME USABLE:
%   1. NGL06 gains an emit-sidecar option (or the Python master saves
%      one by default). Sidecar path convention:
%           <csvDir>/<csvName>_gaze.mat
%      containing at minimum:
%           t         [nFrames x 1]  seconds from session start
%           x, y      [nFrames x 1]  arena coordinates (px or cm)
%           head_dir  [nFrames x 1]  radians (unwrapped or [-pi pi])
%           gaze_az   [nFrames x 1]  radians, per-eye if per-eye render
%      One update on the Python side + a MATLAB parser wrapper unlocks
%      this.
%   2. This function then:
%        a. Loads the sidecar for the session.
%        b. Resamples each covariate onto the TFR trial-time grid.
%        c. Builds a design matrix per (freq, time) bin.
%        d. Calls ft_freqstatistics with
%             cfg.statistic       = 'ft_statfun_depsamplesregrT'
%             cfg.correctm        = 'cluster'
%             cfg.numrandomization = 1000
%           per covariate.
%        e. Returns the cluster-thresholded stat map, cluster p-values,
%           and the covariate values themselves for reference.
%
% USAGE (post-completion):
%   reg = lfpBehaviorRegression(TFR, behavior, opt)
%
% INPUTS (contract; NOT yet enforced):
%   TFR      - ft_freqanalysis output with keeptrials='yes'.
%   behavior - struct from the planned NGL06 sidecar loader with
%                .t, .x, .y, .head_dir, .gaze_az (fields as above).
%   opt      - resolved options:
%                .lfp.behReg.covariates  cellstr of covariate names to
%                                         regress against (default:
%                                         {'x','y','head_dir'}).
%                .lfp.behReg.numrand     default 1000.
%                .lfp.behReg.alpha       default 0.05.
%
% CURRENT BEHAVIOUR:
%   Emits a single NGL:lfpBehaviorRegression:notImplemented warning and
%   returns an empty struct with a `.status = 'notImplemented'` field
%   and a `.blocked_on` note. NGL07_LFPanalysis silently no-ops this
%   analysis until the sidecar contract is fulfilled.
%
% SEE ALSO:
%   NGL06_videoAnalysis (upstream), functions/video/process_gaze.m
%   (would need the sidecar-save hook), NGL07_LFPanalysis (caller).
%
% Last modified 26.06.2026 (Jesus) - Pass 3 stub with upstream contract
%                                     documented; body deferred until
%                                     NGL06 emits the position/HD sidecar.

    warning('NGL:lfpBehaviorRegression:notImplemented', ...
        ['lfpBehaviorRegression is a stub. It becomes usable when ', ...
         'NGL06_videoAnalysis emits a per-frame <csvName>_gaze.mat ', ...
         'sidecar with (t, x, y, head_dir, gaze_az). Add the sidecar ', ...
         'save on the Python side (configfiles/master_gaze.py) then ', ...
         'fill in the body of this function.']);

    reg = struct( ...
        'status',     'notImplemented', ...
        'blocked_on', 'NGL06 gaze.mat sidecar with per-frame position + head_direction', ...
        'planned_output_fields', {{'stat','clusters','pvals','covariates'}} );
end
