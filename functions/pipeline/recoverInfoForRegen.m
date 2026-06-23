function info = recoverInfoForRegen(input, opt)
% recoverInfoForRegen  Build a minimum `info` struct for re-running NGL01
%                      against a preprocessed-only data tree.
%
% PURPOSE:
%   When opt.regenFrom.preproc = true the raw folder is intentionally
%   empty (you moved the curated outputs to a different machine and
%   didn't bring the raw .dat / .rhd files). chckV can't probe the raw
%   folder; this helper produces a minimum-viable `info` from the
%   user's declared system + opt.numChannels + per-system defaults, so
%   that EventProcess, trialdefGen and conditions_script can still run.
%
% USAGE:
%   info = recoverInfoForRegen(input, opt)
%   Called from prepforsession when opt.regenFrom.preproc = true.
%
% INPUTS:
%   input         - struct from set_default. Required:
%                     .run                 (1x2) [subject_idx, session_idx]
%                     .subjects(x).name
%                     .sessions(x).list{y}
%                     .processed           preprocessing root
%                     .analysis            analysis root
%   opt           - resolved options struct. Required:
%                     .numChannels
%                     .regenFrom.system    'INTAN' | 'Deuteron'
%                     .regenFrom.fileformat  ('' -> derive from .system)
%                     .regenFrom.sample_rate ([] -> derive from .system)
%                     .regenFrom.nChannels   ([] -> use opt.numChannels)
%
% OUTPUT (struct, shape mirrors chckV's INTAN return):
%   .fileformat              char (e.g. 'fileperch' / 'DF1')
%   .amplifier_sample_rate   scalar, Hz
%   .nChannels               scalar
%   .HDF5chunkSize           300 * sample_rate
%   .numADCBits / .voltageRes  [] (unset; raw-only consumers don't run)
%   .files                   dir-struct of *_FTcont.mat (or empty);
%                            only LFP_Fieldtrip consumes this in regen
%   .fromPreproc             true (downstream defensive flag)
%   .regenSystem             echo of opt.regenFrom.system
%
% SYSTEM DEFAULTS:
%   'INTAN'    -> fileformat = 'fileperch', sample_rate = 30000
%   'Deuteron' -> fileformat = 'DF1',       sample_rate = 32000
%   These come from the lab's two recorder choices. Override the
%   defaults per-run via opt.regenFrom.fileformat / .sample_rate when
%   you have a rare DT2 / non-standard rate.
%
% ASSERT:
%   Errors with NGL:regen:noEventRecord if the session has no
%   EventRecord.mat to regenerate from.
%
% SEE ALSO:
%   chckV (the raw-folder counterpart this replaces in regen mode),
%   prepforsession, NGL01_Main, savePreprocInfo.
%
% Last modified 18.06.2026 (Jesus)

    subject = input.subjects(input.run(1)).name;
    session = input.sessions(input.run(1)).list{input.run(2)};
    sess    = fullfile(input.processed, subject, session);
    erFile  = fullfile(sess, 'EventRecord.mat');
    assert(isfile(erFile), 'NGL:regen:noEventRecord', ...
        ['opt.regenFrom.preproc=true but session %s has no EventRecord.mat ', ...
         '— nothing to regenerate from. Either disable regenFrom.preproc ', ...
         'or transfer EventRecord.mat into %s.'], ...
        fullfile(subject, session), sess);

    % Per-system defaults; overrides win when non-empty.
    sys = opt.regenFrom.system;
    switch upper(sys)
        case 'INTAN'
            fmt = 'fileperch'; fs = 30000;
        case 'DEUTERON'
            fmt = 'DF1';       fs = 32000;
        otherwise
            error('NGL:regen:badSystem', ...
                'opt.regenFrom.system must be ''INTAN'' or ''Deuteron'' (got ''%s'').', sys);
    end
    if ~isempty(opt.regenFrom.fileformat),  fmt = opt.regenFrom.fileformat; end
    if ~isempty(opt.regenFrom.sample_rate), fs  = opt.regenFrom.sample_rate; end

    if ~isempty(opt.regenFrom.nChannels)
        nCh = opt.regenFrom.nChannels;
    else
        nCh = opt.numChannels;
    end

    info                       = struct();
    info.fileformat            = fmt;
    info.amplifier_sample_rate = fs;
    info.nChannels             = nCh;
    info.HDF5chunkSize         = 300 * fs;
    info.numADCBits            = [];
    info.voltageRes            = [];
    info.fromPreproc           = true;
    info.regenSystem           = sys;

    % info.files is only consumed by LFP_Fieldtrip in regen mode (the
    % raw-side INTAN/Deuteron wrappers are skipped entirely). Auto-
    % populate from <analysis>/<subj>/<sess>/*_FTcont.mat so re-running
    % LFP work just works if the FT continuous file is on disk; empty
    % struct otherwise so isempty(info.files) reads cleanly.
    sessAnalysis = fullfile(input.analysis, subject, session);
    ftHits       = dir(fullfile(sessAnalysis, '*_FTcont.mat'));
    if ~isempty(ftHits)
        info.files = ftHits;
    else
        info.files = struct([]);
    end

    fprintf(['recoverInfoForRegen [%s/%s]: system=%s -> fileformat=%s, ', ...
             'sample_rate=%g Hz, nChannels=%d\n'], ...
            subject, session, sys, info.fileformat, info.amplifier_sample_rate, ...
            info.nChannels);
end