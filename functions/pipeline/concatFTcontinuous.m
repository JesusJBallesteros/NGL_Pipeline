function FT = concatFTcontinuous(ftA, ftB)
% concatFTcontinuous  Merge two FieldTrip "continuous" data structs.
%
% PURPOSE:
%   For NGL_mergeSessionsINTAN. Both inputs are single-trial continuous
%   FT data structures (as produced by ft_preprocessing on a raw INTAN
%   session). Returns a single struct whose trial axis holds [A B]
%   concatenated along time, with a continuous time vector and updated
%   sampleinfo / hdr.nSamples.
%
% ASSERTS:
%   .fsample, .label, channel count must match between A and B.
%
% USAGE:
%   FT = concatFTcontinuous(ftA, ftB);
%
% NOTES:
%   - Only the first cells of .trial / .time are used (FT continuous
%     convention is a single "trial" that spans the whole recording).
%   - .cfg is taken from A; B's preprocessing config is dropped to keep
%     the output struct lightweight. The merge step is recorded in
%     mergeMeta.mat alongside this file.
%
% Last modified 09.06.2026 (Jesus)

    assert(isstruct(ftA) && isstruct(ftB), 'NGL:concatFT', 'Inputs must be FT structs.');
    assert(isfield(ftA,'fsample') && isfield(ftB,'fsample'), ...
        'NGL:concatFT', 'Both inputs must have .fsample.');
    assert(isequal(ftA.label, ftB.label), ...
        'NGL:concatFT', 'Channel labels differ between A and B.');

    nA = size(ftA.trial{1}, 2);
    nB = size(ftB.trial{1}, 2);
    dt = 1 / ftA.fsample;

    FT          = ftA;
    FT.trial{1} = [ftA.trial{1}, ftB.trial{1}];
    FT.time{1}  = (0:(nA + nB - 1)) * dt;

    if isfield(FT,'sampleinfo'), FT.sampleinfo = [1, nA + nB]; end
    if isfield(FT,'hdr') && isstruct(FT.hdr)
        if isfield(FT.hdr,'nSamples'),  FT.hdr.nSamples  = nA + nB; end
        if isfield(FT.hdr,'nTrials'),   FT.hdr.nTrials   = 1;       end
    end

    % Annotate the merge boundary so downstream FT calls can detect it.
    FT.cfg.mergeBoundary = nA;   % last sample of A; B starts at nA+1.
end
