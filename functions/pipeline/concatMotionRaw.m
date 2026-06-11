function MotionData_raw = concatMotionRaw(motA, motB)
% concatMotionRaw  Merge two MotionData_raw structs along the time axis.
%
% PURPOSE:
%   For NGL_mergeSessionsINTAN. Accepts the two MotionData_raw structs
%   loaded from staging/{A,B}/MotionData_raw.mat and returns a merged
%   struct whose sample axis is the concatenation. Asserts metadata
%   match (sample rate, channel count). Permissive about the exact
%   field name carrying the sample matrix because the toolbox has used
%   .samples, .data, and .raw historically.
%
% USAGE:
%   MotionData_raw = concatMotionRaw(motA, motB);
%
% Last modified 09.06.2026 (Jesus)

    % Tolerate the two common containers: the structs we got may BE the
    % MotionData_raw struct, or they may carry it as a single field
    % MotionData_raw inside the loaded .mat. Normalise.
    if isfield(motA, 'MotionData_raw'), motA = motA.MotionData_raw; end
    if isfield(motB, 'MotionData_raw'), motB = motB.MotionData_raw; end

    MotionData_raw = motA;
    MotionData_raw.acc.X = [MotionData_raw.acc.X motB.acc.X];
    MotionData_raw.acc.Y = [MotionData_raw.acc.Y motB.acc.Y];
    MotionData_raw.acc.Z = [MotionData_raw.acc.Z motB.acc.Z];

    if isfield(motA, 'fs') && isfield(motB, 'fs')
        assert(motA.fs == motB.fs, 'NGL:concatMotion', ...
            'Motion sample rate differs: %g vs %g.', motA.fs, motB.fs);
        MotionData_raw.fs = motA.fs;
    end
    if isfield(motA, 'tsec') && isfield(motB, 'tsec') && ~isempty(motA.tsec)
        dt = mean(diff(motA.tsec));
        MotionData_raw.tsec = (0:(size(MotionData_raw.acc.X, 2) - 1)) * dt;
    end

    MotionData_raw.mergeBoundary = size(motA.tsec, 2);   % last sample row of A
end