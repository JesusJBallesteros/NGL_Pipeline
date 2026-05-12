function out = binvec2dec(vec)
% binvec2dec  Convert LSB-first binary vector to decimal integer.
%
% PURPOSE:
%   Converts a binary row vector in LSB-first (least-significant-bit first)
%   order to its decimal equivalent. Used by INTAN_ExtractEvents to convert
%   INTAN digital-input pin states to event codes, following the NGL/Deuteron
%   convention where column 1 is pin 1 (LSB).
%
% USAGE:
%   dec = binvec2dec(vec)
%
% INPUT:
%   vec   - (1 × N) numeric vector; non-zero values are treated as 1.
%           Column 1 = LSB (pin 1), column N = MSB (pin N).
%           Maximum N = 52 (MATLAB bin2dec limit).
%
% OUTPUT:
%   out   - scalar decimal integer in range [0, 2^N - 1].
%
% EXAMPLES:
%   binvec2dec([1 0 0 0])  → 1   (itiOn from Deuteron pin scheme)
%   binvec2dec([1 1 0 0])  → 3   (bhv)
%   binvec2dec([1 1 1 1])  → 15  (end3)
%   binvec2dec([1 1 1 0 1]) → 23
%
% NOTES:
%   - Non-zero values map to 1: [1 2 3 0] is treated as [1 1 1 0].
%   - See also: dec2binvec, bin2dec.
%
%   Original: MP 11-11-98, Copyright 1998-2003 The MathWorks, Inc.
%   Modified for NGL event-coding use.

% Error if B is not defined.
% Non-zero values map to 1.
vec = vec~=0;

% Convert the binvec [0 0 1 1] to a binary string '1100';
h = deblank(num2str(fliplr(vec)'))';

% Convert the binary string to a decimal number.
out = bin2dec(h);