function [data,t,downsmpFactor] = downsampleVolt(data,FsInp,FsOut,dim,t,alignTime)
% downsampleVolt  Downsample voltage data by an integer factor.
%
% PURPOSE:
%   Reduces the sample rate of a multi-dimensional data array by retaining
%   every Nth sample along the specified dimension. No anti-aliasing filter
%   is applied here — the caller is responsible for low-pass filtering before
%   calling this function to prevent aliasing (intan2MAT_wrapper does this).
%
% USAGE:
%   [data, t, downsmpFactor] = downsampleVolt(data, FsInp, FsOut)
%   [data, t, downsmpFactor] = downsampleVolt(data, FsInp, FsOut, dim)
%   [data, t, downsmpFactor] = downsampleVolt(data, FsInp, FsOut, dim, t)
%   [data, t, downsmpFactor] = downsampleVolt(data, FsInp, FsOut, dim, t, alignTime)
%
% INPUTS:
%   data       - array of arbitrary size; data to downsample
%   FsInp      - (scalar, Hz) original sample rate
%   FsOut      - (scalar, Hz) target sample rate; must satisfy:
%                  FsInp/FsOut is an integer (asserted, no fractional resampling)
%   dim        - (optional, scalar) dimension along which to downsample
%                  default = first non-singleton dimension
%                  use dim=2 for [nChannels × nSamples] data
%   t          - (optional, [1 × nSamples]) original time vector to downsample
%   alignTime  - (optional, scalar) reference time point (e.g. 0) to preserve
%                  in the output; if [], starts at t(1)
%
% OUTPUTS:
%   data         - downsampled array (same ndims, 'dim' reduced)
%   t            - downsampled time vector (or [] if t was not supplied)
%   downsmpFactor - integer downsample factor = round(FsInp / FsOut)
%
% NOTES:
%   - Asserts FsInp > FsOut (upsampling not supported).
%   - Asserts FsInp/FsOut is an integer to floating-point tolerance (1e-9).
%   - For INTAN LFP: standard is FsInp=30000, FsOut=937.5 → factor=32.
%
% Originally: downsampleSignal (author unknown)

% todo  add filtering here?
 
  % If original and desired frequencies are equal (up to round-off error), we are done
  % (option may be used to match calls when downsampling is desired or not)
  if abs(FsInp - FsOut) < 1E-9,  return;  end
  
  if (nargin < 6),                  alignTime = []; end
  if (nargin < 5),                  t       = []; end
  if (nargin < 4) || isempty(dim),  dim     = firstNonSingletonDim(data); end

  nDims         = ndims(data);
  inputSize     = size(data);
  nT            = inputSize(dim);

  % Factor to downsample original data by
  downsmpFactor = FsInp / FsOut;
  
  % Handle error conditions in requested output sampling frequency
  assert(downsmpFactor > 1, ...
          ['downsampleSignal: Sampling frequency for final analysis (%d) < original data (%d Hz)\n', ... 
           'Change final sampling freq (or code up upsampling here)'], FsOut, FsInp);
  assert(rem(downsmpFactor,1) <= 1E-9, ...
          ['downsampleSignal: Sampling frequency for final analysis (%d) is not integer multiple of original data (%d Hz)\n' ... 
           'Change final sampling freq (or code up more sophisticated downsampling algorithm)'], FsOut, FsInp);

  downsmpFactor = round(downsmpFactor);             % Note: account for small round-off errors
  
  % Downsample every nth sample from 1:N
  if isempty(alignTime)
    tSmpIdxs  = [1 : downsmpFactor : nT];
    
  % Align downsampling to reference timepoint (eg t = 0)  
  else
    t0Idx     = find(abs(t - alignTime) < 1E-9);      % Note: acct for small round-off errors in t    
    tSmpIdxs  = [fliplr(t0Idx : -downsmpFactor : 1),  (t0Idx+downsmpFactor) : downsmpFactor : nT];              
  end    

  % Rearrange signals so 1st dimension is sampling (eg, time) dimension
  if dim ~= 1
    dimPerm   = [dim setxor([1:nDims],dim)];
    data      = permute(data, dimPerm);
  end
  
  % Downsample data and time-sampling vector
  data        = data(tSmpIdxs,:);
  if ~isempty(t)
    t         = t(tSmpIdxs);   
  end
  
  % Rearrange signals to original size and dimensional order
  if nDims > 2
    nT          = length(tSmpIdxs);
    reshapeSize = [nT inputSize(setxor(1:nDims,dim))];
    data        = reshape(data, reshapeSize);    
  end
  if dim ~= 1
    data        = ipermute(data, dimPerm);
  end
  
end

