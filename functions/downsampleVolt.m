function [data,t,downsmpFactor] = downsampleVolt(data,FsInp,FsOut,dim,t,alignTime)
%DOWNSAMPLESIGNAL Downsamples signal to given sampling frequency
%
% [data,t,downsmpFactor] = downsampleSignal(data,FsInp,FsOut,dim,t,alignTime)
% 
% NOTE: Currently no low-pass filtering is performed here, 
% it's up to user to make sure no aliasing might occur (no signal > Nyquist limit)
% 
% INPUT
% data    Array (any arbitrary size). Data to downsample.
%
% FsInp   Scalar. Original data sampling frequency (Hz)
%
% FsOut   Scalar. Desired final sampling frequency after downsampling (Hz)
%
% dim     Scalar (optional). Sampling (eg, time) dimension of data. Downsampling will
%         be performed along this dim, all other treated as independent time series
%
% t       [1 x nTIn] (optional). Original sampled time points
%
% alignTime Scalar (optional). Can set reference timepoint w/in time vector (eg, t=0) 
%         to align sampling to. This sample will be always be preserved in output, 
%         and samples will be taken relative to it. Set = [] to just start at t(1).
%         Default = [] (no specific reference timepoint; start at t(1))
% 
% OUTPUT
% data    Array (same size as input, but with 'dim' downsampled). Downsampled data.
%
% t       [1 x nTIn]. Downsampled time sampling vector.
%
% downsmpFactor Scalar. Factor to downsample original data by (= FsInp / FsOut)

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

