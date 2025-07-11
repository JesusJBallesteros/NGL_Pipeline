function [spec,b,a] = bandFilter(data, dim, freqBands, smpRate, varargin)
%BANDFILTER 1D frequency band filtering on data array 
%
% [spec,b,a] = bandFilter(data,dim,freqBands,smpRate,varargin)
% 
% Performs phase-preserving 1D frequency-band filtering (using filtfilt()) and 
% (optional) Hilbert transform on time series (or other) data. Though the analysis 
% itself is 1D, it may be performed along a given dimension of multi-dim array,
% with the analysis repeated on every other dimension (ie, treating data array
% as set of multiple time series, eg trials).
% 
% Filter parameters may be set using either predefined set(s) of filter coefficients 
% (eg, via setFilterParams()) or given frequency band(s).
% 
% INPUT
% data      Array of arbitrary size. Set of time series (or other) data to band filter.
%           Can be any arbitrary size -- dimension 'dim' (length nSamples) is treated as 
%           the time samples for each time series to analyze, and all other dimensions 
%           are treated as independent time series (eg, distinct trials, channels, etc.) 
%           to be analyzed separately
% 
% dim       Scalar. Dimension of data array corresponding to time samples of each series.
%           Default: 1st non-singleton dimension
%
% freqBands [nFreqBands x 2]. [low,high]-cut frequencies for each frequency band to use.
%           Note: must input either (freqBands & smpRate) OR (a & b) to give parameters for filtering
%
% smpRate   Scalar. Signal sampling rate (Hz).
%
% OPTIONAL INPUT (given as name/value pairs in varargin)     
% nPad      Scalar. Number of time points to zero-pad time series with before filtering.
%           This helps prevent edge effects at the end of the time series. Default: 0
%
% removeDC  Logical. Remove mean across timepts in each time series before filtering.
%           Default: true
%
% direction String. Direction to run filter on data. Default: 'both'
%           'both'    : Filters in both directions to perform 0-phase filtering (using filtfilt())
%           'forward' : Filters data only in forward direction (using filter())
%           'reverse' : Filters data only in reverse direction (using filter())
%
% signalType String. Type of signal to return. Default: 'real'. Options:
%           'real'    : Band-filtered (real-valued) signal
%           'complex' : Hilbert transform filtered data to get complex analytic signal -- hilbert(dataFilt)
%           'imag'    : Imaginary projection of data (via Hilbert) -- imag(hilbert(dataFilt))
%           'phase'   : Phase of band-filtered data -- angle(hilbert(dataFilt))
%           'power'   : Power of band-filtered data -- hilbert(dataFilt).*conj(hilbert(dataFilt))
%           'magnitude' : Magnitude (envelope) of band-filtered data -- sqrt(hilbert(dataFilt).*conj(hilbert(dataFilt)))
%           'fullwave': Full-wave rectified band-filtered data -- abs(dataFilt)
% 
% b,a       {1 x nFreqBands cell}. Filter coefficients for one or more frequency bands.
%           Each cell element holds filter numerator (b), or denominator (a) of 
%           filter coefficients for each band:
%             [1 x 2*order+1] vector for band-pass filtered bands | 
%             [1 x order+1] vector for hi/lo-pass filtered bands
%           NOTE: If a and b are given, any input values for 'name', 'order', and 'type'
%                 are ignored (these are already implicit in a,b).
%
% (Pass-thru parameters passed from here to setFilterParams())
% name      String. Name of filter to implement. Currently supported options: 
%           'butter'(butterworth) [default],'ellip'(elliptical), 'fir' (FIR filter
%
% order     Scalar. Order of filter to implement. 
%           Default: 25 if name=='fir'; 3 otherwise
%
% type      String | {1 x nFreqBands} cell of strings.
%           Class of filter to implement for each band. Option:
%           'low'   : Low-pass filtering. Default for freqBands(f,1) == 0 or -Inf
%           'high'  : High-pass filtering. Default for freqBands(f,2) == smpRate/2 or +Inf
%           'band'/'bandpass' : Band-pass filtering. Default if above conditions not met
%           'stop'/'bandstop' : Band-stop filtering (eg for notch filtering).
% 
% OUTPUT
% spec      Array of arbitrary size. Band-filtered version of data, of given signalType
%           If nFreqBands == 1, size == same as input;
%             otherwise, frequency dimension (length==nFreqBands) is inserted before time/sampling dim
%             (eg, data [nSamples x nTrials], returns spec [nFreqBands x nSamples x nTrials]).
% 
% b,a       {1 x nFreqBands cell}. Filter coefficients for one or more frequency bands.
%
% From ETALO toolbox. Modified by Jesus.

  %% Process arguments
  if (nargin < 2) || isempty(dim),        dim       = firstNonSingletonDim(data); end  % Default is 1st non-singleton dimension 
  if (nargin < 3) || isempty(freqBands),  freqBands = []; end
  if (nargin < 4) || isempty(smpRate),    smpRate   = []; end  
  
  defaults  = { 'a',          [], ...
                'b',          [], ...
                'nPad',       0, ...                
                'removeDC',   true, ...
                'direction',  'both', ...
                'signalType', 'real' };
  
  [a,b,nPad,removeDC,direction,signalType,passThruArgs] = ...
    processArgs(defaults(1:2:end), defaults(2:2:end), varargin, false);
  direction = lower(direction);
  signalType= lower(signalType);
  
  directionSupported = {'both','forward','reverse'};
  signalTypeSupported= {'real','complex','imag','phase','power','magnitude','fullwave'};
  
  assert(any(strcmp(direction,directionSupported)), ...
          [sprintf('bandFilter: Unknown filtering direction ''%s''. Supported directions:\n ',direction), ...
           sprintf(' ''%s'' ', directionSupported{:})]);
  assert(any(strcmp(signalType,signalTypeSupported)), ...
          [sprintf('bandFilter: Unknown signal type ''%s''. Supported types:\n ',signalType), ...
           sprintf(' ''%s'' ', signalTypeSupported{:})]);
         
  %% Set filter coefficients if not passed in as arg's
  if isempty(a) || isempty(b)
    assert(~isempty(freqBands), 'bandFilter: must input either filter coeff''s a,b or freqBands');
    assert(~isempty(smpRate), 'bandFilter: must set sampling rate to compute coefficients');

    [b,a]   = setFilterParams(freqBands,smpRate,passThruArgs{:});

  end
  
  if ~iscell(b), b = {b}; end
  if ~iscell(a), a = {a}; end  

  nSmp      = size(data,dim);               % Length (#samples) of unpadded time series
  n         = nSmp + nPad;                  % Total length of time series, including padding  
  nF        = length(b);                    % Number of frequency bands to use

  %% Data preprocessing
  % Rearrange data array data -> [nSamples x nDataSeries] matrix, s.t. each time series 
  %  (trial,channel,etc.) to transform is arranged as a column, and each column is a separate time series
  nDims     = ndims(data);
  inputSize	= size(data);                      % Original size of data array
    
  if dim ~= 1
    dimPerm = [dim setxor([1:nDims],dim)];
    data    = permute(data, dimPerm);
  end
  if nDims > 2
    data    = reshape(data, inputSize(dim), []);
  end
  nSeries   = size(data,2);                    % Number of independent time series to transform

  % Remove DC (mean across timepts for each data series, eg trial)
  if removeDC
    data    = bsxfun(@minus, double(data), mean(data,1)); % Note: equiv. to data = data - mean(data)
  end
  % Zero-pad each time series w/ nPad values
  if nPad ~= 0
    data    = [data; zeros(nPad,nSeries)];     
  end
   
  %% Do band filtering & optional signal tranformation (Hilbert transform, etc.) for each frequency band  
  if strcmpi(direction,'both'), filtFunc = @filtfilt;
  else,                         filtFunc = @filter;
  end
    
  % Initialize array for output of filtering (Note: array is permuted below -> [nFreqBands,n,nSeries])
  if strcmp(signalType,'complex'),  spec = complex( zeros(n,nSeries,nF) );
  else,                             spec = zeros(n,nSeries,nF);
  end
    
  for iF = 1:nF
    % Filter to get (real-valued) filtered signal for given frequency band
    if strcmp(direction,'both')
      spec(:,:,iF) = filtfilt(b{iF}, a{iF}, data);
      
    elseif strcmp(direction,'forward')
      spec(:,:,iF) = filter(b{iF}, a{iF}, data);
      
    else % implicit: if strcmpi(direction,'reverse')
      spec(:,:,iF) = flip(filtFunc(b{iF}, a{iF}, flip(data)));
    end
    
    %% Compute output signal type from band-filtered data
    switch signalType
      case 'real';      % do nothing; already have real filtered signal above
      case 'imag';      spec(:,:,iF) = imag(hilbert(spec(:,:,iF)));
      case 'complex';   spec(:,:,iF) = hilbert(spec(:,:,iF));      
      case 'phase';     spec(:,:,iF) = angle(hilbert(spec(:,:,iF)));      
      case 'magnitude'; spec(:,:,iF) = abs(hilbert(spec(:,:,iF)));
      case 'power';     spec(:,:,iF) = abs(hilbert(spec(:,:,iF))).^2;
      case 'fullwave';  spec(:,:,iF) = abs(spec(:,:,iF));                      
    end    
  end
 
  %% Data postprocessing  
  if nPad ~= 0 
    spec  = spec(1:nSmp,:,:);             % Remove padding
  end
  
  % If using > 1 freq band, insert frequency dimension just before time/sampling dim in output data array
  % Otherwise, return data array with same dimensionality/size as input
  if nF > 1
    % Permute dimensions -> [nFreqBands,nSamples,nDataSeries]
    spec = permute(spec, [3 1 2]);
  end
  % Rearrange output variable arrays so in format/size expected based on original input arg's
  if nDims > 2
    if nF > 1,  reshapeSize = [nF nSmp inputSize(setxor(1:nDims,dim))];
    else,       reshapeSize = [nSmp inputSize(setxor(1:nDims,dim))];
    end
    spec = reshape(spec, reshapeSize);
  end
  if dim ~= 1
    if nF > 1,  dimPerm = [2+(1:(dim-1)) 1 2 (dim+2):(nDims+1)]; 
    else,       dimPerm = [1+(1:(dim-1)) 1 (dim+1):nDims];
    end
    spec = permute(spec, dimPerm);
  end  
  
end