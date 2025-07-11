function varargout = setFilterParams(freqBands,smpRate,varargin)
%SETFILTERPARAMS  Sets parameters for high/low/band-pass filtering
%
% varargout = setFilterParams(freqBands,smpRate,varargin)
% 
% Sets filter parameter for band-pass/low-pass/high-pass filtering
% in one or more given frequency band(s), with given filter type and order, 
% with given filter form/represenation, as might be passed into a 
% filter() or filtfilt() call.
%
% INPUT
% freqBands   [nFreqBands x 2]. [Low,High]-cut frequencies of each filter band.
%             For low-pass filtering, set freqBands(f,2) = Inf (or = smpRate/2)
%             For high-pass filtering, set freqBands(f,1) = 0
%             For band-pass filtering, set 0 < freqBands(f,:) < smpRate/2
%
% smpRate     Scalar. Sampling rate (Hz) of signal to filter.
%
% OPTIONAL INPUT (given as name/value pairs in varargin)
% name        String. Name of filter to implement. Default: 'butter'
%             Currently supported options:
%             'butter'(butterworth), 'ellip'(elliptical), 'fir' (FIR filter)
% 
% order       Scalar. Order of filter to implement. 
%             Default: 25 if name=='fir'; 3 otherwise
%
% type        String | {1 x nFreqBands} cell of strings.
%             Class of filter to implement for each band. Option:
%             'low'   : Low-pass filtering. Default for freqBands(f,1) == 0 or -Inf
%             'high'  : High-pass filtering. Default for freqBands(f,2) == smpRate/2 or +Inf
%             'band'/'bandpass' : Band-pass filtering. Default if above conditions not met
%             'stop'/'bandstop' : Band-stop filtering (eg for notch filtering).
%
% form        String. Which form of filter descriptor should we return. Options:
%             'transferFunction'/'b,a'  : Outputs transfer function coefficients [b,a]
%             'zero-pole-gain'/'z,p,k'  : Outputs zero-pole-gain representation [z,p,k]
%             '2ndOrderSections'/'sos,g'  : Outputs 2nd-order sections form [sos,g]
%             Default: 'zero-pole-gain' iff nargout==3, else 'transferFunction'
% 
% OUTPUT      Outputs depend on value of <form>
% For form = 'transferFunction':
% b,a         [1 x (2*)order+1] | {nFreqBands x 1} cell of [1 x (2*)order+1] vectors. 
%             Filter numerator (b) and denominator (a) coefficients based on given parameters.
%             If single freqBand input, a,b returned as pair of vector
%             If multiple freqBands input, a,b returned as cell vectors, each cell holding
%               coefficients (a or b) for the corresponding freqBand
%             For band-pass filtering, each coeff vector is size[1,,2*order+1]
%             For hi/lo-pass filtering, each coeff vector is size[1,,order+1]
%
% For form = 'zero-pole-gain':
% z,p,k       {nFreqBands x 1} cell. Zeros, poles, and gains of filters
%
% For form = '2ndOrderSections'
% sos,g       {nFreqBands x 1} cell. Second-order section representation and system gains of filters

  % Standardizes different "type" specifiers to standardized values (bottom row)
  typeNorm  = containers.Map({'low', 'lowpass', 'high', 'highpass', 'band',     'bandpass', 'stop', 'bandstop'}, ...
                             {'low', 'low',     'high', 'high',     'bandpass', 'bandpass', 'stop', 'stop'});
  typeSupported = unique(typeNorm.values);

  % Standardizes different "form" specifiers to standardized values (bottom row)
  formNorm  = containers.Map({'b,a', 'transferfunction',  'z,p,k', 'zero-pole-gain',  'sos,g',  '2ndordersections'}, ...
                             {'b,a', 'b,a',               'z,p,k', 'z,p,k',           'sos,g', 'sos,g',});
  formSupported = unique(formNorm.values);
  
  %% Process function arguments
  nFreqBands= size(freqBands,1);    
  
  defaults  = { 'name',       'butter', ... % Default: Butterworth filter
                'order',      [], ...   % Default set below based on filter name
                'type',       '', ...   % Defaults set below based on freqBands
                'form',       '', ...   % Default set below based on nargout
              };
  [name,order,type,form] = processArgs(defaults(1:2:end), defaults(2:2:end), varargin, false);    
  
  isTypeGiven = ~isempty(type);
  if ischar(type), type = {type}; end
  if isTypeGiven
    for iF = 1:nFreqBands
      assert(any(strcmpi(type{iF},typeNorm.keys)), ...
              [sprintf('setFilterParams: Unsupported filter type ''%s''. Supported values:\n  ',type{iF}), ...
               sprintf(' ''%s'' ',typeSupported{:})]);              
      type{iF} = typeNorm(lower(type{iF}));
    end
  end
  
  name  = lower(name);
  
  if isempty(order)
    if strcmp(name,'fir'),  order = 25;
    else,                   order = 3;
    end
  end
  
  if isempty(form)
    if nargout == 3,  form = 'z,p,k';
    else,             form = 'b,a';     
    end
  else
    assert(any(strcmpi(form,formNorm.keys)), ...
            [sprintf('setFilterParams: Unknown filter form ''%s''. Supported values:\n  ',form), ...
             sprintf(' ''%s'' ',formSupported{:})]);              
    form = formNorm(lower(form));    
  end
  
  % Parameter error checking
  assert((nargin >= 2) && ~isempty(smpRate) && ~isempty(freqBands), ...
    'setFilterParams: arguments must be given for both ''freqBands'' and ''smpRate''');  
  assert(all(freqBands(:,2) > freqBands(:,1)), 'setFilterParams: lo-cut freq must be < hi-cut freq');
  assert(all((freqBands(:,2) <= smpRate/2) | (freqBands(:,2) == Inf)), ...
        'setFilterParams: hi-cut freq must be < smpRate/2 (Nyquist freq = %d) or Inf', smpRate/2);
  assert(any(strcmp(name,{'butter','ellip','fir'})), ...
        sprintf('setFilterParams: Unknown filter name ''%s''. Supported filters: ''butter'',''ellip'',''fir''', name));
  if strcmp(name,'fir')
    assert(strcmp(form,'b,a'), ...
        sprintf('setFilterParams: Must set form = ''b,a'' for fir filter (''%s'' input)', form));
  end
      
  fNyquist      = smpRate/2;            % Nyquist frequency = 1/2 sampling rate
  freqBandsNorm = freqBands./fNyquist;  % Normalized frequency range (by Nyquist freq)

  if nFreqBands > 1
    switch form
    case 'b,a';
      b = cell(nFreqBands,1);
      a = cell(nFreqBands,1);    
    case 'z,p,k';
      z = cell(nFreqBands,1);
      p = cell(nFreqBands,1);    
      k = cell(nFreqBands,1);    
    case 'sos,g';
      sos = cell(nFreqBands,1);
      g = cell(nFreqBands,1);      
    end      
  end
  
  %% Step thru each frequency band requested
  for iF = 1:nFreqBands
    % Low-pass filtering -- f(1) == 0 or -Inf
    if (freqBandsNorm(iF,1) == 0) || (freqBandsNorm(iF,1) == -Inf)
      fBandTmp    = freqBandsNorm(iF,2);
      if isTypeGiven, typeTmp = type{iF};
      else,           typeTmp = 'low';
      end
      
    % High-pass filtering -- f(2) == fNyquist or +Inf
    elseif (freqBandsNorm(iF,2) == 1) || (freqBandsNorm(iF,2) == Inf)  
      fBandTmp    = freqBandsNorm(iF,1);
      if isTypeGiven, typeTmp = type{iF};
      else,           typeTmp = 'high';
      end
      
    % Band-pass (default for rest of cases) or band-stop  filtering  
    else
      fBandTmp    = freqBandsNorm(iF,:);
      if isTypeGiven, typeTmp = type{iF};
      else,           typeTmp = 'bandpass';
      end      
    end
    
    if strcmp(typeTmp,'bandpass'),  extraArgs = {};
    else,                           extraArgs = {typeTmp};
    end

    %% Set filter coefficients based on above parameters
    switch lower(name)
    % Butterworth filter
    case 'butter';
      if strcmp(form,'b,a')
        [bTmp,aTmp] = butter(order, fBandTmp(:), extraArgs{:});
      else
        [zTmp,pTmp,kTmp] = butter(order, fBandTmp(:), extraArgs{:});
      end

    % Elliptical filter
    case 'ellip';
      if strcmp(form,'b,a')
        [bTmp,aTmp] = ellip(order, 3, 50, fBandTmp(:), extraArgs{:});
      else
        [zTmp,pTmp,kTmp] = ellip(order, 3, 50, fBandTmp(:), extraArgs{:});
      end
      
    % Finite impulse response filter  
    case 'fir';
      [bTmp,aTmp] = fir1(order, fBandTmp(:), extraArgs{:});
    end
    
    % Convert zero-pole-gain filter parameters to second-order sections form
    if strcmp(form,'sos,g')
      [sosTmp,gTmp] = zp2sos(zTmp,pTmp,kTmp);
    end
    
    % Save to output variables
    switch form
    case 'b,a';
      if nFreqBands > 1,  b{iF} = bTmp;   a{iF} = aTmp;
      else,               b     = bTmp;   a     = aTmp;      
      end
    case 'z,p,k';
      if nFreqBands > 1,  z{iF} = zTmp;   p{iF} = pTmp;   k{iF} = kTmp;
      else,               z     = zTmp;   p     = pTmp;   k     = kTmp; 
      end
    case 'sos,g';
      if nFreqBands > 1,  sos{iF} = sosTmp;   g{iF} = gTmp;
      else,               sos     = sosTmp;   g     = gTmp;      
      end
    end
  end
  
  % Copy to varargout
  switch form
  case 'b,a';     varargout = {b,a};
  case 'z,p,k';   varargout = {z,p,k};
  case 'sos,g';   varargout = {sos,g};
  end  
end


