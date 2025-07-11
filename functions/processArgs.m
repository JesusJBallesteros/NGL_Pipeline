function [varargout] = processArgs(names,defaults,args,errForMissing)
%PROCESSARGS  Processes function arguments w.r.t default values 
%
% [varargout] = processArgs(names,defaults,args,errForMissing)
% 
% Parses and procesess function arguments given as list of name/value pairs 
% with respect to a list of default values -- overwrites default value
% for any parameters included in input list; otherwise return defaults. 
% 
% As in Matlab built-ins, argument names are case insensitive ('argName'=='ArgName')
% 
% This function merges lists of individual-variable override values and defaults
% In contrast, mergeFields() merges a list of individual-variable override values 
%   with default values contained in a struct or object
% 
% Based on Matlab private function parseArgs.m
% 
% INPUT
% names     {1 x nDefaults cell} of strings. Names of parameters that might be set in arg's
% 
% defaults  {1 x nDefaults cell} of arbitrary types. Default values for each corrsponding param name.
% 
% args      {1 x 2*nArgs cell} of name/value (string/arbitrary type) pairs. 
%           List of arguments used to overwrite default values, 
%           usually obtained from varargin of caller function arguments.
%
% errForMissing Logical. If set=true, an error is raised if any parameters missing 
%           from <names> are encountered [default]. Otherwise, any parameters in <args> 
%           not present in original list of names are just ignored and returned to output.
% 
% OUTPUT
% varargout {1 x nParams cell} of arbitrary types. Set of values for each
%            parameter in names, w/ values modified by name/value pairs in args
%           If errForMissing==false, args name/value pairs that were not matched 
%           are appended to the end of varargout as a cell vector.

  if (nargin < 4) || isempty(errForMissing), errForMissing = true; end  % WAS 12-1: false
  
  % Initialize some variables
  nArgs     = length(args);
  varargout = defaults;               % Save default values to output

  % Check the number of input arguments -- Must be given as name/value pairs
  assert(mod(nArgs,2)==0 && iscellstr(args(1:2:end)), ...
        'processArgs: Options must be given as name/value pairs');  

  % If no error for parameters missing from defaults, keep track of which ones actually matched
  if ~errForMissing, missingParamIdxs = false(1,length(args)); end
  
  % Process name/value pairs
  for iArg = 1:2:nArgs
    name    = args{iArg};        
    idx     = strcmpi(names, name);   % Find argument in default list
       
    % If parameter was found, assign modified value to appropriate element of output
    if any(idx)
      varargout{idx}  = args{iArg+1}; 

    % Otherwise, raise an error or log that it was missing
    else
      if errForMissing
        error('processArgs: Argument ''%s'' not present in default list', name);
      else
        missingParamIdxs(iArg:(iArg+1)) = true;
      end
    end    
  end
  
  % Append unmatched parameters to varargout and return
  if ~errForMissing, varargout{end+1} = args(missingParamIdxs); end
end


