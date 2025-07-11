function  dim = firstNonSingletonDim(data)
%FIRSTNONSINGLETONDIM   Find the first non-singleton dimension of an array
%
% function  dim = firstNonSingletonDim(data)
% 
% INPUTS
% data    Array of any arbitary size and type
% 
% OUTPUTS
% dim     Scalar. First non-singleton dimension of data.
%         Returns 1 if no non-singleton dimensions in data

  dim   = find(size(data) ~= 1, 1);
  
  % If no non-singleton dimensions, return 1
  if isempty(dim); dim = 1; end
   
end