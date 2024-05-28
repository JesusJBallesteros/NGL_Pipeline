function [sterr] = nanste(dat,flagIn,dimIn)

% function [sterr] = nanste(dat,flagIn,dimIn)
%
% calculates the sandard error of the mean based on not NAN entries
% passed matrices are processed column by column
% function is based on ste.m
%
% Author : Jonas
% Ver    : 1.2     
% Date   : 10/25/2013

% 1.1    : forced nanstd into first dimension
% 1.2    : added optional input arguments

% defaults
flag = 0;
dim  = 1;

if nargin>1
    flag    = flagIn;
    dim     = dimIn;
end

sterr = nanstd(dat,flag,dim) ./ sqrt(sum(~isnan(dat),dim));

% % old
% bins  = size(dat,2);  
% sterr = zeros(bins,1);
% 
% for i=1:bins
%     sterr(i) = nanstd(dat(:,i),0,1) / sqrt(sum(~isnan(dat(:,i))));
% end

