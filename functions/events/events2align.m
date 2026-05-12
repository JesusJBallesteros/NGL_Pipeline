function [t0] = events2align(opt)
% events2align  Convert opt.alignto event names to {name, decimal} pairs.
%
% PURPOSE:
%   Resolves the user-supplied opt.alignto (event name strings) to their
%   decimal integer codes from opt.eventdef, producing the t0 cell array
%   consumed by trialdefGen. Supports a single char array or a cell of chars.
%
% USAGE:
%   t0 = events2align(opt)
%
% INPUT:
%   opt.alignto   - char (single event name, e.g. 'itiOn')
%                   OR cell of chars (e.g. {'itiOn', 'rwd', 'stimOn1'})
%   opt.eventdef  - struct from eventDefinitions; field names are event names,
%                   values are decimal integers
%
% OUTPUT:
%   t0  - (n × 2) or (n × 3) cell array:
%           column 1: event name (char)
%           column 2: decimal event code (double, from opt.eventdef)
%           column 3: qualifier string '0' (or last char of name for bhv variants)
%
% NOTES:
%   - For bhvN-style variants (e.g. 'bhv1', 'bhv2'), the name is truncated to
%     'bhv' for eventdef lookup and the digit becomes column 3.
%   - Output row count equals numel(opt.alignto).
%
% Jesus 30.05.2024

%% Create field 't0' for event/s to align trials to
% this field 't0' will end up being a numerical array i.e. [0 8 12 45 ...] 
    inputclass = class(opt.alignto);
    switch inputclass
        case 'char'
            % single event given as char array
            t0       = {opt.alignto, opt.eventdef.(opt.alignto)};
        case 'cell'
            % multiple events given as a cell array of characters
            t0       = {};
            for i=1:size(opt.alignto,2)
                try t0(i,:)  = {opt.alignto{i}, opt.eventdef.(opt.alignto{i}), '0'};
                catch
                    t0(i,:)  = {opt.alignto{i}(1:3), opt.eventdef.(opt.alignto{i}(1:3)), opt.alignto{i}(end)};
                end
            end
    end
    
% if isfield(opt, 'newEvent')
%     opt.alignto = [opt.alignto opt.newEvent(1)];
% end
% 
end

