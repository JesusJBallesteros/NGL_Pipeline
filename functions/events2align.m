function [t0] = events2align(opt)
% This function is used to adquire equivalences between given event names
% to align to, and their decimal value automatically, following the event
% definitions for the project.
%
% INPUT
% 'opt.alignto': a char array with the name of the event, e.g. 'itiOn'
%          a decimal value that represents an event, e.g. [8]
%          a cell array of strings with several event definitions, e.g. {'itiOn', 'rwd'}
%
% OUTPUT
% 't0': cell array containing the equivalences decimal value and 'event_name'
%       to be used consequently.
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

