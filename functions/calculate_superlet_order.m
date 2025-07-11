function [Order] = calculate_superlet_order(foi, orderRange)
% As per Moca VV et al. "Time-Frequency super-resolution with superlets"
% the linearly increasing order magnitude A(f) is a function of the central
% frequency (F), with integer values 'round()':
%   Eq.7:  A(f) = Omin + round((Omax - Omin) * (F-Fmin)/(Fmax-Fmin))

% INPUTS
% foi,  array or cell array of central frequency of interest (integers),
%       e.g foi = 1:1:10 or {[1:1:20]; [21:2:80]}. 
%       Default: foi = 1:1:100;
% orderRange,   array or cell array of desired range of Order to use (integers),
%               e.g orderRange = [1 20] or {[1 2]; [2 20]}. 
%               Default: orderRange = [1 20];

%% Defaults
if isempty(foi), foi = 1:1:100; end
if isempty(orderRange), orderRange = [1 20]; end

%% Define the vectors of linearly increasing integer values for A(f):
% For cell, more than one foi ranges
if iscell(foi) && iscell(orderRange)
    AF = cell(size(foi));
    
    if all(size(foi) == size(orderRange))
        for ff = 1:size(foi,1)
            % get max and min Order values
            Omin = min(orderRange{ff});
            Omax = max(orderRange{ff});
        
            % get max and min Frequency values
            Fmin = min(foi{ff});
            Fmax = max(foi{ff});
    
            for f = 1:size(foi{ff},2)
                AF{ff}(f) = Omin + round((Omax - Omin) * ((foi{ff}(f)-Fmin) / (Fmax-Fmin)));
            end
        end
    else
        warning('When cell arrays, both ''foi'' and ''orderRange'' inputs must have the same size.')
    end

% for array, a single foi range and single orderRange
elseif ~iscell(foi) && ~iscell(orderRange)
    AF = zeros(size(foi));

    % get max and min Order values
    Omin = min(orderRange);
    Omax = max(orderRange);

    % get max and min Frequency values
    Fmin = min(foi);
    Fmax = max(foi);

    for f = 1:size(foi,2)
        AF(1,f) = Omin + round((Omax - Omin) * ((foi(f)-Fmin) / (Fmax-Fmin)));
    end

else
    warning('Both ''foi'' and ''orderRange'' inputs must be same, arrays or cell arrays.')
    AF = [];
end

Order = AF;

end