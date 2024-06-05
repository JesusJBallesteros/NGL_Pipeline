function prettify(varargin)
% A few common instructions to prettify basic plots.
% INPUT:    plotops, struct with a number of fields related to figure attributes,
%                    of which the very basic are defaulted here.

%% Default

% if nargin == 0
plotops = struct('xlabel', {'time'}, 'ylabel', {'metric'}, ...
                 'xticks', [],      'yticks', [], ...
                 'xticklabels', {}, 'yticklabels', {});
isitiON =0;
    
if nargin == 1
    plotops = varargin{1};

elseif nargin == 2
    plotops = varargin{1};
    isitiON = varargin{2};
end

gcf;

box("off");
xline(0,'--k');
ylabel(plotops.ylabel);
xlabel(plotops.xlabel);
yticks(plotops.ytick);

% time axis special case for itiON
if isitiON, xticks('auto');        xticklabels('auto');
else,       xticks(plotops.xtick); xticklabels(plotops.xticklabels{1}); end

yticklabels(plotops.yticklabels{1});

if isfield(plotops, 'title')
    title(plotops.title);
    subtitle(plotops.subtitle);
end

end