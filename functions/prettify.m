function prettify(varargin)
% A few common instructions to prettify basic plots.
% INPUT:    plotops, struct with a number of fields related to figure attributes,
%                    of which the very basic are defaulted here.

%% Default
if nargin == 0
    plotops = struct('xlabel', {'time'}, 'ylabel', {'metric'}, ...
                     'xticks', [],      'yticks', [], ...
                     'xticklabels', {}, 'yticklabels', {});
    isitiON = 0;   
elseif nargin == 1
    plotops = varargin{1};
    isitiON = 0;   
elseif nargin == 2
    plotops = varargin{1};
    isitiON = varargin{2};
end

% Get Current Figure
gcf;

% Modify attributes
box("off");
xline(0,'--k');
ylabel(plotops.ylabel);
xlabel(plotops.xlabel);
yticks(plotops.ytick);

% X axis, special cases
if isitiON == 1, xticks('auto');        xticklabels('auto');
else,            xticks(plotops.xtick); xticklabels(plotops.xticklabels{1}); end

% Y axis
yticklabels(plotops.yticklabels{1});

% Title
if isfield(plotops, 'title')
    title(plotops.title);
    subtitle(plotops.subtitle);
    set(get(gca, 'Title'), 'FontSize', 18);
end

% Fonts
set(gca, 'FontSize', 8);
set(get(gca, 'XLabel'), 'FontSize', 10);
set(get(gca, 'YLabel'), 'FontSize', 10);

end