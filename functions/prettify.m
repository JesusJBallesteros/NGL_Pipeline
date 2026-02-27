function prettify(varargin)
% A few common instructions to prettify basic plots.
% INPUT:    plotops, struct with a number of fields related to figure attributes,
%                    of which the very basic are defaulted here.
%           loop variable, integer to differentiate itiOn vs other alignments

% Default
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

% Modify/add attributes
box("off");
ylabel(plotops.ylabel);
xlabel(plotops.xlabel);
yticks(plotops.yticks);

% % for raster alignment
% if ~plotops.ispsh
%     xline(0,'--k');
% end

% X axis
if isitiON == 1 % iti alignment
    xticks('auto'); 
    xticklabels('auto');
    xlim([plotops.xticks(1) plotops.xticks(end)])
else
    xticks(plotops.xticks);      
    xticklabels(plotops.xticklabels);
    if ~strcmp(plotops.xticks, 'auto')    
        xlim([plotops.xticks(1) plotops.xticks(end)])
    end
end

%  % driftmap
%     xlim([plotops.xtick(1) plotops.xtick(end)])
% end        

% Y axis, normally defined
yticklabels(plotops.yticklabels);
ylim([plotops.yticks(1) plotops.yticks(end)])

% Fonts, for all
set(gca, 'FontSize', 10);
set(get(gca, 'XLabel'), 'FontSize', 12);
set(get(gca, 'YLabel'), 'FontSize', 12);

end