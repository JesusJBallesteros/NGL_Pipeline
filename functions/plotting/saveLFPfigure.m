function figFile = saveLFPfigure(fig, kind, opt, varargin)
% saveLFPfigure  Write an LFP figure where the rest of the pipeline expects it.
%
% PURPOSE:
%   Same naming, same resolution, same folder as the .mat the figure came
%   from, so a figure can always be traced back to its data file by name.
%
% USAGE:
%   f = saveLFPfigure(fig, 'TFR', opt, 'area', 'NCL', 'align', 'stim2')
%   f = saveLFPfigure(fig, 'TFRcontrast', opt, 'tags', {'NCL','correct-vs-error'})
%
% INPUTS:
%   fig  - figure handle.
%   kind - char, the analysis ('TFR', 'TFRcontrast', 'CSD', 'PAC', ...).
%   opt  - options struct; uses opt.analysis (or opt.FolderProcDataMat) as the
%          destination, opt.SavFileName as the session stem, and
%          opt.lfp.plot.Resolution.
%   Name/value pairs:
%     'area' / 'align'  appended to the name when given
%     'tags'            cellstr of extra name parts
%     'folder'          destination override
%     'name'            full stem override (skips the naming rule)
%
% OUTPUT:
%   figFile - absolute path written.
%
% NAMING:
%   <SavFileName>_LFP_<kind>[_<area>][_<align>][_<tags...>].png
%   matching saveLFPresult, so figure and data files sort next to each other.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 0).

    p = inputParser;
    p.addParameter('area',   '');
    p.addParameter('align',  '');
    p.addParameter('tags',   {});
    p.addParameter('folder', '');
    p.addParameter('name',   '');
    p.parse(varargin{:});
    a = p.Results;

    folder = a.folder;
    if isempty(folder)
        folder = localOpt(opt, {'analysis'}, localOpt(opt, {'FolderProcDataMat'}, pwd));
    end
    if ~isfolder(folder), mkdir(folder); end

    stem = a.name;
    if isempty(stem)
        stem = lfpResultName(kind, opt, 'area', a.area, 'align', a.align, 'tags', a.tags);
    end
    figFile = fullfile(folder, [stem '.png']);
    exportgraphics(fig, figFile, 'Resolution', localOpt(opt, {'lfp','plot','Resolution'}, 300));
end

function v = localOpt(opt, path, default)
    v = default;
    s = opt;
    for k = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{k}), return; end
        s = s.(path{k});
    end
    if ~isempty(s) || ischar(s), v = s; end
end
