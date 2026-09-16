function outFile = saveLFPresult(payload, kind, input, opt, varargin)
% saveLFPresult  Write an LFP result file with its provenance attached.
%
% PURPOSE:
%   One writer for every LFP analysis, so each output carries the same
%   breadcrumbs: which FT file it came from, the opt that produced it, and how
%   the power was normalised and tested. A file that cannot say how it was
%   made cannot be trusted six months later, and this makes saying it the
%   default rather than an act of discipline.
%
% USAGE:
%   f = saveLFPresult(tfr, 'TFR', input, opt, 'area','NCL', 'align','stim2')
%   f = saveLFPresult(struct('stat',stat), 'TFRcontrast', input, opt, ...
%                     'norm', normInfo, 'stats', statInfo, ...
%                     'tags', {'correct-vs-error'})
%
% INPUTS:
%   payload - struct whose fields become the variables in the .mat.
%   kind    - char, the analysis ('TFR', 'TFRcontrast', 'CSD', 'PAC', ...).
%   input   - NGL input struct (recorded for the session it belongs to).
%   opt     - resolved options struct.
%   Name/value pairs:
%     'area' / 'align' / 'tags'  naming parts, as lfpResultName
%     'sourceFT'  path to the FT file this derives from (for provenance)
%     'norm'      info struct from normalizeTFR
%     'stats'     info struct from lfpClusterStats
%     'extra'     struct of further provenance fields
%     'folder'    destination override (default opt.analysis)
%
% OUTPUT:
%   outFile - absolute path written. The .mat holds every field of `payload`
%             plus `provenance`.
%
% NOTES:
%   Saved with '-v7.3' so a TFR with trials kept is not silently truncated at
%   the 2 GB limit of the older format.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 0).

    p = inputParser;
    p.addParameter('area',     '');
    p.addParameter('align',    '');
    p.addParameter('tags',     {});
    p.addParameter('sourceFT', '');
    p.addParameter('norm',     struct());
    p.addParameter('stats',    struct());
    p.addParameter('extra',    struct());
    p.addParameter('folder',   '');
    p.parse(varargin{:});
    a = p.Results;

    assert(isstruct(payload) && ~isempty(fieldnames(payload)), 'saveLFPresult:payload', ...
        'payload must be a non-empty struct; its fields become the saved variables.');

    folder = a.folder;
    if isempty(folder)
        folder = localOpt(opt, {'analysis'}, localOpt(opt, {'FolderProcDataMat'}, pwd));
    end
    if ~isfolder(folder), mkdir(folder); end

    extra = a.extra;
    extra.kind = char(kind);
    if ~isempty(a.area),  extra.area  = char(a.area);  end
    if ~isempty(a.align), extra.align = char(a.align); end
    if ~isempty(fieldnames(a.norm)),  extra.normalize = a.norm;  end
    if ~isempty(fieldnames(a.stats)), extra.statistics = a.stats; end
    if isstruct(input) && isfield(input, 'run'), extra.run = input.run; end
    payload.provenance = buildLFPProvenance(a.sourceFT, opt, extra);

    stem = lfpResultName(kind, opt, 'area', a.area, 'align', a.align, 'tags', a.tags);
    outFile = fullfile(folder, [stem '.mat']);
    save(outFile, '-struct', 'payload', '-v7.3');
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
