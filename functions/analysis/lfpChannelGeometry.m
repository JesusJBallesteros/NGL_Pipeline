function shanks = lfpChannelGeometry(FT_data, input, opt, varargin)
% lfpChannelGeometry  Where each LFP channel sits on the probe, grouped by shank.
%
% PURPOSE:
%   Anything spatial - current source density above all - needs to know which
%   contact is where. The channel map already holds it (xcoords, ycoords in
%   micrometres, kcoords per shank), so this reads it, matches it to the FT
%   channels and hands back one entry per shank with its contacts in depth
%   order.
%
% USAGE:
%   shanks = lfpChannelGeometry(FT_data, input, opt)
%   shanks = lfpChannelGeometry(FT_data, input, opt, 'chanMap', 'C:\...\map.mat')
%
% INPUTS:
%   FT_data - FieldTrip data whose .label lists the channels in acquisition
%             order (optionally .chanArea, carried through to the output).
%   input   - NGL input struct; uses input.areaMap.chanMapPath and
%             input.analysisCode to find the map.
%   opt     - options; uses opt.KSchanMapFile as a fallback name.
%   Name/value pairs:
%     'chanMap' explicit path, overriding the search
%
% OUTPUT (struct array, one per shank):
%   .shank     kcoords value
%   .labels    channel labels, ordered superficial -> deep
%   .idx       their row indices into FT_data.label
%   .depth     [n x 1] ycoords in micrometres, ascending
%   .spacing   median contact spacing (µm)
%   .uniform   true when every gap equals the spacing (CSD needs this)
%   .x         xcoords, for the record
%   .area      the area tag of the shank when chanArea is present
%
% NOTES:
%   * Channels are matched to map rows **by position**: the i-th FT channel is
%     the i-th row of the map. That is the same assumption Kilosort makes when
%     it reads the .bin, so a mismatch here means the sort was wrong too - but
%     it is an assumption, and the function refuses when the counts differ
%     rather than pairing the first N and hoping.
%   * ycoords increases with depth in this pipeline's maps (0 at the tip-most
%     contact). The output is sorted ascending and nothing is flipped: which
%     end is the surface is the user's to know.
%   * A shank with non-uniform spacing is returned with .uniform = false.
%     Standard CSD assumes equal spacing; computeCSD refuses without it.
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 2, CSD).

    p = inputParser;
    p.addParameter('chanMap', '');
    p.parse(varargin{:});

    mapPath = localFindMap(p.Results.chanMap, input, opt);
    m = load(mapPath);
    for f = {'xcoords', 'ycoords', 'kcoords'}
        assert(isfield(m, f{1}), 'lfpChannelGeometry:field', ...
            'channel map %s has no %s.', mapPath, f{1});
    end
    x = m.xcoords(:); y = m.ycoords(:); k = m.kcoords(:);

    nLabel = numel(FT_data.label);
    assert(numel(y) == nLabel, 'lfpChannelGeometry:count', ...
        ['the channel map describes %d contacts but the data has %d channels. ', ...
         'They are matched by position, so a mismatch has to be resolved before ', ...
         'anything spatial is computed (a reduced map from reduceChanMap is the ', ...
         'usual cause - point ''chanMap'' at the reduced one).'], numel(y), nLabel);

    shanks = struct('shank', {}, 'labels', {}, 'idx', {}, 'depth', {}, ...
                    'spacing', {}, 'uniform', {}, 'x', {}, 'area', {});
    for s = unique(k)'
        sel = find(k == s);
        [depth, order] = sort(y(sel), 'ascend');
        idx = sel(order);
        gaps = diff(depth);
        spacing = median(gaps);
        entry = struct();
        entry.shank   = s;
        entry.labels  = FT_data.label(idx)';
        entry.idx     = idx(:)';
        entry.depth   = depth(:);
        entry.spacing = spacing;
        entry.uniform = isempty(gaps) || all(abs(gaps - spacing) < 1e-6);
        entry.x       = x(idx);
        entry.area    = localArea(FT_data, idx);
        shanks(end+1) = entry; %#ok<AGROW>
    end
end

% ---------------- helpers ----------------
function mapPath = localFindMap(explicit, input, opt)
    if ~isempty(explicit)
        assert(isfile(explicit), 'lfpChannelGeometry:noMap', ...
            'channel map not found: %s', explicit);
        mapPath = explicit;
        return
    end
    cand = {};
    if isstruct(input) && isfield(input, 'areaMap') && ~isempty(input.areaMap) ...
            && isfield(input.areaMap, 'chanMapPath')
        cand{end+1} = input.areaMap.chanMapPath;
    end
    if isstruct(opt) && isfield(opt, 'KSchanMapFile') && ~isempty(opt.KSchanMapFile)
        if isstruct(input) && isfield(input, 'analysisCode')
            cand{end+1} = fullfile(input.analysisCode, opt.KSchanMapFile);
        end
        cand{end+1} = opt.KSchanMapFile;
    end
    for c = cand
        if ~isempty(c{1}) && isfile(c{1})
            mapPath = c{1};
            return
        end
    end
    error('lfpChannelGeometry:noMap', ...
        ['no channel map found (looked at: %s). Spatial analyses need one: set ', ...
         'opt.KSchanMapFile, or pass ''chanMap'' explicitly. A linear-array ', ...
         'session run without a map has no geometry to work from.'], ...
        strjoin(cand(~cellfun(@isempty, cand)), '; '));
end

function area = localArea(FT_data, idx)
    area = '';
    if ~isfield(FT_data, 'chanArea') || isempty(FT_data.chanArea), return; end
    tags = unique(string(FT_data.chanArea(idx)), 'stable');
    area = char(strjoin(tags, '+'));     % '+' flags a shank spanning two areas
end
