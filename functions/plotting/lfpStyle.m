function st = lfpStyle(opt, kind)
% lfpStyle  Resolve the look of every LFP figure from opt, in one place.
%
% PURPOSE:
%   Colour map, colour limits, fonts and figure size decided once, so a power
%   map, a comodulogram and a CSD from the same session can be read side by
%   side. Anything the user set in opt.lfp.plot.* wins; the rest falls back to
%   defaults chosen per plot kind.
%
% USAGE:
%   st = lfpStyle(opt)                 % sequential (power)
%   st = lfpStyle(opt, 'diverging')    % contrasts, CSD: signed around zero
%
% INPUTS:
%   opt  - options struct. Reads opt.lfp.plot.colormap / .divergingColormap /
%          .zlim / .interp / .visible / .Resolution / .figSize / .fontSize.
%   kind - 'sequential' (default) or 'diverging'.
%
% OUTPUT (struct):
%   .colormap    colormap matrix or name
%   .diverging   logical; true = limits are symmetric about zero
%   .zlim        [] (decide per panel from the data) or a fixed [lo hi]
%   .interp      'bilinear' | 'nearest' | 'none'
%   .visible     'on' | 'off' for new figures
%   .resolution  dpi for exported PNG
%   .figSize     [width height] in pixels
%   .fontSize    base font size
%
% NOTES:
%   A signed quantity on a sequential colour map hides the sign, and zero
%   lands at an arbitrary colour - hence the separate 'diverging' kind, whose
%   limits are symmetric so white/centre always means "no change".
%
% Last modified 16.09.2026 (Jesus) - new helper (LFP analysis Phase 0).

    if nargin < 2 || isempty(kind), kind = 'sequential'; end
    kind = lower(char(kind));
    assert(ismember(kind, {'sequential', 'diverging'}), 'lfpStyle:kind', ...
        'kind must be ''sequential'' or ''diverging'' (got ''%s'').', kind);

    st.diverging  = strcmp(kind, 'diverging');
    if st.diverging
        st.colormap = localOpt(opt, {'lfp','plot','divergingColormap'}, localRdBu());
    else
        st.colormap = localOpt(opt, {'lfp','plot','colormap'}, 'parula');
    end
    st.zlim       = localOpt(opt, {'lfp','plot','zlim'}, []);
    st.interp     = localOpt(opt, {'lfp','plot','interp'}, 'bilinear');
    st.visible    = localOpt(opt, {'lfp','plot','visible'}, 'off');
    st.resolution = localOpt(opt, {'lfp','plot','Resolution'}, 300);
    st.figSize    = localOpt(opt, {'lfp','plot','figSize'}, [1100 700]);
    st.fontSize   = localOpt(opt, {'lfp','plot','fontSize'}, 10);
end

function cm = localRdBu()
% Blue-white-red, so zero is white and the sign is readable at a glance.
% Written out rather than taken from a toolbox: MATLAB ships no diverging map
% before R2023b, and the pipeline must look the same on older releases.
    anchor = [ 0.13 0.30 0.60;
               0.26 0.52 0.78;
               0.65 0.81 0.89;
               1.00 1.00 1.00;
               0.99 0.72 0.60;
               0.84 0.38 0.30;
               0.60 0.10 0.15];
    x  = linspace(0, 1, size(anchor, 1));
    xi = linspace(0, 1, 256);
    cm = [interp1(x, anchor(:,1), xi)', ...
          interp1(x, anchor(:,2), xi)', ...
          interp1(x, anchor(:,3), xi)'];
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
