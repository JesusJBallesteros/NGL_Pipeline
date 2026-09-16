function csd = computeCSD(FT_data, shank, opt, varargin)
% computeCSD  Current source density down a shank: where current enters and leaves.
%
% PURPOSE:
%   The LFP at one contact is a mixture of everything conducting to it, so a
%   depth profile of voltage says little about where the current actually
%   moved. The second spatial derivative removes what is common to
%   neighbouring contacts and leaves the local sinks and sources - the
%   laminar signature of an input arriving.
%
%       CSD = -sigma * d2(phi)/dz2
%
%   Sign convention: a **sink** (current entering cells, the extracellular
%   negativity under an active input) is NEGATIVE; a **source** is positive.
%
% USAGE:
%   csd = computeCSD(FT_data, shank, opt)
%   csd = computeCSD(FT_data, shank, opt, 'trials', mask, 'smoothPasses', 2)
%
% INPUTS:
%   FT_data - trial-parsed FieldTrip data for ONE alignment (trials in
%             .trial, a common time axis). The average over trials is what
%             gets differentiated, unless 'keeptrials' is set.
%   shank   - one entry from lfpChannelGeometry: .idx .depth .spacing .labels
%   opt     - options; reads opt.lfp.csd.*
%   Name/value pairs:
%     'trials'       logical mask over trials (default: all)
%     'conductivity' tissue conductivity S/m (default 0.3)
%     'smoothPasses' 3-point Hamming smoothing across CONTACTS before
%                    differentiating (default 1; 0 = none)
%     'vaknin'       duplicate the end contacts so CSD is defined at every
%                    depth instead of losing the outermost two (default true)
%     'keeptrials'   return per-trial CSD as well (default false)
%
% OUTPUT (struct):
%   .csd      [nDepth x nTime] in µA/mm³ (µV, µm and S/m folded in)
%   .lfp      [nDepth x nTime] the averaged LFP it came from, µV
%   .depth    [nDepth x 1] µm
%   .time     [1 x nTime] s
%   .trialCSD [nTrial x nDepth x nTime] when 'keeptrials'
%   .info     settings, nTrials, edge flag
%
% NOTES:
%   * Differentiating amplifies noise - the second derivative of white noise
%     is whiter and larger. Hence the smoothing pass across contacts, which is
%     standard practice and also the main knob: too little and the map is
%     speckled, too much and neighbouring sinks merge into one. Change it
%     deliberately and record what you used (it lands in the provenance).
%   * Uniform contact spacing is required and enforced. With an irregular
%     array the second difference weights gaps unequally and the result is not
%     a CSD; an inverse method (iCSD) would be needed instead.
%   * Vaknin's extension assumes the potential is flat beyond the probe. It
%     buys the two edge depths, but those two are an assumption, not a
%     measurement: .info.edgeEstimated marks them.
%   * The CSD is computed on the trial AVERAGE. Trial-by-trial CSD (keeptrials)
%     is far noisier and mainly useful as input to a statistic, not to look at.
%
% Last modified 16.09.2026 (Jesus) - new (LFP analysis Phase 2, CSD).

    p = inputParser;
    p.addParameter('trials',       []);
    p.addParameter('conductivity', localOpt(opt, {'lfp','csd','conductivity'}, 0.3));
    p.addParameter('smoothPasses', localOpt(opt, {'lfp','csd','smoothPasses'}, 1));
    p.addParameter('vaknin',       localOpt(opt, {'lfp','csd','vaknin'}, true));
    p.addParameter('keeptrials',   false);
    p.parse(varargin{:});
    a = p.Results;

    assert(isstruct(shank) && isfield(shank, 'idx') && isfield(shank, 'depth'), ...
        'computeCSD:shank', 'shank must be one entry from lfpChannelGeometry.');
    assert(numel(shank.idx) >= 3, 'computeCSD:tooFewContacts', ...
        ['CSD needs at least 3 contacts on a shank; this one has %d. A second ', ...
         'derivative cannot be taken from fewer.'], numel(shank.idx));
    assert(shank.uniform, 'computeCSD:spacing', ...
        ['shank %g has uneven contact spacing (median %.1f um). The second ', ...
         'difference assumes equal gaps; an inverse method (iCSD) is needed ', ...
         'for irregular arrays.'], shank.shank, shank.spacing);

    nTrials = numel(FT_data.trial);
    mask = a.trials;
    if isempty(mask), mask = true(nTrials, 1); end
    mask = logical(mask(:));
    assert(numel(mask) == nTrials, 'computeCSD:trials', ...
        'the trial mask covers %d trials but the data holds %d.', numel(mask), nTrials);
    assert(any(mask), 'computeCSD:noTrials', 'the trial mask selects no trials.');

    time = FT_data.time{find(mask, 1)};
    nT = numel(time);
    use = find(mask)';
    stack = nan(numel(use), numel(shank.idx), nT);
    for k = 1:numel(use)
        trl = FT_data.trial{use(k)};
        assert(size(trl, 2) == nT, 'computeCSD:ragged', ...
            ['trial %d has %d samples, the first has %d. CSD needs a common ', ...
             'time axis - use the trial-parsed file for one alignment.'], ...
            use(k), size(trl, 2), nT);
        stack(k, :, :) = trl(shank.idx, :);
    end
    lfp = squeeze(mean(stack, 1, 'omitnan'));

    csd = struct();
    csd.csd   = localCSD(lfp, shank.spacing, a);
    csd.lfp   = lfp;
    csd.depth = shank.depth(:);
    csd.time  = time;
    if a.keeptrials
        csd.trialCSD = nan(size(stack));
        for k = 1:size(stack, 1)
            csd.trialCSD(k, :, :) = localCSD(squeeze(stack(k, :, :)), shank.spacing, a);
        end
    end
    csd.info = struct('conductivity', a.conductivity, 'smoothPasses', a.smoothPasses, ...
                      'vaknin', logical(a.vaknin), 'spacing_um', shank.spacing, ...
                      'shank', shank.shank, 'area', shank.area, ...
                      'labels', {shank.labels}, 'nTrials', sum(mask), ...
                      'units', 'uA/mm^3', 'sign', 'negative = sink', ...
                      'edgeEstimated', logical(a.vaknin));
end

% ---------------- helpers ----------------
function C = localCSD(phi, spacing_um, a)
% phi [nDepth x nTime] in µV, spacing in µm -> CSD in µA/mm³.
    for s = 1:a.smoothPasses
        phi = localSmoothContacts(phi);
    end
    if a.vaknin
        phi = [phi(1, :); phi; phi(end, :)];
    end
    d2 = phi(1:end-2, :) - 2 * phi(2:end-1, :) + phi(3:end, :);
    % µV and µm: -sigma * d2[µV] * 1e-6 / (h[µm] * 1e-6)^2 gives A/m³;
    % A/m³ -> µA/mm³ is a further factor 1e-3, leaving 1e3 / h².
    C = -a.conductivity * d2 * 1e3 / (spacing_um^2);
end

function y = localSmoothContacts(x)
% 3-point Hamming across contacts, edges held by repeating the end contact.
    pad = [x(1, :); x; x(end, :)];
    y = 0.25 * pad(1:end-2, :) + 0.5 * pad(2:end-1, :) + 0.25 * pad(3:end, :);
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
