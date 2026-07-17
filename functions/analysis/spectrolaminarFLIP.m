function flip = spectrolaminarFLIP(FT_data, opt)
% spectrolaminarFLIP  vFLIP spectrolaminar mapping wrapper.
%
% PURPOSE:
%   Wraps vFLIP_NGL (functions/video/vFLIP_NGL.m) as a first-class LFP
%   analysis with the same provenance / per-area conventions the rest
%   of NGL07 uses. vFLIP identifies the superficial / deep cortical
%   layers of a linear probe from the crossover of low-freq and high-
%   freq power along the shank (Mendoza-Halliday et al. 2024).
%
% USAGE:
%   flip = spectrolaminarFLIP(FT_data, opt);
%
% INPUTS:
%   FT_data - continuous FieldTrip struct with .label, .trial{1},
%             .time{1}, .fsample. When chanArea is present and
%             opt.lfp.tfrAreaFilter is set, the analysis is run on the
%             restricted channel subset (typical: one area = one probe).
%   opt     - resolved options. Consumes:
%               .lfp.flip.laminaraxis   default 0:0.05:1.55
%               .lfp.flip.freqaxis      default 1:150
%               .lfp.flip.setfreqbool   0 = vFLIP (default), 1 = default FLIP
%               .lfp.tfrAreaFilter      optional per-area subselection
%
% OUTPUT (struct):
%   .FLIP        - per-probe FLIP struct (fields from vFLIP_NGL:
%                    startinglowfreq/endinglowfreq/startinghighfreq/
%                    endinghighfreq/goodnessvalue/superficialchannel/
%                    deepchannel/highfreqmaxchannel/lowfreqmaxchannel/
%                    crossoverchannel/laminaraxis/freqaxis).
%   .relpow      - per-probe relative power maps (freq x channel).
%   .abspow      - per-probe absolute power maps.
%   .provenance  - buildLFPProvenance snapshot.
%
% NOTES:
%   * vFLIP_NGL groups channels into "banks" by the FIRST CHARACTER of
%     the FieldTrip label ('cellfun(@(x) x(1), FT_data.hdr.label)').
%     For INTAN-native labels ('A-000', 'B-001', ...) that works out to
%     "one bank per headstage". For projects where that heuristic breaks,
%     re-label channels or extend vFLIP_NGL to accept an explicit bank
%     mapping.
%
% SEE ALSO:
%   vFLIP_NGL, FLIPAnalysis (toolboxes/vFLIP), NGL07_LFPanalysis.
%
% Last modified 26.06.2026 (Jesus) - new wrapper (LFP Pass 3).

    laminaraxis  = localOptField(opt, {'lfp','flip','laminaraxis'},   0:0.05:1.55);
    freqaxis     = localOptField(opt, {'lfp','flip','freqaxis'},      1:150);
    setfreqbool  = localOptField(opt, {'lfp','flip','setfreqbool'},   0);

    %% Optional per-area subselection.
    keepChanIdx = 1:numel(FT_data.label);
    if isfield(opt,'lfp') && isfield(opt.lfp,'tfrAreaFilter') ...
            && ~isempty(opt.lfp.tfrAreaFilter) && isfield(FT_data,'chanArea')
        keepChanIdx = find(ismember(FT_data.chanArea, cellstr(opt.lfp.tfrAreaFilter)));
        if isempty(keepChanIdx)
            error('NGL:spectrolaminarFLIP:noChans', ...
                'opt.lfp.tfrAreaFilter matched 0 channels.');
        end
    end
    if numel(keepChanIdx) ~= numel(FT_data.label)
        selCfg = []; selCfg.channel = FT_data.label(keepChanIdx);
        FT_data = ft_selectdata(selCfg, FT_data);
    end

    if ~isfield(FT_data, 'hdr') || ~isfield(FT_data.hdr, 'label')
        FT_data.hdr = struct('label', {FT_data.label});
    end

    [FLIP, relpow, abspow] = vFLIP_NGL(FT_data, laminaraxis, freqaxis, setfreqbool);

    flip = struct();
    flip.FLIP       = FLIP;
    flip.relpow     = relpow;
    flip.abspow     = abspow;
    flip.provenance = buildLFPProvenance('', opt, struct( ...
        'analysis',    'spectrolaminarFLIP', ...
        'nChans',      numel(FT_data.label), ...
        'setfreqbool', setfreqbool, ...
        'laminaraxis', laminaraxis, ...
        'freqaxis',    [min(freqaxis) max(freqaxis)]));
end


% =======================================================================
function v = localOptField(opt, path, dflt)
    v = dflt;  cursor = opt;
    for k = 1:numel(path)
        if isstruct(cursor) && isfield(cursor, path{k})
            cursor = cursor.(path{k});
        else
            return
        end
    end
    v = cursor;
end
