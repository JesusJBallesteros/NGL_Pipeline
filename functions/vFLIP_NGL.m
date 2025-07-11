function [FLIP, relpow, abspow] = vFLIP_NGL(FT_data, laminaraxis, freqaxis, setfreqbool)
% Based on the resources for vFLIP, available online at
% (https://datadryad.org/stash/dataset/doi:10.5061/dryad.9w0vt4bnp)
% Published Nov 07, 2023; Updated Jan 17, 2024 on Dryad.
%
% Adaptation to NGL processing pipeline by Jesus J. Ballesteros. 
% 1st version on Jul 10, 2024. 
%
% Related publication: 'Mendoza-Halliday, D, Major, A. et al. A ubiquitous
% spectrolaminar motif of local field potential power across the primate
% cortex. 2024.'
% 
% INPUTS:
%           FT_data:      FT formatted data
%           laminaraxis:  vector i.e. [0:0.1:3]. length == nchans
%           freqaxis:     vector i.e. [1:250].
%           setfreqbool:  bool 0 = vFLIP, or 1 = default FLIP.
% OUTPUTS:
%
%
% Requires FT toolbox.
% Jesus J Ballesteros 10.07.2024, v1.0

%% Meta info
banks =  cellfun(@(x) x(1), FT_data.hdr.label, UniformOutput=true);
banklbl = unique(banks);
nProbes = size(banklbl,1);

%% Generate the spectrolaminar pattern (a.k.a. relative power map) from LFP data
for p = 1:nProbes
chn = strcmp(banks, string(banklbl(p)));
    cfg = [];
        cfg.channel = FT_data.label(chn);
        cfg.method = 'mtmfft';
        cfg.taper = 'hanning';
        cfg.output = 'pow';
        cfg.keeptrials = 'yes';
        cfg.foi = freqaxis;
        cfg.timwin = [0 10];
        cfg.pad = 'nextpow2';
    
    pow = ft_freqanalysis(cfg, FT_data);
    
    abspow = squeeze(mean(pow.powspctrm));
    relpow.(banks(p)) = abspow ./ max(abspow);
end  
clear pow cfg

%% Find Optimal Frequency Bin Values with vFLIP
for p = 1:nProbes
    % relpow is relative power matrix calculated from raw LFP using FT, above.
    powdata = squeeze(relpow.(banks(p)));
    powdata(any(isnan(powdata), 2), : ) = [];
    
    [startinglowfreq, endinglowfreq, startinghighfreq, ...
     endinghighfreq, goodnessvalue, superficialchannel, ...
     deepchannel, highfreqmaxchannel, lowfreqmaxchannel, ...
     crossoverchannel] = FLIPAnalysis(powdata, laminaraxis, freqaxis, setfreqbool);
    
    % Collect FLIP outputs
    FLIP.(banks(p)).startinglowfreq    = startinglowfreq;
    FLIP.(banks(p)).endinglowfreq      = endinglowfreq;
    FLIP.(banks(p)).startinghighfreq   = startinghighfreq;
    FLIP.(banks(p)).endinghighfreq     = endinghighfreq;
    FLIP.(banks(p)).goodnessvalue      = goodnessvalue;
    FLIP.(banks(p)).superficialchannel = superficialchannel;
    FLIP.(banks(p)).deepchannel        = deepchannel;
    FLIP.(banks(p)).highfreqmaxchannel = highfreqmaxchannel;
    FLIP.(banks(p)).lowfreqmaxchannel  = lowfreqmaxchannel;
    FLIP.(banks(p)).crossoverchannel   = crossoverchannel;
    FLIP.(banks(p)).laminaraxis        = laminaraxis;
    FLIP.(banks(p)).freqaxis           = freqaxis;
end

end