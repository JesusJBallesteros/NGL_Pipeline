function figFile = plotContinuousTFR(TFR, param, opt)
% plotContinuousTFR  Render the continuous-mode TFR heatmap PNG.
%
% PURPOSE:
%   Plot half of the pre-26.06.2026 continous_MTspectrogram. Consumes a
%   TFR struct produced by computeContinuousTFR (or loaded from
%   <SavFileName>_TFR_continuous.mat) and writes the same PNG the old
%   monolithic function produced. Kept minimal - all the styling
%   decisions (colormap, zlim, tick spacing) are inherited from the old
%   code path.
%
% USAGE:
%   figFile = plotContinuousTFR(TFR, param, opt);
%   figFile = plotContinuousTFR([], param, opt);   % load TFR from disk
%
% INPUTS:
%   TFR    - struct returned by computeContinuousTFR (needs .dB with the
%            ft_freqanalysis output). Empty -> load the .mat from disk.
%   param  - plotting params. Uses .visible, .size ('adaptive' or [w h]),
%            .Resolution.
%   opt    - resolved options; .analysis, .SavFileName.
%
% OUTPUT:
%   figFile - path of the written PNG.
%
% FILE:
%   <opt.analysis>/plots/TFR/Cont_allCh.png
%
% SEE ALSO:
%   computeContinuousTFR (compute half), NGL02_LFP (caller).
%
% Last modified 26.06.2026 (Jesus) - Pass 1 split of continous_MTspectrogram.

    %% Load TFR from disk if the caller didn't hand it in.
    if isempty(TFR)
        inFile = fullfile(opt.analysis, [opt.SavFileName '_TFR_continuous.mat']);
        assert(isfile(inFile), 'NGL:plotContinuousTFR:noTFR', ...
            'TFR .mat not found: %s. Run computeContinuousTFR first.', inFile);
        S = load(inFile, 'TFR');
        TFR = S.TFR;
    end
    assert(isstruct(TFR) && isfield(TFR, 'dB'), 'NGL:plotContinuousTFR:badTFR', ...
        'TFR must be a struct with a .dB freq field.');

    %% Figure size (matches the pre-split behaviour).
    if strcmpi('adaptive', param.size)
        param.screen.size   = get(0, 'ScreenSize');
        param.screen.width  = param.screen.size(3);
        param.screen.height = param.screen.size(4);
    else
        param.screen.width  = param.size(1);
        param.screen.height = param.size(2);
    end

    %% Plot dB scaled TFR.
    cfg2               = [];
    cfg2.colormap      = hot;
    cfg2.xlim          = 'maxmin';
    if isfield(TFR,'cfg') && isfield(TFR.cfg,'foi')
        cfg2.ylim = [0 max(TFR.cfg.foi)];
    else
        cfg2.ylim = [0 40];
    end
    cfg2.zlim          = [0 500];
    cfg2.colorbartext  = 'Norm. Pow. (dB)';
    cfg2.interactive   = 'no';
    cfg2.fontsize      = 12;
    cfg2.title         = 'allCh';

    fig = figure('Visible', param.visible, 'Position', [0 0 700 400]);
    set(fig, 'Position', [0, 0, round(param.screen.width), round(param.screen.height)]);
        ft_singleplotTFR(cfg2, TFR.dB); hold on
        ylim([0 40]);
        xticks(TFR.dB.time(1):600:TFR.dB.time(end));
        xticklabels(TFR.dB.time(1)/60:10:TFR.dB.time(end)/60);
        xlabel('Time (min)'); ylabel('Frequency (Hz)');
        hold off

    if ~isfolder(fullfile(opt.analysis,'plots','TFR'))
        mkdir(fullfile(opt.analysis,'plots','TFR'));
    end
    figFile = fullfile(opt.analysis, 'plots', 'TFR', 'Cont_allCh.png');
    exportgraphics(fig, figFile, 'Resolution', param.Resolution);
    close all hidden
    fprintf('plotContinuousTFR: wrote %s\n', figFile);
end
