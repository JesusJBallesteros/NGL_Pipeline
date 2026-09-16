function [epochs, info] = nftEpochs(FT_data, window, baseFreq, varargin)
% nftEpochs  Cut a tagging block into epochs holding a whole number of cycles.
%
% PURPOSE:
%   Frequency tagging lives or dies on where the FFT bins fall. An epoch
%   holding a whole number of stimulation cycles puts the tagging frequency
%   exactly on a bin, and its energy stays there; a fractional epoch smears it
%   across neighbours, which both lowers the peak and raises the baseline the
%   peak is measured against - the two errors push the response down together.
%
%   So epochs are trimmed to an integer cycle count rather than to a round
%   number of seconds, and what that cost in samples is reported.
%
% USAGE:
%   [epochs, info] = nftEpochs(FT_data, [10 190], 1.3)
%   [epochs, info] = nftEpochs(FT_data, win, 1.3, 'epochSeconds', 20)
%   [epochs, info] = nftEpochs(FT_data, win, 1.3, 'mode', 'whole')
%
% INPUTS:
%   FT_data  - continuous FieldTrip data (one trial spanning the recording).
%   window   - [t0 t1] of the block, in seconds on FT_data.time{1}.
%   baseFreq - stimulation frequency of this block, Hz (e.g. 1.3).
%   Name/value pairs:
%     'epochSeconds' target epoch length; trimmed down to whole cycles
%                    (default 20). Longer = finer frequency resolution,
%                    fewer epochs for phase consistency.
%     'mode'         'split' (default) or 'whole' - one epoch per block,
%                    the finest resolution the block allows, but then there
%                    is nothing to measure phase consistency across.
%     'overlap'      fraction of epoch overlap, 0 (default) to < 1. Overlapping
%                    epochs are not independent; use only to stabilise a
%                    spectrum, never to inflate a count entering a test.
%
% OUTPUT:
%   epochs - FieldTrip raw structure whose trials are the epochs.
%   info   - .baseFreq .cyclesPerEpoch .epochSeconds .nEpochs .resolution
%            (Hz per FFT bin) .droppedSeconds .binOfBase (which bin the tag
%            lands on) .window
%
% NOTES:
%   * resolution = 1/epochSeconds. The tag sits on bin cyclesPerEpoch+1 of the
%     one-sided spectrum, so cyclesPerEpoch is also how many bins separate the
%     tag from DC - worth knowing when choosing the neighbourhood in
%     taggingResponse.
%   * An epoch shorter than ~5 cycles makes the neighbourhood too coarse to
%     estimate a baseline; the function warns below that.
%
% Last modified 16.09.2026 (Jesus) - new (NFT project).

    p = inputParser;
    p.addParameter('epochSeconds', 20);
    p.addParameter('mode', 'split');
    p.addParameter('overlap', 0);
    p.parse(varargin{:});
    a = p.Results;

    assert(isstruct(FT_data) && isfield(FT_data, 'trial') && ~isempty(FT_data.trial), ...
        'nftEpochs:input', 'FT_data must be a FieldTrip raw structure.');
    assert(isnumeric(baseFreq) && isscalar(baseFreq) && baseFreq > 0, ...
        'nftEpochs:baseFreq', 'baseFreq must be a positive scalar (Hz).');
    assert(isnumeric(window) && numel(window) == 2 && window(1) < window(2), ...
        'nftEpochs:window', 'window must be [t0 t1] with t0 < t1.');
    assert(a.overlap >= 0 && a.overlap < 1, 'nftEpochs:overlap', ...
        'overlap must be in [0, 1).');

    fs = localFs(FT_data);
    time = FT_data.time{1};
    inWin = time >= window(1) & time <= window(2);
    assert(any(inWin), 'nftEpochs:emptyWindow', ...
        ['block window [%g %g] s holds no samples of the recording (%g to %g s). ', ...
         'Check the block boundaries against this session.'], ...
        window(1), window(2), time(1), time(end));
    first = find(inWin, 1, 'first');
    last  = find(inWin, 1, 'last');
    blockSamples = last - first + 1;
    blockSeconds = blockSamples / fs;

    % Whole cycles per epoch, then back to samples.
    if strcmpi(a.mode, 'whole')
        cycles = floor(blockSeconds * baseFreq);
    else
        cycles = floor(a.epochSeconds * baseFreq);
    end
    assert(cycles >= 1, 'nftEpochs:tooShort', ...
        ['%g s of epoch holds less than one cycle of %g Hz. Raise epochSeconds ', ...
         'above %.2f s.'], a.epochSeconds, baseFreq, 1 / baseFreq);
    if cycles < 5
        warning('nftEpochs:fewCycles', ...
            ['%d cycle(s) per epoch at %g Hz gives %.3f Hz bins - too coarse for a ', ...
             'neighbourhood baseline. Raise epochSeconds towards %.0f s.'], ...
            cycles, baseFreq, baseFreq / cycles, 10 / baseFreq);
    end
    epochSamples = round(cycles / baseFreq * fs);
    assert(epochSamples <= blockSamples, 'nftEpochs:epochTooLong', ...
        ['one epoch of %d cycles needs %.1f s but the block is %.1f s. Lower ', ...
         'epochSeconds.'], cycles, epochSamples / fs, blockSeconds);

    step = max(1, round(epochSamples * (1 - a.overlap)));
    starts = first : step : (last - epochSamples + 1);
    assert(~isempty(starts), 'nftEpochs:noEpochs', 'no epoch fits in the block.');

    trl = [starts(:), starts(:) + epochSamples - 1, zeros(numel(starts), 1)];
    epochs = ft_redefinetrial(struct('trl', trl), FT_data);

    info = struct();
    info.baseFreq       = baseFreq;
    info.cyclesPerEpoch = cycles;
    info.epochSeconds   = epochSamples / fs;
    info.nEpochs        = numel(starts);
    info.resolution     = fs / epochSamples;
    info.binOfBase      = cycles + 1;          % one-sided spectrum, DC = bin 1
    info.droppedSeconds = blockSeconds - (starts(end) + epochSamples - 1 - first + 1) / fs;
    info.window         = window;
    info.overlap        = a.overlap;
    info.mode           = lower(char(a.mode));
end

function fs = localFs(FT_data)
    if isfield(FT_data, 'fsample') && ~isempty(FT_data.fsample)
        fs = FT_data.fsample;
    else
        t = FT_data.time{1};
        fs = 1 / median(diff(t));
    end
end
