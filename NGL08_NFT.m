%% NGL08_NFT. Neural frequency tagging (project stage).
%
% PURPOSE:
%   Blocks of periodic stimulation analysed in the frequency domain: is there
%   a response at the frequency the experiment imposed, on which harmonics,
%   and is it locked to the stimulus rather than merely rhythmic?
%
%   Project-specific by design (Q1 convention): the analyses live in
%   functions/analysis/projects/NFT/ and only run behind opt.nft.do.
%
% USAGE:
%   Do NOT run directly. Configure and run from NGL_SetAndRunMe.m, after
%   NGL01 has produced the continuous FieldTrip file.
%
% REQUIRES:
%   <SavFileName>_FTcont.mat   continuous LFP (NGL01)
%   opt.nft.blocks             struct array: .name, .base (Hz), .window [t0 t1]
%
%   Block windows are study-specific: they come from the stimulation log or
%   from event codes, and this stage does not guess them. Fill them in when
%   configuring the study; a session whose blocks are unknown is skipped with
%   a warning rather than analysed against invented boundaries.
%
% PIPELINE (per session x block x area):
%   nftEpochs              cut the block into whole-cycle epochs
%   computeTaggingSpectrum amplitude + complex spectrum over epochs
%   taggingResponse        SNR / z against neighbouring bins, harmonics
%   computeITPC            phase consistency across epochs
%   plotTaggingSpectrum    the figure (opt.nft.plot)
%
% OUTPUTS (per session, in opt.analysis):
%   <SavFileName>_LFP_NFT_<area>_<block>_<base>Hz.mat   and .png
%
% Last modified 16.09.2026 (Jesus) - new stage (NFT project).

%% 00. Scaffolding.
NGL00_Prep
[input, opt] = set_default(input, opt);
input.sessions = findSessions(input);

if ~isfield(opt, 'nft') || ~opt.nft.do
    fprintf('NGL08_NFT: opt.nft.do is false; nothing to do.\n');
    return
end
blocks = opt.nft.blocks;
if isempty(blocks)
    error('NGL08:noBlocks', ...
        ['opt.nft.blocks is empty. Define the tagging blocks as a struct array ', ...
         'with .name, .base (Hz) and .window [t0 t1] in seconds.']);
end
assert(all(isfield(blocks, {'name', 'base', 'window'})), 'NGL08:blockFields', ...
    'each entry of opt.nft.blocks needs .name, .base and .window.');

%% 01. Loop over (subject, session).
for x = 1:input.nsubjects
    for y = 1:input.sessions(x).nsessions
        input.run = [x y];
        subject = input.subjects(x).name;
        session = input.sessions(x).list{y};
        fprintf('\nNGL08_NFT: ===== %s / %s =====\n', subject, session);

        [input, opt] = prepforsession(input, opt);
        ftFile = fullfile(opt.FolderProcDataMat, [opt.SavFileName '_FTcont.mat']);
        if ~isfile(ftFile)
            warning('NGL08:noFT', '[%s/%s] no FTcont at %s; skipping.', ...
                    subject, session, ftFile);
            continue
        end
        S = load(ftFile, '-mat', 'FT_data');
        FT_data = S.FT_data;
        if isfield(FT_data, 'FT_data'), FT_data = FT_data.FT_data; end
        areaMapForBackfill = [];
        if isfield(input, 'areaMap'), areaMapForBackfill = input.areaMap; end
        FT_data = ensureChanArea(FT_data, areaMapForBackfill);

        areas = nftAreaList(FT_data, opt.nft.areas);

        for b = 1:numel(blocks)
            blk = blocks(b);
            fprintf('NGL08: block %s, %g Hz, %g-%g s\n', blk.name, blk.base, ...
                    blk.window(1), blk.window(2));
            try
                [epochs, einfo] = nftEpochs(FT_data, blk.window, blk.base, ...
                    'epochSeconds', opt.nft.epochSeconds, 'overlap', opt.nft.overlap);
            catch ME
                warning('NGL08:epochFail', '[%s/%s] block %s: %s', ...
                        subject, session, blk.name, ME.message);
                continue
            end
            fprintf('NGL08:   %d epochs x %.2f s (%d cycles), %.4f Hz bins\n', ...
                    einfo.nEpochs, einfo.epochSeconds, einfo.cyclesPerEpoch, einfo.resolution);

            for aI = 1:numel(areas)
                areaName = areas(aI).name;
                try
                    ep = ft_selectdata(struct('channel', {areas(aI).labels}), epochs);
                    spec = computeTaggingSpectrum(ep, opt);
                    resp = taggingResponse(spec, blk.base, opt);
                    itpc = computeITPC(spec, resp.harmonics);

                    payload = struct('response', resp, 'itpc', itpc, ...
                                     'epochs', einfo, 'block', blk);
                    outFile = saveLFPresult(payload, 'NFT', input, opt, ...
                        'area', areaName, 'tags', {blk.name, sprintf('%gHz', blk.base)}, ...
                        'sourceFT', ftFile, ...
                        'extra', struct('block', blk, 'epochs', einfo));
                    [best, bestCh] = max(resp.sumAmp);
                    fprintf('NGL08:   %-6s strongest %s, base z = %.1f, ITPC = %.2f, sum = %.3g\n', ...
                        areaName, resp.label{bestCh}, resp.z(bestCh, 1), ...
                        itpc.itpc(bestCh, 1), best);
                    fprintf('NGL08:   wrote %s\n', outFile);

                    if opt.nft.plot
                        [fig, figFile] = plotTaggingSpectrum(resp, opt, 'itpc', itpc, ...
                            'block', blk.name, 'area', areaName);
                        close(fig);
                        fprintf('NGL08:   wrote %s\n', figFile);
                    end
                catch ME
                    warning('NGL08:blockFail', '[%s/%s] block %s, area %s: %s', ...
                            subject, session, blk.name, areaName, ME.message);
                end
            end
        end
        clear FT_data S epochs spec resp itpc
    end
end

%% Local helpers
function areas = nftAreaList(FT_data, wanted)
% One entry per area, with the channels it owns; '' = everything together.
    areas = struct('name', {}, 'labels', {});
    if ~isfield(FT_data, 'chanArea') || isempty(FT_data.chanArea)
        areas(1) = struct('name', '', 'labels', {FT_data.label(:)'});
        return
    end
    tags = string(FT_data.chanArea(:));
    names = unique(tags, 'stable');
    if ~isempty(wanted), names = names(ismember(names, string(wanted))); end
    for k = 1:numel(names)
        areas(end+1) = struct('name', char(names(k)), ...
                              'labels', {FT_data.label(tags == names(k))'}); %#ok<AGROW>
    end
    if isempty(areas)
        areas(1) = struct('name', '', 'labels', {FT_data.label(:)'});
    end
end
