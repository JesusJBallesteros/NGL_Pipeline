%% NGL07_LFPanalysis.
% PURPOSE:
%   Stage 7: consumes NGL01's continuous FT LFP file, NGL02_postPhy's
%   spikes, and (when available) NGL06_videoAnalysis's behaviour, and
%   produces the research-grade LFP outputs that NGL02_LFP's quick-look
%   pass intentionally does not cover:
%
%     (a) Trial-parsed TFR per opt.alignto              opt.lfp.session.tfr
%     (b) Oscillation / burst detection per band         opt.lfp.session.bursts
%     (c) Continuous phase + envelope per band           opt.lfp.session.phase
%     (d) Spike-field coupling (PPC + coherence)         opt.lfp.session.spikeField
%     (e) LFP x behaviour regression (STUB; needs NGL06 sidecar)
%                                                         opt.lfp.session.behReg
%     * (f) Spectrolaminar (vFLIP) mapping                  opt.lfp.session.flip
%       (g) Event-centered power contrasts, cluster-corrected
%                                                          opt.lfp.session.contrast
%
%   Each analysis is opt-gated so users pick what they need per project.
%   Multi-area handling is via FT_data.chanArea + opt.lfp.tfrAreaFilter -
%   the SAME continuous FT file is loaded once and each analysis
%   restricts to its area subset at the ft_selectdata step (Q3 (b)
%   design decision, LFP Pass 2).
%
% DEPENDENCIES (upstream):
%   NGL01_Main             -> <SavFileName>_FTcont.mat with chanArea tags
%   NGL02_postPhy          -> spike.mat / neurons.mat / condition.mat /
%                             trialdef.mat / events.mat  (required for
%                             opt.lfp.session.spikeField).
%   NGL06_videoAnalysis    -> optional gaze.mat sidecar per session
%                             (required for opt.lfp.session.behReg)
%
% USAGE:
%   Do NOT edit this script. Configure through NGL_SetAndRunMe (section
%   9) and invoke it from there.
%
% REQUIRED WORKSPACE VARIABLES (set by NGL_SetAndRunMe -> NGL00_Prep):
%   datadrive, studyname, subjects, dates, opt
%
% OUTPUTS (per session, under <opt.analysis>/<subject>/<session>/):
%   LFP_TFR_<align>.mat           cell of per-band FT freq structs
%   LFP_bursts.mat                per-band burst episode table
%   LFP_phase.mat                 per-band phase + envelope matrices
%   LFP_spikeField.mat            PPC + spike-field coherence cube
%   LFP_behReg_<covariate>.mat    (stub) regression stat map
%   LFP_FLIP.mat                  spectrolaminar map per bank
%
%   All files carry a `provenance` struct via buildLFPProvenance.
%
% Last modified 26.06.2026 (Jesus) - new stage (LFP Pass 3).

%% 00. Standard scaffolding.
NGL00_Prep

% Areas recovery (same fallback pattern as NGL02_LFP / NGL04).
if ~isfield(input,'Areas') || isempty(input.Areas)
    analysisCodePath = fullfile(input.datadrive, input.studyName, 'analysisCode');
    if ~contains(analysisCodePath, ':\') && ~isempty(input.datadrive)
        analysisCodePath = fullfile([input.datadrive(1) ':\'], input.studyName, 'analysisCode');
    end
    try
        masterInfo = loadPreprocInfo(analysisCodePath, 'master');
        if isfield(masterInfo,'Areas') && ~isempty(masterInfo.Areas)
            input.Areas = masterInfo.Areas;
        end
    catch ME
        if ~strcmp(ME.identifier, 'NGL:loadPreprocInfo:notFound')
            warning('NGL07:preflight', 'Could not read master preprocInfo: %s', ME.message);
        end
    end
end

[input, opt] = set_default(input, opt);
input.sessions = findSessions(input);

%% 00b. Master gate.
if ~localGate(opt, {'lfp','session','do'}, false)
    warning('NGL07:nothingToDo', ...
        'opt.lfp.session.do is false; NGL07_LFPanalysis has nothing to do.');
    return
end

% Warn if no per-analysis gate is set, the driver would then run only
% its scaffolding.
anyAnalysis = localGate(opt, {'lfp','session','tfr'},        false) || ...
              localGate(opt, {'lfp','session','bursts'},     false) || ...
              localGate(opt, {'lfp','session','phase'},      false) || ...
              localGate(opt, {'lfp','session','spikeField'}, false) || ...
              localGate(opt, {'lfp','session','behReg'},     false) || ...
              localGate(opt, {'lfp','session','flip'},       false);
if ~anyAnalysis
    warning('NGL07:noAnalysis', ...
        ['opt.lfp.session.do is true but every per-analysis gate is off. ', ...
         'Turn on at least one of opt.lfp.session.{tfr,bursts,phase,spikeField,behReg,flip}.']);
    return
end

%% 01. Loop over (subject, session).
for x = 1:input.nsubjects
    for y = 1:input.sessions(x).nsessions
        input.run = [x y];
        subject   = input.subjects(x).name;
        session   = input.sessions(x).list{y};
        fprintf('\nNGL07_LFPanalysis: ===== %s / %s =====\n', subject, session);

        try
            [input, opt] = prepSession(input, opt);
        catch ME
            warning('NGL07:prepFail', ...
                '[%s/%s] prepSession failed: %s. Skipping.', subject, session, ME.message);
            continue
        end

        ftFile = fullfile(opt.FolderProcDataMat, [opt.SavFileName '_FTcont.mat']);
        if ~isfile(ftFile)
            fprintf('NGL07: no FTcont at %s; skipping session.\n', ftFile);
            continue
        end
        fprintf('NGL07: loading %s\n', ftFile);
        load(ftFile, '-mat', 'FT_data');
        if isfield(FT_data,'FT_data'), FT_data = FT_data.FT_data; end
        FT_data.cfg.continuous = 'yes';
        areaMapForBackfill = [];
        if isfield(input, 'areaMap'), areaMapForBackfill = input.areaMap; end
        FT_data = ensureChanArea(FT_data, areaMapForBackfill);

        % Condition / trialdef / spike loaded on demand (each analysis
        % states its own prerequisites).
        condition = localLoadIfExists(fullfile(opt.trialSorted, 'condition.mat'), 'condition');
        trialdef  = localLoadIfExists(fullfile(opt.trialSorted, 'trialdef.mat'),  'trialdef');

        %% (a) Trial-parsed TFR (delegates to computeTrialparsedTFR).
        if localGate(opt, {'lfp','session','tfr'}, false)
            if ~iscell(opt.alignto) || isempty(opt.alignto)
                warning('NGL07:tfrNoAlignto', ...
                    '[%s/%s] opt.alignto is empty; skipping trial-parsed TFR.', subject, session);
            else
                % opt.lfp.alignSubset restricts the LFP path only; keeps
                % spike-side alignment handling untouched.
                alignsToRun = opt.alignto;
                if isfield(opt,'lfp') && isfield(opt.lfp,'alignSubset') ...
                        && ~isempty(opt.lfp.alignSubset)
                    keep = ismember(opt.alignto, opt.lfp.alignSubset);
                    if ~any(keep)
                        warning('NGL07:emptyAlignSubset', ...
                            ['[%s/%s] opt.lfp.alignSubset did not intersect opt.alignto; ', ...
                             'skipping trial-parsed TFR for this session.'], subject, session);
                        alignsToRun = {};
                    else
                        alignsToRun = opt.alignto(keep);
                    end
                end
                for k = 1:numel(alignsToRun)
                    alignName = alignsToRun{k};
                    ftAlignFile = fullfile(opt.trialSorted, [opt.SavFileName '_' alignName '.mat']);
                    if ~isfile(ftAlignFile)
                        warning('NGL07:noAlignFT', ...
                            '[%s/%s] no trial-parsed FT for alignment ''%s''; skipping.', ...
                            subject, session, alignName);
                        continue
                    end
                    S = load(ftAlignFile, '-mat', 'FT_data');
                    FTal = S.FT_data;
                    if isfield(FTal,'FT_data'), FTal = FTal.FT_data; end
                    FTal = ensureChanArea(FTal, areaMapForBackfill);
                    param = struct();  % kept for signature symmetry
                    try
                        [TFR, cfg] = computeTrialparsedTFR(FTal, condition, param, opt, alignName);
                        out = fullfile(opt.analysis, [opt.SavFileName '_LFP_TFR_' alignName '.mat']);
                        provenance = buildLFPProvenance(ftAlignFile, opt, struct( ...
                            'analysis','NGL07_tfr','align',alignName));
                        save(out, 'TFR', 'cfg', 'provenance', '-v7.3');
                        fprintf('NGL07: wrote %s\n', out);
                    catch ME
                        warning('NGL07:tfrFail', ...
                            '[%s/%s] TFR (%s) failed: %s', subject, session, alignName, ME.message);
                    end
                    clear FTal S TFR cfg provenance
                end
            end
        end

        %% (b) Burst detection.
        if localGate(opt, {'lfp','session','bursts'}, false)
            try
                bursts = detectBursts(FT_data, opt);
                out = fullfile(opt.analysis, [opt.SavFileName '_LFP_bursts.mat']);
                save(out, 'bursts', '-v7.3');
                fprintf('NGL07: wrote %s (%d episodes total)\n', out, height(bursts.episodes));
            catch ME
                warning('NGL07:burstsFail', ...
                    '[%s/%s] detectBursts failed: %s', subject, session, ME.message);
            end
        end

        %% (c) Phase + envelope per band.
        if localGate(opt, {'lfp','session','phase'}, false)
            try
                bp = hilbertBandpass(FT_data, opt);
                out = fullfile(opt.analysis, [opt.SavFileName '_LFP_phase.mat']);
                save(out, 'bp', '-v7.3');
                fprintf('NGL07: wrote %s\n', out);
            catch ME
                warning('NGL07:phaseFail', ...
                    '[%s/%s] hilbertBandpass failed: %s', subject, session, ME.message);
            end
        end

        %% (d) Spike-field coupling.
        if localGate(opt, {'lfp','session','spikeField'}, false)
            spikeFile = fullfile(opt.spikeSorted, 'spike.mat');
            if ~isfile(spikeFile)
                warning('NGL07:noSpike', ...
                    '[%s/%s] spike.mat missing at %s; skipping spike-field.', ...
                    subject, session, spikeFile);
            else
                Sspk = load(spikeFile, 'spike');
                try
                    sfc = spikeFieldCoupling(FT_data, Sspk.spike, opt); 
                    out = fullfile(opt.analysis, [opt.SavFileName '_LFP_spikeField.mat']);
                    save(out, 'sfc', '-v7.3');
                    fprintf('NGL07: wrote %s\n', out);
                catch ME
                    warning('NGL07:spikeFieldFail', ...
                        '[%s/%s] spikeFieldCoupling failed: %s', subject, session, ME.message);
                end
                clear Sspk sfc
            end
        end

        %% (e) LFP x behaviour regression (STUB).
        if localGate(opt, {'lfp','session','behReg'}, false)
            % Currently a stub - see lfpBehaviorRegression.m for the
            % upstream contract (NGL06 gaze.mat sidecar with per-frame
            % x/y/head_dir/gaze_az).
            try
                reg = lfpBehaviorRegression([], [], opt); 
                % Not saved on stub path; when the stub is filled in,
                % each covariate will be persisted as
                %   LFP_behReg_<covariate>.mat
            catch ME
                warning('NGL07:behRegFail', ...
                    '[%s/%s] lfpBehaviorRegression failed: %s', subject, session, ME.message);
            end
        end

        %% (g) Event-centered power: condition contrasts, cluster-corrected.
        % Reuses the TFR from (a) through computeTrialparsedTFR's own cache,
        % so turning this on without (a) costs one TFR, not two.
        if localGate(opt, {'lfp','session','contrast'}, false)
            pairs = localGate(opt, {'lfp','contrast','pairs'}, {});
            if isempty(pairs)
                warning('NGL07:noContrasts', ...
                    ['[%s/%s] opt.lfp.session.contrast is on but ', ...
                     'opt.lfp.contrast.pairs is empty; nothing to compare.'], ...
                    subject, session);
            elseif isempty(condition)
                warning('NGL07:noCondition', ...
                    ['[%s/%s] condition.mat not found in %s; contrasts need it ', ...
                     'to know which trials are which.'], subject, session, opt.trialSorted);
            else
                localRunContrasts(input, opt, condition, areaMapForBackfill, ...
                                  subject, session);
            end
        end

        %% (f) vFLIP spectrolaminar mapping.
        if localGate(opt, {'lfp','session','flip'}, false)
            try
                flip = spectrolaminarFLIP(FT_data, opt);
                out = fullfile(opt.analysis, [opt.SavFileName '_LFP_FLIP.mat']);
                save(out, 'flip', '-v7.3');
                fprintf('NGL07: wrote %s\n', out);
            catch ME
                warning('NGL07:flipFail', ...
                    '[%s/%s] spectrolaminarFLIP failed: %s', subject, session, ME.message);
            end
        end

        clear FT_data condition trialdef bursts bp flip
    end
end

%% Local helpers
function v = localGate(opt, path, dflt)
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

function localRunContrasts(input, opt, condition, areaMapForBackfill, subject, session)
% Every (alignment x area x contrast) for this session, each written as one
% .mat and one figure. Failures are per-combination: a contrast naming a
% missing condition field must not cost the others.
    aligns = opt.alignto;
    subset = localGate(opt, {'lfp','alignSubset'}, {});
    if ~isempty(subset), aligns = aligns(ismember(aligns, subset)); end
    pairs  = localGate(opt, {'lfp','contrast','pairs'}, {});
    if ischar(pairs), pairs = {pairs}; end
    wantAreas = localGate(opt, {'lfp','contrast','areas'}, {});
    if ischar(wantAreas) && ~isempty(wantAreas), wantAreas = {wantAreas}; end
    doPlot = localGate(opt, {'lfp','contrast','plot'}, true);

    for k = 1:numel(aligns)
        alignName = aligns{k};
        ftAlignFile = fullfile(opt.trialSorted, [opt.SavFileName '_' alignName '.mat']);
        if ~isfile(ftAlignFile)
            warning('NGL07:contrastNoFT', ...
                '[%s/%s] no trial-parsed FT for ''%s''; skipping its contrasts.', ...
                subject, session, alignName);
            continue
        end
        S = load(ftAlignFile, '-mat', 'FT_data');
        FTal = S.FT_data;
        if isfield(FTal, 'FT_data'), FTal = FTal.FT_data; end
        FTal = ensureChanArea(FTal, areaMapForBackfill);

        % chanArea lives on the FT data; ft_freqanalysis does not carry it
        % through, so the label lists are taken here and used to select
        % channels in the TFR.
        areas = localAreaMap(FTal, wantAreas);
        try
            TFR = computeTrialparsedTFR(FTal, condition, struct(), opt, alignName);
        catch ME
            warning('NGL07:contrastTFRfail', ...
                '[%s/%s] TFR for ''%s'' failed: %s', subject, session, alignName, ME.message);
            continue
        end

        for p = 1:numel(pairs)
            for aI = 1:numel(areas)
                areaName = areas(aI).name;
                try
                    res = cell(numel(TFR), 1);
                    for b = 1:numel(TFR)
                        band = TFR{b};
                        if ~isempty(areaName)
                            band = ft_selectdata(struct('channel', {areas(aI).labels}), band);
                        end
                        spec = parseTrialContrast(pairs{p}, condition, size(band.powspctrm, 1));
                        res{b} = computeTFRcontrast(band, spec, opt, ...
                                    'area', areaName, 'align', alignName);
                    end
                    payload = struct('contrast', {res});
                    outFile = saveLFPresult(payload, 'TFRcontrast', input, opt, ...
                        'area', areaName, 'align', alignName, ...
                        'tags', {res{1}.spec.label}, 'sourceFT', ftAlignFile, ...
                        'norm', res{1}.norm, 'stats', res{1}.stats);
                    fprintf('NGL07: wrote %s (%s, %d/%d trials, %d cluster(s))\n', ...
                        outFile, res{1}.spec.request, res{1}.nA, res{1}.nB, ...
                        res{1}.stats.nPos + res{1}.stats.nNeg);
                    if doPlot
                        [fig, figFile] = plotTFRcontrast(res, opt);
                        close(fig);
                        fprintf('NGL07: wrote %s\n', figFile);
                    end
                catch ME
                    warning('NGL07:contrastFail', ...
                        '[%s/%s] contrast ''%s'' (%s, %s) failed: %s', ...
                        subject, session, pairs{p}, areaName, alignName, ME.message);
                end
            end
        end
        clear FTal S TFR
    end
end

function areas = localAreaMap(FTal, wanted)
% One entry per area to analyse, with the channel labels it owns. Empty name
% = every channel together, which is also the fallback when the data carry no
% area tagging.
    areas = struct('name', {}, 'labels', {});
    if ~isfield(FTal, 'chanArea') || isempty(FTal.chanArea)
        areas(1) = struct('name', '', 'labels', {FTal.label(:)'});
        return
    end
    tags = string(FTal.chanArea(:));
    names = unique(tags, 'stable');
    if ~isempty(wanted)
        names = names(ismember(names, string(wanted)));
    end
    for k = 1:numel(names)
        sel = tags == names(k);
        areas(end+1) = struct('name', char(names(k)), ...
                              'labels', {FTal.label(sel)'}); %#ok<AGROW>
    end
    if isempty(areas)
        areas(1) = struct('name', '', 'labels', {FTal.label(:)'});
    end
end

function out = localLoadIfExists(fpath, varName)
    out = struct();
    if ~isfile(fpath), return; end
    S = load(fpath);
    if isfield(S, varName)
        out = S.(varName);
    elseif numel(fieldnames(S)) == 1
        fn = fieldnames(S);
        out = S.(fn{1});
    end
end
