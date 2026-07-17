function [TFR, cfg] = computeTrialparsedTFR(FT_data, condition, param, opt, alignName)
% computeTrialparsedTFR  Generic project-agnostic trial-parsed TFR core.
%
% PURPOSE:
%   The project-neutral half of the pre-Pass-2 trialparsed_MTspectrogram.
%   Runs ft_freqanalysis on the whole trial pool for one alignment, with
%   cfg.keeptrials = 'yes' so downstream code (compareByBlock, planned
%   NGL07_LFPanalysis regressions and spike-field analyses) has per-trial
%   power to work with. No hardcoded conditions, no block halving, no
%   NS-FS subtraction - those live in projects/socialLearning/
%   computeTrialparsedTFR_ASL.m and are only invoked when the project
%   flag opt.proj_socialLearning is on.
%
% USAGE:
%   [TFR, cfg] = computeTrialparsedTFR(FT_data, condition, param, opt);
%   [TFR, cfg] = computeTrialparsedTFR(FT_data, condition, param, opt, alignName);
%
% INPUTS:
%   FT_data   - FieldTrip trial-parsed data for ONE alignment. May carry
%               .chanArea (Pass 2 tagging); if present, opt.tfr.areaFilter
%               restricts the analysis to a subset of channels via
%               ft_selectdata before ft_freqanalysis.
%   condition - condition struct (from NGL02_postPhy). Only used for
%               provenance and downstream contrast calls; not read here.
%   param     - unused here (kept in the signature for symmetry with
%               the ASL wrapper).
%   opt       - resolved options. Consumes:
%                 .freqInterest      cell of freq vectors, one per band
%                 .TFRmethod         'wavelet'|'mtmconvol'|'superlet'
%                 .timeResol         TFR toi step (s)
%                 .superletOrder     cell (per band; when method=superlet)
%                 .width             cell (per band; when method=superlet)
%                 .combine           char  (when method=superlet)
%                 .analysis          output folder root
%                 .SavFileName       session-name prefix
%                 .lfp.tfrCacheDir   (optional) explicit cache folder
%                 .lfp.tfrAreaFilter (optional) char/cellstr of areas
%                                     to keep before ft_freqanalysis
%   alignName - alignment tag; baked into the cache key and save filename.
%
% OUTPUT:
%   TFR - cell {nBands, 1} of ft_freqanalysis output structs
%         (cfg.keeptrials = 'yes'; TFR{f}.powspctrm is
%         [Ntrials x Nchan x Nfreq x Ntime]).
%   cfg - cell {nBands, 1} of the cfg structs used per band.
%
% CACHE:
%   tfrCacheKey / loadTFRcache / saveTFRcache. Cache key encodes
%   (alignment, method, freqSignature, area-filter tag). Staleness by
%   mtime of the FT source file, matching the spike-side idiom.
%
% PROVENANCE:
%   Every saved .mat carries a `provenance` struct from
%   buildLFPProvenance covering source_FTfile / mtime / opt snapshot /
%   MATLAB / FT version / caller / timestamp.
%
% SEE ALSO:
%   trialparsed_MTspectrogram (thin shim that dispatches based on
%   opt.proj_socialLearning), computeTrialparsedTFR_ASL (the ASL
%   project wrapper that uses this core as scaffolding).
%
% Last modified 26.06.2026 (Jesus)

    if nargin < 5 || isempty(alignName), alignName = ''; end
    assert(isstruct(opt) && isfield(opt, 'freqInterest'), ...
        'NGL:computeTrialparsedTFR:noFreqInterest', ...
        'opt.freqInterest is required (cell of freq vectors, one per band).');
    assert(iscell(opt.freqInterest) && ~isempty(opt.freqInterest), ...
        'NGL:computeTrialparsedTFR:emptyFreqInterest', ...
        'opt.freqInterest must be a non-empty cell of numeric vectors.');

    %% Optional per-area subselection (Q3 answer, LFP Pass 2).
    areaTag = 'all';
    if isfield(opt, 'lfp') && isfield(opt.lfp, 'tfrAreaFilter') ...
            && ~isempty(opt.lfp.tfrAreaFilter) && isfield(FT_data, 'chanArea')
        keepAreas = cellstr(opt.lfp.tfrAreaFilter);
        keepMask  = ismember(FT_data.chanArea, keepAreas);
        if ~any(keepMask)
            error('NGL:computeTrialparsedTFR:noChans', ...
                'opt.lfp.tfrAreaFilter = {%s} matched 0 of %d channels (chanArea unique = {%s}).', ...
                strjoin(keepAreas, ', '), numel(FT_data.chanArea), ...
                strjoin(unique(FT_data.chanArea), ', '));
        end
        selCfg = []; selCfg.channel = FT_data.label(keepMask);
        FT_data = ft_selectdata(selCfg, FT_data);
        areaTag = strjoin(sort(keepAreas), '+');
    end

    %% Resolve the auto knobs (defaulting for older opt structs).
    autoDs      = localOptField(opt, {'lfp','autoDownsample'},       true);
    autoDsFac   = localOptField(opt, {'lfp','autoDownsampleFactor'}, 4);
    autoMet     = localOptField(opt, {'lfp','autoMethod'},           true);
    autoMetHz   = localOptField(opt, {'lfp','autoMethodThresholdHz'},30);
    autoFoi     = localOptField(opt, {'lfp','autoFoi'},              true);
    autoFoiStep = localOptField(opt, {'lfp','autoFoiStep'},          1/4);
    explicitDs  = localOptField(opt, {'lfp','downsampleFs'},         []);

    %% Compute the per-band effective (foi, method, dsFs) UP FRONT.
    % Deterministic function of opt + input bands; used both for logging
    % and for the cache key so a change in any auto knob invalidates.
    nBands   = numel(opt.freqInterest);
    effFoi   = cell(nBands, 1);
    effMeth  = cell(nBands, 1);
    effDsFs  = zeros(nBands, 1);
    for f = 1:nBands
        raw     = opt.freqInterest{f};
        rmin    = min(raw);  rmax = max(raw);
        % (3) Auto foi: quarter-octave log-spaced across [rmin rmax].
        if autoFoi
            steps    = 0:autoFoiStep:(log2(rmax/rmin));
            effFoi{f} = rmin * 2.^steps;
            % Ensure the last point sits at rmax to guarantee band coverage.
            if effFoi{f}(end) < rmax * 0.999
                effFoi{f}(end + 1) = rmax;
            end
        else
            effFoi{f} = raw;
        end
        % (2) Auto method: mtmconvol+hanning below the threshold.
        if autoMet && rmax <= autoMetHz
            effMeth{f} = 'mtmconvol';
        else
            effMeth{f} = opt.TFRmethod;
        end
        % (1) Auto downsample: target = autoDsFac * rmax per band.
        if ~isempty(explicitDs)
            effDsFs(f) = explicitDs;                 % user override
        elseif autoDs
            effDsFs(f) = ceil(autoDsFac * rmax);     % per-band Nyquist*factor
        else
            effDsFs(f) = 0;                          % no downsample
        end
    end

    %% Cache key + fast path.
    % Auto flags + per-band effective settings baked into extraKey so
    % flipping any auto knob (or changing the factor) invalidates cleanly.
    extraKey = struct();
    if isfield(opt,'timeResol'), extraKey.timeRes = opt.timeResol; end
    if isfield(opt,'combine'),   extraKey.combine = opt.combine;   end
    extraKey.autoDs   = autoDs;
    extraKey.autoDsFc = autoDsFac;
    extraKey.autoMet  = autoMet;
    extraKey.autoMetHz = autoMetHz;
    extraKey.autoFoi  = autoFoi;
    extraKey.autoStep = autoFoiStep;
    extraKey.toiRange = opt.toi;
    key = tfrCacheKey(alignName, opt.TFRmethod, opt.freqInterest, areaTag, extraKey);

    cacheDir = localCacheDir(opt);
    sourceFT = localSourceFT(opt, alignName);
    [TFR, cfg, hit] = loadTFRcache(cacheDir, key, sourceFT);
    if hit
        fprintf('computeTrialparsedTFR: cache hit (%s).\n', key);
        return
    end

    %% Announce the auto choices per band so the user sees what's on.
    for f = 1:nBands
        raw = opt.freqInterest{f};
        fprintf(['computeTrialparsedTFR: band %d ', ...
                 'raw=[%g:%g:%g]Hz -> foi=[%g..%g] %d pts | method=%s | dsFs=%s\n'], ...
                f, min(raw), (max(raw)-min(raw))/max(numel(raw)-1,1), max(raw), ...
                min(effFoi{f}), max(effFoi{f}), numel(effFoi{f}), ...
                effMeth{f}, ...
                localFmtDsFs(effDsFs(f), FT_data.fsample));
    end

    %% Parallel mode setup.
    global ft_default;
    ft_default.notification.warning = [];

    parMode = 'none';
    if isfield(opt,'lfp') && isfield(opt.lfp,'parallel'), parMode = lower(opt.lfp.parallel); end
    poolobj = [];
    if strcmp(parMode, 'trials')
        try
            poolobj = gcp('nocreate');
            if isempty(poolobj), poolobj = parpool; end
        catch ME
            warning('NGL:computeTrialparsedTFR:noPool', ...
                'Requested opt.lfp.parallel=''trials'' but could not open parpool (%s). Falling back to serial.', ME.message);
            parMode = 'none';
        end
    end

    %% Pre-build per-band cfgs (each with its own effective foi + method).
    toiVec = opt.toi(1):opt.timeResol:opt.toi(2);
    cfgs   = cell(nBands, 1);
    for f = 1:nBands
        c            = struct();
        c.method     = effMeth{f};
        c.output     = 'pow';
        c.pad        = 'nextpow2';
        c.keeptrials = 'yes';
        c.foi        = effFoi{f};

        switch lower(effMeth{f})
            case 'wavelet'
                c.width = 7;
                c.toi   = toiVec;
            case 'mtmconvol'
                c.taper     = 'hanning';
                c.t_ftimwin = 3 ./ c.foi;
                c.toi       = toiVec;
                c.tapsmofrq = 4;
            case 'superlet'
                assert(isfield(opt,'width') && numel(opt.width) >= f, ...
                    'NGL:computeTrialparsedTFR:noWidth', ...
                    'opt.width{%d} is required for superlet method.', f);
                assert(isfield(opt,'superletOrder') && numel(opt.superletOrder) >= f, ...
                    'NGL:computeTrialparsedTFR:noOrder', ...
                    'opt.superletOrder{%d} is required for superlet method.', f);
                c.toi     = toiVec;
                c.width   = opt.width{f};
                c.combine = opt.combine;
                c.order   = calculate_superlet_order(c.foi, opt.superletOrder{f});
        end

        if strcmp(parMode, 'trials') && ~isempty(poolobj)
            c.parallel = poolobj;
        end

        cfgs{f} = c;
    end

    %% Per-band downsample-then-freqanalysis. Serial or parfor over bands.
    % Downsample happens INSIDE the per-band block (not once up front) so
    % each band picks its own target sample rate. The downsampled FT_data
    % is band-local and doesn't leak.
    TFR = cell(nBands, 1);
    if strcmp(parMode, 'bands') && nBands >= 2
        fprintf('computeTrialparsedTFR: parfor over %d bands.\n', nBands);
        parfor f = 1:nBands
            FT_band = localMaybeDownsample(FT_data, effDsFs(f));
            TFR{f}  = ft_freqanalysis(cfgs{f}, FT_band);
        end
    else
        for f = 1:nBands
            FT_band = localMaybeDownsample(FT_data, effDsFs(f));
            TFR{f}  = ft_freqanalysis(cfgs{f}, FT_band);
        end
    end
    cfg = cfgs;

    %% Persist with provenance.
    provenance = buildLFPProvenance(sourceFT, opt, struct( ...
        'area',      areaTag, ...
        'align',     alignName, ...
        'nBands',    nBands, ...
        'method',    opt.TFRmethod));
    saveTFRcache(cacheDir, key, TFR, cfg, provenance);
    fprintf('computeTrialparsedTFR: cache miss -> wrote %s\n', ...
            fullfile(cacheDir, [key '.mat']));
end

%% Local helpers
function d = localCacheDir(opt)
% Prefer opt.lfp.tfrCacheDir if set; else <opt.analysis>/cache/lfp_tfr/.
    if isfield(opt,'lfp') && isfield(opt.lfp,'tfrCacheDir') && ~isempty(opt.lfp.tfrCacheDir)
        d = opt.lfp.tfrCacheDir;
    else
        d = fullfile(opt.analysis, 'cache', 'lfp_tfr');
    end
end

function p = localSourceFT(opt, alignName)
% Best-effort source FT file path for staleness / provenance.
% Trial-parsed FT files live at <trialSorted>/<SavFileName>_<align>.mat.
    if isfield(opt,'trialSorted') && isfield(opt,'SavFileName') && ~isempty(alignName)
        p = fullfile(opt.trialSorted, [opt.SavFileName '_' alignName '.mat']);
    else
        p = '';
    end
end

function FT_out = localMaybeDownsample(FT_in, targetFs)
% Downsample-a-working-copy helper. targetFs = 0 or >= FT_in.fsample -> no-op.
% Uses ft_resampledata (which handles anti-aliasing via a lowpass filter).
    if targetFs <= 0 || targetFs >= FT_in.fsample - 1e-6
        FT_out = FT_in;
        return
    end
    rsCfg           = [];
    rsCfg.resamplefs = targetFs;
    rsCfg.detrend    = 'no';
    rsCfg.demean     = 'no';
    FT_out          = ft_resampledata(rsCfg, FT_in);
end

function s = localFmtDsFs(targetFs, sourceFs)
    if targetFs <= 0 || targetFs >= sourceFs - 1e-6
        s = sprintf('none (%g Hz source)', sourceFs);
    else
        s = sprintf('%g -> %g Hz (%.1fx)', sourceFs, targetFs, sourceFs / targetFs);
    end
end

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
