function results = detect_jerkEvents(data, fs, varargin)
% DETECT_EVENTS  KS4-inspired two-pass detector for sharp, large-amplitude events in a single-channel continuous 
%   signal, with template learning, PCA-based clustering, and waveform/PC visualisation.
%   Implements the full Kilosort 4 detection logic adapted to 1-D signals:
%       Pass 1 — universal template bank catches all candidate events above a liberal threshold.
%       Learning — PCA + k-means on Pass-1 waveforms learns the actual event shapes present in this signal
%       Pass 2 — learned templates re-detect events at a tighter threshold, assigning each detection to a
%                cluster (event type).
%   USAGE
%       results = detect_events(data, fs)
%       results = detect_events(data, fs, 'Param', Value, ...)
%
%   REQUIRED INPUTS
%       data            [N×1] or [1×N] double — continuous signal
%       fs              scalar — sampling rate in Hz
%
%   OPTIONAL PARAMETERS  (name-value pairs)
%   'CustomTemplates'   [Ntemplates × Nsamples] numeric matrix of waveforms to use as the Pass-1 template bank instead of the
%                       built-in synthetic bank. Each row is one template; columns are time samples. Typical size: 600 × 100.
%                       Templates are unit-norm normalised internally, so raw or pre-normalised waveforms both work.
%                       When supplied, MaxDurMs is automatically set to match the template width (Nsamples / fs * 1000 ms).
%                       Set [] to use the built-in synthetic bank. [default: []]
%   'ThresholdSD'       Pass-1 threshold in multiples of estimated noise SD. Deliberately liberal to maximise recall.
%   'ThresholdSD2'      Pass-2 threshold. Tighter because templates are now signal-specific. 
%   'MaxDurMs'          Maximum event duration in ms (≤50 ms). Determines snippet length and longest template.
%                       [default: 50]
%   'RefractoryMs'      Minimum gap between accepted events (ms). Applied independently after each pass.
%                       [default: 50]
%   'Polarity'          'pos' | 'neg' | 'both' [default: 'both']
%   'ReduceBank'        true | false — cluster the Pass-1 template bank down to NumPrototypes representative prototypes before
%                       running template matching. Reduces peak memory from (N × Ntemplates) to (N × NumPrototypes).
%                       Automatically set to true when CustomTemplates is provided; false otherwise.
%                       [default: true]
%   'NumPrototypes'     Number of prototype templates to keep after bank reduction.  Ignored when ReduceBank is false.
%                       [default: 60]
%   'NumClusters'       Number of k-means clusters (= event types). 'auto' = estimate with elbow method.
%                       [default: 'auto']
%   'MaxClusters'       Upper bound when NumClusters='auto'. [default: 8]
%   'NumPCs'            Principal components retained for clustering and visualisation. [default: 3]
%   'MinEventsLearn'    Minimum Pass-1 events needed for learning. If fewer are found, Pass-2 reuses universal templates.
%                       [default: 20]
%   'PlotTrace'         Full-trace overview figure.     [default: true]
%   'PlotWaveforms'     Per-cluster waveform figure.    [default: true]
%   'PlotPCs'           PC scatter figure.              [default: true]
%   'MaxWaveOverlay'    Max individual waveforms overlaid per cluster panel. [default: 200]
%
%   OUTPUT results — struct with fields
%       .times_s           [K×1]  event times (seconds)
%       .times_samp        [K×1]  event times (samples, 1-indexed)
%       .amplitudes        [K×1]  signed peak amplitude (data units)
%       .clusters          [K×1]  cluster index (1..NumClusters)
%       .waveforms         [K×S]  raw waveform snippets (S = snip samples)
%       .waveforms_norm    [K×S]  unit-norm waveforms used for PCA
%       .pc_scores         [K×P]  projections onto learned PCs
%       .pc_loadings       [S×P]  PC eigenvectors
%       .pc_variance       [P×1]  fraction of variance per PC
%       .match_scores_p1   [K1×1] Pass-1 match scores (all candidates)
%       .match_scores_p2   [K×1]  Pass-2 match scores (final events)
%       .templates_univ    struct  universal template bank
%       .templates_learned struct  learned templates
%       .noise_sd          scalar  MAD noise estimate
%       .threshold_p1      scalar  Pass-1 threshold (data units)
%       .threshold_p2      scalar  Pass-2 threshold (data units)
%       .params            struct  all parameters used

    % PARSE INPUTS
    p = inputParser;
    p.addRequired ('data',                    @(x) isnumeric(x) && isvector(x));
    p.addRequired ('fs',                      @(x) isscalar(x) && x > 0);
    p.addParameter('CustomTemplates',[],      @(x) isempty(x)||(isnumeric(x)&&ismatrix(x)));
    p.addParameter('ThresholdSD',    20,      @(x) isscalar(x) && x > 0);
    p.addParameter('ThresholdSD2',   40,      @(x) isscalar(x) && x > 0);
    p.addParameter('MaxDurMs',       35,      @(x) isscalar(x) && x > 0 && x <= 50);
    p.addParameter('RefractoryMs',   100,     @(x) isscalar(x) && x > 0);
    p.addParameter('Polarity',       'pos',   @(x) any(strcmp(x,{'pos','neg','both'})));
    % p.addParameter('BandpassHz',     [],      @(x) isempty(x)||(isnumeric(x)&&numel(x)==2));
    p.addParameter('ReduceBank',     'auto',  @(x) (ischar(x)&&strcmp(x,'auto'))||islogical(x));
    p.addParameter('NumPrototypes',  20,      @(x) isscalar(x) && x >= 1);
    p.addParameter('NumClusters',    'auto',  @(x) (ischar(x)&&strcmp(x,'auto'))||(isscalar(x)&&x>=1));
    p.addParameter('MaxClusters',    5,       @(x) isscalar(x) && x >= 2);
    p.addParameter('NumPCs',         3,       @(x) isscalar(x) && x >= 1);
    p.addParameter('MinEventsLearn', 20,      @(x) isscalar(x) && x > 0);
    p.addParameter('PlotTrace',      true,    @islogical);
    p.addParameter('PlotWaveforms',  true,    @islogical);
    p.addParameter('PlotPCs',        true,    @islogical);
    p.addParameter('MaxWaveOverlay', 200,     @(x) isscalar(x) && x > 0);
    p.parse(data, fs, varargin{:});
    prm = p.Results;

    data = double(data(:));
    % N    = numel(data);
    fprintf('\n=== detect_events (two-pass, KS4-inspired) ===\n');

    %  1.  OPTIONAL BANDPASS FILTER
    data_filt = data;
    % if ~isempty(prm.BandpassHz)
    %     lo = prm.BandpassHz(1);  hi = prm.BandpassHz(2);
    %     assert(lo > 0 && lo < hi && hi < fs/2, ...
    %         'BandpassHz must satisfy 0 < lo < hi < %.1f Hz (Nyquist)', fs/2);
    %     [b, a]    = butter(4, [lo hi] / (fs/2), 'bandpass');
    %     data_filt = filtfilt(b, a, data);
    %     fprintf('[Filter]    Bandpass %.1f – %.1f Hz applied\n', lo, hi);
    % end

    %  2.  NOISE ESTIMATION  (MAD, robust to large transients — same as KS4)
    noise_sd  = mad(data_filt,1); % median(abs(data_filt)) / 0.6745;
    thresh_p1 = prm.ThresholdSD  * noise_sd;
    % thresh_p2 = prm.ThresholdSD2 * noise_sd;
    fprintf('[Noise]     SD=%.5g | Thr-P1=%.5g (%.1fσ)\n', ...
            noise_sd, thresh_p1, prm.ThresholdSD);

    % Shared geometry  — may be overridden by CustomTemplates below
    refr_samp = round(prm.RefractoryMs * fs / 1000);

    %  3. TEMPLATE BANK  (custom or built-in synthetic)
    if ~isempty(prm.CustomTemplates)
        % Custom waveform bank
        CT = double(prm.CustomTemplates);
        assert(size(CT,1) >= 1, 'CustomTemplates must have at least 1 row.');
        assert(size(CT,2) >= 4, 'CustomTemplates must have at least 4 columns (samples).');

        n_ctmpl   = size(CT,1);
        tmpl_samp = size(CT,2); % template width in samples

        % Snippet geometry is driven by the template width
        snip_samp = tmpl_samp;
        half_snip = floor(snip_samp / 2);
        % Recompute MaxDurMs so downstream figures use correct time axis
        prm.MaxDurMs = snip_samp / fs * 1000;

        % Convert matrix rows → cell array of column vectors (internal format)
        tmpl_cells = cell(n_ctmpl, 1);
        for ti = 1:n_ctmpl
            w  = CT(ti,:)';
            nr = norm(w);
            if nr > 0, w = w / nr; end
            tmpl_cells{ti} = w;
        end

        univ.waveforms = tmpl_cells;
        univ.n         = n_ctmpl;
        univ.meta      = struct('source', repmat({'custom'},n_ctmpl,1));
        univ.is_custom = true;

        fprintf('[Templates] Custom bank: %d templates  (%d samples = %.1f ms each)\n', ...
                n_ctmpl, tmpl_samp, prm.MaxDurMs);
    else
        error('CustomTemplates is empty. Need templates to match in the signal.');
    end

    %  3b.  OPTIONAL BANK REDUCTION
    %       Cluster the template bank to NumPrototypes representatives.
    %       Eliminates the O(N × Ntemplates) score_map entirely by
    %       reducing nt before template_match is ever called.
    do_reduce = prm.ReduceBank;
    if ischar(do_reduce) && strcmp(do_reduce,'auto')
        do_reduce = univ.is_custom;   % auto = true only for custom banks
    end

    if do_reduce && univ.n > prm.NumPrototypes
        fprintf('[Reduce]    Clustering %d templates → %d prototypes ...\n', ...
                univ.n, prm.NumPrototypes);
        univ = reduce_template_bank(univ, prm.NumPrototypes);
        fprintf('[Reduce]    Done. Bank size: %d\n', univ.n);
    elseif do_reduce
        fprintf('[Reduce]    Bank already has %d templates (≤ NumPrototypes=%d); skipping.\n', ...
                univ.n, prm.NumPrototypes);
    end

    %  4.  PASS 1 — UNIVERSAL DETECTION
    fprintf('[Pass 1]    Running universal template matching ...\n');

    [score_p1, best_uid] = template_match(data_filt, univ.waveforms, prm.Polarity);
    [cand_p1, sc_p1]     = find_peaks(score_p1, thresh_p1);
    tid_p1               = best_uid(cand_p1);    % best universal template ID

    fprintf('[Pass 1]    Raw candidates: %d\n', numel(cand_p1));
    [cand_p1, sc_p1, ~] = deduplicate(cand_p1, sc_p1, tid_p1, refr_samp);
    fprintf('[Pass 1]    After deduplication: %d\n', numel(cand_p1));

    wav_p1 = extract_waveforms(data_filt, cand_p1, snip_samp, half_snip);

    %  5.  TEMPLATE LEARNING
    %      PCA on unit-norm Pass-1 waveforms → k-means → learned templates
    n_p1 = numel(cand_p1);
    use_learning = n_p1 >= prm.MinEventsLearn;

    if ~use_learning
        warning('[Learning]  Only %d Pass-1 events (min=%d). Skipping; Pass-2 uses universal templates.', ...
                n_p1, prm.MinEventsLearn);
        learned = struct();
        learned.templates   = univ.waveforms;
        learned.n           = univ.n;
        learned.labels      = ones(n_p1, 1);
        learned.pc_loadings = [];
        learned.pc_variance = [];
        learned.cluster_wav = {};
        n_clust = univ.n;
        % pc_scores_p1 = zeros(n_p1, prm.NumPCs);
    else
        fprintf('[Learning]  Running PCA + k-means on %d waveforms ...\n', n_p1);
        [learned, ~, n_clust] = learn_templates( ...
            wav_p1, prm.NumPCs, prm.NumClusters, prm.MaxClusters);
        fprintf('[Learning]  %d clusters learned\n', n_clust);
        for ci = 1:n_clust
            ni = sum(learned.labels == ci);
            fprintf('            Cluster %d: %d waveforms (%.0f%%)\n', ci, ni, 100*ni/n_p1);
        end
        pve = learned.pc_variance * 100;
        fprintf('[Learning]  PC variance: ');
        fprintf('PC%d=%.1f%% ', [1:numel(pve); pve']); fprintf('\n');
    end

    %  6.  PASS 2 — LEARNED-TEMPLATE DETECTION
    fprintf('[Pass 2]    Running learned-template matching ...\n');

    [score_p2, best_cid_map] = template_match(data_filt, learned.templates, prm.Polarity);
    noise_sd  = mad(score_p2,1);
    thresh_p2 = prm.ThresholdSD2 * noise_sd;
    fprintf('[Noise]     SD=%.5g | Thr-P2=%.5g (%.1fσ)\n', ...
            noise_sd, thresh_p2, prm.ThresholdSD2);

    [cand_p2, sc_p2]         = find_peaks(score_p2, thresh_p2);
    cid_p2_raw               = best_cid_map(cand_p2);   % cluster at each peak

    fprintf('[Pass 2]    Raw candidates: %d\n', numel(cand_p2));
    [cand_p2, sc_p2, cid_p2] = deduplicate(cand_p2, sc_p2, cid_p2_raw, refr_samp);
    fprintf('[Pass 2]    After deduplication: %d\n', numel(cand_p2));

    wav_p2 = extract_waveforms(data_filt, cand_p2, snip_samp, half_snip);
    amp_p2 = zeros(numel(cand_p2), 1);
    for k = 1:numel(cand_p2)
        amp_p2(k) = data_filt(cand_p2(k));
    end

    fprintf('[Pass 2]    Final events: %d\n', numel(cand_p2));

    %  7.  PCA PROJECTION OF PASS-2 WAVEFORMS
    n_p2 = numel(cand_p2);
    if n_p2 > 0 && ~isempty(learned.pc_loadings)
        wav_norm_p2  = normalise_waveforms(wav_p2);
        n_pc_avail   = size(learned.pc_loadings, 2);
        pc_scores_p2 = wav_norm_p2 * learned.pc_loadings(:, 1:n_pc_avail);
        pc_loadings  = learned.pc_loadings;
        pc_variance  = learned.pc_variance;
    elseif n_p2 > 0
        % Learning was skipped — compute PCA directly on Pass-2 waveforms
        wav_norm_p2 = normalise_waveforms(wav_p2);
        n_pc = min(prm.NumPCs, min(n_p2-1, snip_samp));
        [~, S, V]    = svd(wav_norm_p2, 'econ');
        sv           = diag(S);
        pc_loadings  = V(:, 1:n_pc);
        pc_variance  = sv(1:n_pc).^2 / sum(sv.^2);
        pc_scores_p2 = wav_norm_p2 * pc_loadings;
    else
        wav_norm_p2  = zeros(0, snip_samp);
        pc_scores_p2 = zeros(0, prm.NumPCs);
        pc_loadings  = [];
        pc_variance  = [];
        cid_p2       = zeros(0, 1);
        amp_p2       = zeros(0, 1);
    end

    %  8.  PACK RESULTS
    results.times_s           = cand_p2 / fs;
    results.times_samp        = cand_p2;
    results.amplitudes        = amp_p2;
    results.clusters          = cid_p2;
    results.waveforms         = wav_p2;
    results.waveforms_norm    = wav_norm_p2;
    results.pc_scores         = pc_scores_p2;
    results.pc_loadings       = pc_loadings;
    results.pc_variance       = pc_variance;
    results.match_scores_p1   = sc_p1;
    results.match_scores_p2   = sc_p2;
    results.templates_univ    = univ;
    results.templates_learned = learned;
    results.noise_sd          = noise_sd;
    results.threshold_p1      = thresh_p1;
    results.threshold_p2      = thresh_p2;
    results.params            = prm;
    results.fs                = fs;        % stored for downstream use

    %  9.  VISUALISATION
    cmap = cluster_cmap(n_clust);

    if prm.PlotTrace
        fig_trace(data_filt, fs, results, cmap, n_clust, prm);
    end
    if prm.PlotWaveforms && n_p2 > 0
        fig_waveforms(results, fs, cmap, n_clust, snip_samp, prm);
    end
    if prm.PlotPCs && n_p2 > 0 && ~isempty(pc_scores_p2) && size(pc_scores_p2,2) >= 2
        fig_pcs(results, fs, cmap, n_clust, snip_samp, prm);
    end

    fprintf('=== detect_events complete: %d events in %d cluster(s) ===\n\n', ...
            n_p2, n_clust);
end

function [combined, best_id] = template_match(sig, templates, polarity)
    %  Cross-correlates the signal with each template one at a time, updating
    %  running maximum scores.  Peak memory is 4 × N bytes (one single vector)
    %  regardless of the number of templates — compared to 4 × N × nt bytes
    %  in the previous matrix version.
    N  = numel(sig);
    nt = numel(templates);

    combined = zeros(N, 1, 'single');   % running max score   — O(N)
    best_id  = ones(N,  1, 'uint16');   % best template index — O(N)
    row      = zeros(N, 1, 'single');   % scratch buffer      — O(N)

    for ti = 1:nt
        tmpl = single(templates{ti});
        tmpl = tmpl / (norm(tmpl) + eps('single'));
        row(:) = single(conv(double(sig), flipud(double(tmpl)), 'same'));

        switch polarity
            case 'pos',  s = row;
            case 'neg',  s = -row;
            case 'both', s = abs(row);
        end

        better = s > combined;
        combined(better) = s(better);
        best_id(better)  = uint16(ti);
    end

    combined = double(combined);
    best_id  = double(best_id);
end

function univ_out = reduce_template_bank(univ, n_proto)
    %  reduce_template_bank
    %  Cluster Ntemplates waveforms to NumPrototypes via k-means in PCA space.
    %  Each prototype is the centroid of its cluster projected back to waveform
    %  space and unit-normalised.
    nt = univ.n;
    S  = numel(univ.waveforms{1});

    % Stack all templates into a matrix [nt × S]
    W = zeros(nt, S);
    for ti = 1:nt
        w = univ.waveforms{ti}(:)';
        W(ti,:) = w / (norm(w) + eps);
    end

    % PCA: reduce to min(n_proto-1, S) dimensions for stable k-means
    n_pc = min(n_proto - 1, min(nt - 1, S));
    [~, ~, V] = svd(W, 'econ');
    % sv  = diag(Sv);
    L   = V(:, 1:n_pc);   % [S × n_pc]  loadings
    X   = W * L;          % [nt × n_pc] PC scores

    % k-means on PC scores
    rng(0);
    try
        opts = statset('MaxIter',300,'Display','off');
        [labels, centres] = kmeans(X, n_proto, ...
            'Replicates',3, 'Options',opts, 'Distance','sqeuclidean');
    catch
        [labels, centres] = simple_kmeans(X, n_proto);
    end

    % Project centres back to waveform space and normalise
    proto_wav = centres * L';     % [n_proto × S]
    proto_cells = cell(n_proto,1);
    for ci = 1:n_proto
        w  = proto_wav(ci,:)';
        nr = norm(w);
        if nr > 0, w = w / nr; end
        proto_cells{ci} = w;
    end

    % Report worst-case approximation error (mean cosine distance)
    cos_err = zeros(nt,1);
    for ti = 1:nt
        ci = labels(ti);
        cos_err(ti) = 1 - dot(W(ti,:), proto_wav(ci,:)) / ...
                          (norm(W(ti,:)) * norm(proto_wav(ci,:)) + eps);
    end
    fprintf('[Reduce]    Mean cosine error: %.4f  |  Max: %.4f\n', ...
            mean(cos_err), max(cos_err));

    univ_out           = univ;
    univ_out.waveforms = proto_cells;
    univ_out.n         = n_proto;
    univ_out.meta      = struct('source', repmat({'prototype'},n_proto,1));
    univ_out.labels_orig = labels;   % which original template → which proto
end

function [peaks, scores] = find_peaks(score, thresh)
    %  One peak per contiguous above-threshold region (local maximum).
    above = score >= thresh;
    peaks  = [];
    scores = [];
    N = numel(score);
    s = 1;
    while s <= N
        if above(s)
            e = s;
            while e < N && above(e+1), e = e+1; end
            [sc, rel] = max(score(s:e));
            peaks  = [peaks;  s + rel - 1];
            scores = [scores; sc];
            s = e + 1;
        else
            s = s + 1;
        end
    end
end

function [idx_out, sc_out, id_out] = deduplicate(idx_in, sc_in, id_in, refr)
    %  Greedy refractory suppression — keep highest score within window.
    if isempty(idx_in)
        idx_out = []; sc_out = []; id_out = []; return
    end
    [idx_in, ord] = sort(idx_in);
    sc_in  = sc_in(ord);
    id_in  = id_in(ord);
    n    = numel(idx_in);
    kept = true(n,1);
    for i = 1:n
        if ~kept(i), continue; end
        for j = i+1:n
            if ~kept(j), continue; end
            if idx_in(j) - idx_in(i) > refr, break; end
            if sc_in(j) > sc_in(i)
                kept(i) = false; break;
            else
                kept(j) = false;
            end
        end
    end
    idx_out = idx_in(kept);
    sc_out  = sc_in(kept);
    id_out  = id_in(kept);
end

function W = extract_waveforms(sig, peaks, snip_samp, half_snip)
    %  Cut snip_samp-long snippets centred on each peak (zero-padded at edges).
    N = numel(sig);
    K = numel(peaks);
    W = zeros(K, snip_samp);
    for k = 1:K
        i1 = peaks(k) - half_snip + 1;
        i2 = i1 + snip_samp - 1;
        ss = max(i1,1);  se = min(i2,N);
        ds = ss - i1 + 1;
        de = ds + (se - ss);
        W(k, ds:de) = sig(ss:se)';
    end
end

function Wn = normalise_waveforms(W)
    %  Unit-L2-norm each waveform row (shape preserved, amplitude removed).
    nr = sqrt(sum(W.^2, 2));
    Wn = W;
    nz = nr > 0;
    Wn(nz,:) = W(nz,:) ./ nr(nz);
end

function [learned, pc_scores, n_clust] = learn_templates(W, n_pc_req, nc_param, max_clust)
    %  PCA → k-means → learned templates, one per cluster.
    Wn   = normalise_waveforms(W);
    K    = size(Wn,1);
    S    = size(Wn,2);
    n_pc = min(n_pc_req, min(K-1, S));

    % PCA via compact SVD
    [~, Sv, V] = svd(Wn, 'econ');
    sv         = diag(Sv);
    total_var  = sum(sv.^2);
    loadings   = V(:, 1:n_pc);          % [S × n_pc]
    var_exp    = sv(1:n_pc).^2 / total_var;
    pc_scores  = Wn * loadings;         % [K × n_pc]

    % Number of clusters
    if ischar(nc_param) && strcmp(nc_param, 'auto')
        n_clust = elbow_clusters(pc_scores, min(max_clust, floor(K/5)));
    else
        n_clust = min(round(nc_param), floor(K/5));
        n_clust = max(n_clust, 1);
    end

    % k-means in PC space
    rng(42);
    try
        opts = statset('MaxIter',300,'Display','off');
        [labels, centres] = kmeans(pc_scores, n_clust, ...
            'Replicates',5, 'Options',opts, 'Distance','sqeuclidean');
    catch
        [labels, centres] = simple_kmeans(pc_scores, n_clust);
    end

    % Project centres back to waveform space → learned templates
    tmpl_wav = centres * loadings';   % [n_clust × S]
    tmpl_cell = cell(n_clust,1);
    for ci = 1:n_clust
        t = tmpl_wav(ci,:)';
        nr = norm(t);
        if nr > 0, t = t / nr; end
        tmpl_cell{ci} = t;
    end

    % Collect raw waveforms per cluster (for waveform figure)
    cluster_wav = cell(n_clust,1);
    for ci = 1:n_clust
        mask = labels == ci;
        cluster_wav{ci} = W(mask,:);
    end

    learned.templates   = tmpl_cell;
    learned.n           = n_clust;
    learned.labels      = labels;
    learned.centres_pc  = centres;
    learned.pc_loadings = loadings;
    learned.pc_variance = var_exp;
    learned.cluster_wav = cluster_wav;
end

function k_opt = elbow_clusters(X, k_max)
    %  elbow_clusters  — greatest second-difference of within-cluster SS
    k_max = max(2, min(k_max, size(X,1)-1));
    wcss  = zeros(k_max,1);
    rng(0);
    for k = 1:k_max
        try
            opts = statset('MaxIter',200,'Display','off');
            [~, c] = kmeans(X, k, 'Replicates',3,'Options',opts);
        catch
            [~, c] = simple_kmeans(X, k);
        end
        mn = min(pdist2(X, c,'squaredeuclidean'), [], 2);
        wcss(k) = sum(mn);
    end
    if k_max >= 3
        d2 = diff(diff(wcss));
        [~, i] = max(d2);
        k_opt = max(1, i + 1);
    else
        k_opt = 2;
    end
end

function [labels, centres] = simple_kmeans(X, k)
    %  simple_kmeans  — fallback without Statistics Toolbox
    [n, ~] = size(X);
    idx  = randperm(n, k);
    centres = X(idx,:);
    labels  = ones(n,1);
    for iter = 1:300
        D = zeros(n, k);
        for ci = 1:k
            diff_ = X - centres(ci,:);
            D(:,ci) = sum(diff_.^2, 2);
        end
        [~, new_labels] = min(D, [], 2);
        if isequal(new_labels, labels) && iter > 1, break; end
        labels = new_labels;
        for ci = 1:k
            m = X(labels==ci,:);
            if ~isempty(m), centres(ci,:) = mean(m,1); end
        end
    end
end

function cmap = cluster_cmap(n)
    %  cluster_cmap  — distinct colours for up to 12 clusters
    base = [0.22 0.49 0.82;  0.89 0.35 0.13;  0.18 0.70 0.48;
            0.75 0.20 0.60;  0.93 0.70 0.12;  0.14 0.62 0.70;
            0.85 0.32 0.42;  0.43 0.73 0.25;  0.55 0.38 0.70;
            0.95 0.50 0.22;  0.30 0.30 0.30;  0.70 0.70 0.20];
    cmap = base(mod((1:n)-1, size(base,1))+1, :);
end

function fig_trace(sig, fs, res, cmap, n_clust, prm)
    %  fig_trace — overview: full signal + Pass-1 candidates + Pass-2 events
    N  = numel(sig);
    t  = (0:N-1) / fs;
    K  = numel(res.times_s);

    figure('Name','detect_events — Signal trace','Color','w', ...
           'Position',[60 60 1200 520]);

    % Top panel: full trace
    ax1 = subplot(3,1,[1 2]);
    plot(ax1, t, sig, 'Color',[0.42 0.42 0.42],'LineWidth',0.5);
    hold(ax1,'on');

    yline(ax1,  res.threshold_p1,'--','Color',[0.80 0.25 0.25],'LineWidth',1.2,...
          'Label',sprintf('P1 (%.1fσ)', prm.ThresholdSD),'LabelHorizontalAlignment','left');
    yline(ax1, -res.threshold_p1,'--','Color',[0.80 0.25 0.25],'LineWidth',1.2);
    yline(ax1,  res.threshold_p2,'--','Color',[0.20 0.50 0.80],'LineWidth',1.2,...
          'Label',sprintf('P2 (%.1fσ)', prm.ThresholdSD2),'LabelHorizontalAlignment','left');
    yline(ax1, -res.threshold_p2,'--','Color',[0.20 0.50 0.80],'LineWidth',1.2);

    leg_h = [];
    for ci = 1:n_clust
        mask = res.clusters == ci;
        if ~any(mask), continue; end
        h = plot(ax1, res.times_s(mask), res.amplitudes(mask), ...
            'v','MarkerSize',9,'MarkerFaceColor',cmap(ci,:), ...
            'MarkerEdgeColor','w','LineWidth',0.5, ...
            'DisplayName',sprintf('Cluster %d  (n=%d)', ci, sum(mask)));
        leg_h = [leg_h h];
    end
    if ~isempty(leg_h), legend(ax1, leg_h,'Location','northeast','FontSize',9); end
    xlabel(ax1,'Time (s)'); ylabel(ax1,'Amplitude');
    title(ax1, sprintf('Two-pass detection — %d events, %d cluster(s)', K, n_clust));
    axis(ax1,'tight'); grid(ax1,'on');

    % Bottom panel: match-score trace (Pass-2)
    ax2 = subplot(3,1,3);
    if ~isempty(res.match_scores_p2)
        % Reconstruct a sparse score trace for display
        score_disp = zeros(K,1);
        for k2 = 1:K
            score_disp(k2) = res.match_scores_p2(k2);
        end
        stem(ax2, res.times_s, res.match_scores_p2, ...
             'Marker','none','LineWidth',0.8,'Color',[0.5 0.5 0.5]);
        hold(ax2,'on');
        for ci = 1:n_clust
            mask = res.clusters == ci;
            if ~any(mask), continue; end
            stem(ax2, res.times_s(mask), res.match_scores_p2(mask), ...
                 'Marker','o','MarkerSize',5,'MarkerFaceColor',cmap(ci,:), ...
                 'Color',cmap(ci,:),'LineWidth',0.8);
        end
        yline(ax2, res.threshold_p2,'--','Color',[0.20 0.50 0.80],'LineWidth',1.2);
        xlabel(ax2,'Time (s)'); ylabel(ax2,'Match score (P2)');
        title(ax2,'Pass-2 match scores at detected events');
        axis(ax2,'tight'); grid(ax2,'on');
    end

    sgtitle(sprintf('detect\\_events  |  P1 thr=%.1fσ  P2 thr=%.1fσ  refractory=%.0f ms  MaxDur=%.0f ms', ...
            prm.ThresholdSD, prm.ThresholdSD2, prm.RefractoryMs, prm.MaxDurMs), ...
            'FontSize',10);
end

function fig_waveforms(res, fs, cmap, n_clust, snip_samp, prm)
    %  fig_waveforms — per-cluster waveform panels + learned template overlay
    t_ms = ((0:snip_samp-1) / fs * 1000) - prm.MaxDurMs/2;

    cols  = min(n_clust, 4);
    n_rows_clust = ceil(n_clust / cols);
    n_rows_total = n_rows_clust + 1;   % +1 row for learned templates

    fig = figure('Name','detect_events — Waveforms','Color','w', ...
                 'Position',[80 80 290*cols 230*n_rows_total]);

    % Per-cluster panels
    for ci = 1:n_clust
        ax = subplot(n_rows_total, cols, ci, 'Parent', fig);
        mask = res.clusters == ci;
        W    = res.waveforms(mask, :);
        n_ci = size(W,1);

        if n_ci == 0
            text(ax,0.5,0.5,'(empty)','Units','normalized', ...
                 'HorizontalAlignment','center','Color',cmap(ci,:));
            axis(ax,'off');
            title(ax, sprintf('Cluster %d — empty',ci),'Color',cmap(ci,:));
            continue
        end

        % Thin individual traces (subsampled if large)
        n_ov = min(n_ci, prm.MaxWaveOverlay);
        idx  = unique(round(linspace(1, n_ci, n_ov)));
        col_lt = min(cmap(ci,:) * 0.45 + 0.55, 1);  % lighter tint
        plot(ax, t_ms, W(idx,:)', 'Color',[col_lt 0.20], 'LineWidth',0.5);
        hold(ax,'on');

        % ±1 SD band + mean
        mu = mean(W,1);
        sg = std(W,0,1);
        fill(ax, [t_ms, fliplr(t_ms)], [mu+sg, fliplr(mu-sg)], ...
             cmap(ci,:), 'FaceAlpha',0.18, 'EdgeColor','none');
        plot(ax, t_ms, mu, 'Color',cmap(ci,:), 'LineWidth',2.2);
        xline(ax, 0, ':k', 'LineWidth',0.8);

        xlabel(ax,'Time (ms)'); ylabel(ax,'Amplitude');
        title(ax, sprintf('Cluster %d  (n=%d)', ci, n_ci), ...
              'Color',cmap(ci,:), 'FontWeight','bold');
        axis(ax,'tight'); grid(ax,'on');
    end

    % Bottom row: learned templates overlaid
    ax_t = subplot(n_rows_total, cols, ...
                   (n_rows_total-1)*cols+1 : n_rows_total*cols, ...
                   'Parent', fig);
    hold(ax_t,'on');
    for ci = 1:n_clust
        if ci > numel(res.templates_learned.templates), continue; end
        tmpl   = res.templates_learned.templates{ci};
        n_tmpl = numel(tmpl);
        t_t    = linspace(t_ms(1), t_ms(end), n_tmpl);
        tmpl_n = tmpl / (max(abs(tmpl)) + eps);
        plot(ax_t, t_t, tmpl_n, 'Color',cmap(ci,:), 'LineWidth',2.0, ...
             'DisplayName', sprintf('Cluster %d', ci));
    end
    xline(ax_t, 0, ':k', 'LineWidth',0.8);
    xlabel(ax_t,'Time relative to peak (ms)'); ylabel(ax_t,'Norm. amplitude');
    title(ax_t,'Learned templates (unit-normalised, overlaid)');
    legend(ax_t,'show','Location','northeast','FontSize',9);
    grid(ax_t,'on'); axis(ax_t,'tight');

    sgtitle(fig, sprintf('Waveform gallery  |  %d clusters  |  snippet=%.0f ms', ...
            n_clust, prm.MaxDurMs), 'FontSize',11);
end

function fig_pcs(res, fs, cmap, n_clust, snip_samp, prm)
    % fig_pcs — PC scatter plots + amplitude vs PC1 + PC loadings shape
    scores  = res.pc_scores;
    var_exp = res.pc_variance * 100;
    n_pc    = size(scores,2);
    t_ms    = ((0:snip_samp-1) / fs * 1000) - prm.MaxDurMs/2;

    % Layout: [PC1v2] [PC1v3 or IEI]  [Amp v PC1]
    %         [PC loadings ————————————————————————]
    n_top = 3;   % always 3 columns on top row
    fig = figure('Name','detect_events — PC space','Color','w', ...
                 'Position',[100 100 1150 680]);

    % Helper: draw one scatter panel
    function draw_scatter(ax, xi, yi)
        hold(ax,'on');
        for ci = 1:n_clust
            mask = res.clusters == ci;
            if ~any(mask), continue; end
            scatter(ax, scores(mask,xi), scores(mask,yi), 30, ...
                    cmap(ci,:), 'filled', 'MarkerFaceAlpha',0.65, ...
                    'DisplayName', sprintf('Cluster %d  (n=%d)',ci,sum(mask)));
        end
        % Centroids
        for ci = 1:n_clust
            mask = res.clusters == ci;
            if sum(mask) < 2, continue; end
            plot(ax, mean(scores(mask,xi)), mean(scores(mask,yi)), ...
                 '+','Color',cmap(ci,:)*0.5,'MarkerSize',14,'LineWidth',2.5, ...
                 'HandleVisibility','off');
        end
        xlabel(ax, sprintf('PC%d  (%.1f%%)', xi, var_exp(xi)));
        ylabel(ax, sprintf('PC%d  (%.1f%%)', yi, var_exp(yi)));
        title(ax, sprintf('PC%d vs PC%d', xi, yi));
        legend(ax,'show','Location','best','FontSize',8);
        grid(ax,'on'); axis(ax,'tight');
    end

    % Panel 1 — PC1 vs PC2
    ax1 = subplot(2, n_top, 1, 'Parent', fig);
    draw_scatter(ax1, 1, 2);

    % Panel 2 — PC1 vs PC3 (or PC2 vs PC3 if only 2 PCs)
    ax2 = subplot(2, n_top, 2, 'Parent', fig);
    if n_pc >= 3
        draw_scatter(ax2, 1, 3);
    else
        draw_scatter(ax2, 2, 1);  % duplicate with axes flipped for inspection
        title(ax2,'PC2 vs PC1 (only 2 PCs available)');
    end

    % Panel 3 — Amplitude vs PC1
    ax3 = subplot(2, n_top, 3, 'Parent', fig);
    hold(ax3,'on');
    for ci = 1:n_clust
        mask = res.clusters == ci;
        if ~any(mask), continue; end
        scatter(ax3, scores(mask,1), abs(res.amplitudes(mask)), 25, ...
                cmap(ci,:), 'filled', 'MarkerFaceAlpha',0.55, ...
                'DisplayName', sprintf('Cluster %d',ci));
    end
    xlabel(ax3, sprintf('PC1  (%.1f%%)', var_exp(1)));
    ylabel(ax3,'|Amplitude| (data units)');
    title(ax3,'Amplitude vs PC1');
    legend(ax3,'show','Location','best','FontSize',8);
    grid(ax3,'on'); axis(ax3,'tight');

    % Panel 4 — PC loadings (bottom row, spanning all columns)
    ax4 = subplot(2, n_top, n_top+1 : 2*n_top, 'Parent', fig);
    hold(ax4,'on');
    n_pc_show = min(n_pc, 5);
    lc = lines(n_pc_show);
    offset_step = 2.5;
    for pi = 1:n_pc_show
        ld = res.pc_loadings(:, pi);
        ld_n  = ld / (max(abs(ld)) + eps);
        off   = (n_pc_show - pi) * offset_step;
        plot(ax4, t_ms, ld_n + off, 'Color',lc(pi,:), 'LineWidth',1.8, ...
             'DisplayName', sprintf('PC%d  (%.1f%%)', pi, var_exp(pi)));
        % Zero line for this PC
        plot(ax4, [t_ms(1) t_ms(end)], [off off], ':', ...
             'Color',lc(pi,:)*0.5,'LineWidth',0.6,'HandleVisibility','off');
    end
    xline(ax4, 0, ':k','LineWidth',0.8,'HandleVisibility','off');
    xlabel(ax4,'Time relative to detection peak (ms)');
    ylabel(ax4,'Loading (offset for clarity)');
    title(ax4,'PC loadings — temporal waveform decomposition');
    legend(ax4,'show','Location','northeast','FontSize',9);
    grid(ax4,'on'); axis(ax4,'tight');

    sgtitle(fig, sprintf('PC space  |  top %d PCs (%.1f%% total var.)  |  %d clusters  |  %d events', ...
            n_pc, sum(var_exp), n_clust, numel(res.times_s)), 'FontSize',11);
end