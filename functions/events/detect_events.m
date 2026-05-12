function collision_events = detect_events(templates, jerk, par, see_examples)
    %% PARAMETERS
    templates    = single(templates);
    jerk         = single(jerk(:));
    Kcentroids   = round(size(templates,1)/10);
    Kpc          = 6;
    blocksize    = 2^18;
    tol_samples  = 200;
    votes        = 2;

    % MAD-based score threshold. Raise to be stricter on shape match
    % par.kMAD = 8;
    % par.minAmplFraction = 0.5; % Lower amplitude gate
    % par.maxAmplFraction = 1.5; % Upper amplitude gate
    templatePeaks   = max(abs(templates), [], 2);
    meanTemplatePeak = mean(templatePeaks);
    minAmplitude    = par.minAmplFraction * meanTemplatePeak;   
    maxAmplitude    = par.maxAmplFraction * meanTemplatePeak;

    fprintf('Amplitude gate:  [%.2e  to  %.2e]  (mean template peak: %.2e)\n', ...
            minAmplitude, maxAmplitude, meanTemplatePeak);

    %% Reduce templates
    [centroids, PCs] = reduce_templates(templates, Kpc, Kcentroids);

    %% Centroid-based detection
    centroidNorms  = sqrt(sum(centroids.^2, 2))';
    dets_centroids = batch_match_and_detect(jerk, centroids, centroidNorms, ...
                                            blocksize, par.kMAD, minAmplitude, maxAmplitude);
    fprintf('Centroid detections: %d\n', numel(dets_centroids));

    %% PCA-based detection
    pcKernels = single(PCs');
    pcKernels = pcKernels - mean(pcKernels, 2);
    pcNorms   = sqrt(sum(pcKernels.^2, 2))';
    dets_pca  = batch_match_and_detect(jerk, pcKernels, pcNorms, ...
                                       blocksize, par.kMAD, minAmplitude, maxAmplitude);
    fprintf('PCA detections: %d\n', numel(dets_pca));

    %% Consensus
    r = compare_detections(dets_centroids, dets_pca, tol_samples);
    fprintf('Centroid vs PCA matches (tol=%d): %d\n', tol_samples, r.matches);

    collision_events = build_consensus(dets_centroids, dets_pca, tol_samples, votes);

    if see_examples
        visualize_examples(jerk, collision_events, templates, 6, 500);

        % Amplitude distribution diagnostic
        if ~isempty(dets_centroids)
            figure;
            histogram([dets_centroids.amplitude], 30); hold on;
            xline(minAmplitude, 'b--', 'LineWidth', 2, 'Label', 'min threshold');
            xline(maxAmplitude, 'r--', 'LineWidth', 2, 'Label', 'max threshold');
            xlabel('Peak amplitude'); title('Centroid detection amplitudes');
        end
    end
end

%% -----------------------------------------------------------------------
function [centroids, PCs] = reduce_templates(templates, Kpc, Kcentroids)
    % Zero-center each template row
    templates_c = templates - mean(templates, 2);   % N x L

    % PCA in sample space (SVD of L x N matrix; L is small ~100)
    [U,~,~] = svd(templates_c', 'econ');
    PCs = U(:, 1:Kpc);                              % L x Kpc

    % Project and cluster
    scores = (PCs' * templates_c')';                % N x Kpc
    opts = statset('UseParallel', false, 'MaxIter', 200);
    idx  = kmeans(scores, Kcentroids, 'Replicates', 3, 'Options', opts);

    % Cluster centroids in original sample space
    L = size(templates, 2);
    centroids = zeros(Kcentroids, L, 'single');
    for k = 1:Kcentroids
        members = templates(idx == k, :);
        if ~isempty(members)
            centroids(k,:) = mean(single(members), 1);
        end
    end
    centroids = centroids - mean(centroids, 2);     % center each centroid row
end

%% -----------------------------------------------------------------------
function detections = batch_match_and_detect(signal, kernels, kernelNorms, blocksize, kMAD, minAmplitude, maxAmplitude)
    signal = signal(:);
    N = numel(signal);
    [K, L] = size(kernels);
    halfL  = floor(L / 2);

    fftlen = 2^nextpow2(blocksize + L - 1);
    Kfft   = fft(flip(kernels, 2)', fftlen);

    % Causal sliding-window energy
    cumE = [0; cumsum(double(signal.^2))];
    windowEnergy = zeros(N, 1, 'double');
    windowEnergy(L:N)   = cumE(L+1:N+1) - cumE(1:N-L+1);
    windowEnergy(1:L-1) = cumE(2:L);

    % Causal sliding-window peak amplitude
    absSignal  = double(abs(signal));
    windowPeak = movmax(absSignal, [L-1, 0]);   % cleaner than the for-loop; causal

    detections = struct('time', {}, 'kernel', {}, 'score', {});
    detCount   = 0;

    nblocks = ceil(N / blocksize);
    for b = 1:nblocks
        startIdx = (b-1)*blocksize + 1;
        endIdx   = min(b*blocksize, N);
        numValid = endIdx - startIdx + 1;

        if startIdx > L
            seg = signal(startIdx - L + 1 : endIdx);
        else
            seg = [zeros(L - startIdx, 1, 'single'); signal(1:endIdx)];
        end

        segFft  = fft(single(seg), fftlen);
        convOut = real(ifft(bsxfun(@times, segFft, Kfft)));

        valid    = convOut(L : L + numValid - 1, :);
        timeIdxs = (startIdx : endIdx)';

        denom  = sqrt(max(1e-12, windowEnergy(timeIdxs)));
        scores = bsxfun(@rdivide, valid,  denom);
        scores = bsxfun(@rdivide, scores, kernelNorms);

        for k = 1:K
            v = scores(:, k);

            % MAD-based threshold (replaces prctile)
            % Robust to rare-event contamination; consistent across blocks
            medV   = median(v);
            madV   = median(abs(v - medV)) / 0.6745;  % estimate of sigma
            thresh = medV + kMAD * madV;

            [~, locs] = findpeaks(v, 'MinPeakHeight', thresh, 'MinPeakDistance', L);

            for i = 1:numel(locs)
                tRaw  = timeIdxs(locs(i));
                tCent = max(1, tRaw - halfL);
                amp   = windowPeak(tRaw);

                % Lower amplitude gate
                if amp < minAmplitude, continue; end

                % NEW: Upper amplitude gate
                if amp > maxAmplitude, continue; end

                detCount = detCount + 1;
                detections(detCount).time      = tCent;
                detections(detCount).kernel    = k;
                detections(detCount).score     = v(locs(i));
                detections(detCount).amplitude = amp;
            end
        end
    end
end

%% -----------------------------------------------------------------------
function report = compare_detections(dA, dB, tol)
    if isempty(dA) || isempty(dB)
        report.matches = 0;  report.totalA = numel(dA);
        report.totalB  = numel(dB);
        report.uniqueA = numel(dA); report.uniqueB = numel(dB);
        return
    end
    timesA   = [dA.time]';
    timesB   = [dB.time]';
    matchedB = false(numel(timesB), 1);
    matchCount = 0;
    for i = 1:numel(timesA)
        j = find(abs(timesB - timesA(i)) <= tol, 1);
        if ~isempty(j)
            matchCount  = matchCount + 1;
            matchedB(j) = true;
        end
    end
    report.matches = matchCount;
    report.totalA  = numel(timesA);
    report.totalB  = numel(timesB);
    report.uniqueA = numel(timesA) - matchCount;
    report.uniqueB = numel(timesB) - sum(matchedB);
end

%% -----------------------------------------------------------------------
function events = build_consensus(d1, d2, tol, votes)
    t1   = [d1.time]';
    t2   = [d2.time]';
    allT = sort([t1; t2]);
    events = [];
    k = 0;
    for i = 1:numel(allT)
        t = allT(i);
        if k == 0 || abs(t - events(k).time) > tol
            k = k + 1;
            events(k).time     = t;
            events(k).centroid = any(abs(t1 - t) <= tol);
            events(k).pca      = any(abs(t2 - t) <= tol);
        else
            events(k).centroid = events(k).centroid | any(abs(t1 - t) <= tol);
            events(k).pca      = events(k).pca      | any(abs(t2 - t) <= tol);
        end
    end
    for i = 1:numel(events)
        events(i).votes = events(i).centroid + events(i).pca;
    end
    events([events.votes] < votes) = [];
end

%% -----------------------------------------------------------------------
function visualize_examples(signal, events, templates, nExamples, window)
    if isempty(events), disp('No consensus events to display.'); return; end
    L    = size(templates, 2);
    N    = numel(signal);
    nShow = min(nExamples, numel(events));
    idx  = randperm(numel(events), nShow);
    templateMean = mean(templates, 1);      % 1 x L

    figure;
    tiledlayout('flow');                    % FIX 7: argument required
    for i = 1:nShow
        t = events(idx(i)).time;
        a = max(1, t - window);
        b = min(N, t + window - 1);

        nexttile;
        plot(a:b, signal(a:b), 'k'); hold on;

        % Overlay mean template centred on detected event
        xTpl = (t - floor(L/2)) : (t + ceil(L/2) - 1);
        if xTpl(1) >= 1 && xTpl(end) <= N
            plot(xTpl, templateMean, 'r', 'LineWidth', 2);
        end
        title(sprintf('t=%d ms  votes=%d', t, events(idx(i)).votes));
        xlabel('Sample (ms)'); grid on;
    end
end