function neuralDynamics = calculate_neural_dynamics(neurons, fireRate, trialdef, opt)
% calculate_neural_dynamics  Low-dimensional embedding of trial-level
%                            population activity via PCA, t-SNE, or UMAP.
%
% PURPOSE:
%   Builds a per-trial population feature vector by binning each cluster's
%   spike train across the trial window, then projects every trial to a
%   3-D space using PCA, t-SNE, or UMAP. The intent is to expose any
%   coarse structure across trials (drift, condition clusters, stimulus
%   geometry) at the level of the whole recorded population.
%
%   NOTE on naming: the file is called "neural_dynamics" but in its
%   current form it produces "neural states per trial" (each trial is one
%   point in the embedding) rather than time-resolved trajectories (a
%   line per trial / condition). See AUDIT NOTES below.
%
% USAGE:
%   neuralDynamics = calculate_neural_dynamics(neurons, fireRate, trialdef, opt)
%
% INPUTS:
%   neurons  - struct from sort2trials (standard mode). Used field:
%                .itiOn  cell {Nclust, 1} of {Ntrials, 1} spike vectors.
%   fireRate - struct from calculate_fireRate_general. Used field:
%                .sps   cell {Nclust, 1} of [Ntrials x Nbins] matrices.
%              Currently used only to derive min/max firing rate for the
%              synthetic surrogate; the rate matrix itself is NOT used as
%              the embedding input (the function rebins from neurons).
%   trialdef - cell from trialdefGen. trialdef{2,1}(:,1:2) (start/end ms)
%              is used to compute median trial duration.
%   opt      - resolved options struct. Used subfields:
%                .neurDyn.binSize  bin width in SECONDS (default 0.1)
%                .neurDyn.method   'PCA' | 'tSNE' | 'UMAP' (default tSNE)
%                .neurDyn.synth    if true, bypass real data and embed a
%                                  Poisson surrogate instead (default true)
%
% OUTPUT:
%   neuralDynamics - currently unused / unassigned. The function only
%                    produces figures and does not return data. See AUDIT
%                    NOTES.
%
% METHODS:
%   'PCA'  - mean-centred eigendecomposition; first 3 components.
%   'tSNE' - 3-D t-SNE with perplexity=10, Standardize=true, default rng.
%   'UMAP' - requires UMAP for MATLAB (FileExchange #71902); calls
%            run_umap with n_components=3.
%
% AUDIT NOTES (open issues; see docs/audit_calculate_neural_dynamics.md
% for the full breakdown and recommended fixes):
%   1. SYNTH IS THE DEFAULT PATH. opt.neurDyn.synth defaults to true,
%      so the function embeds a homogeneous-rate Poisson surrogate rather
%      than the real spikes. The real-data branch is also broken (see #2).
%   2. REAL-DATA BRANCH IS DEAD. When synth=false, spTimes = neurons,
%      then the binning loop does spTimes{neuron, trial}. But neurons is
%      a struct with .itiOn, not a 2-D cell array. This errors.
%   3. UNITS MIX. The synth surrogate uses seconds; neurons stores
%      timestamps in ms. The histcounts edges (timeBins, in seconds)
%      do not match the real-data unit.
%   4. NOT MULTI-AREA AWARE. Hardcodes neurons.itiOn. In multi-area runs
%      neurons has shape neurons.<Area>.itiOn (since the recent rework),
%      so this function will fail in multi-area mode.
%   5. NOT ALIGNMENT-GENERIC. Hardcodes the itiOn alignment; ignores
%      opt.alignto.
%   6. TRAJECTORY PLOT DOESN'T DRAW TRAJECTORIES. The second figure
%      plots one point per trial with plot3, not connected lines through
%      time — there's no time dimension in reducedData by construction.
%   7. FIREATE IS UNDERUSED. fireRate.sps is already binned/rate-
%      normalised; this function rebins from raw spike times and only
%      uses fireRate for min/max scaling of the surrogate.
%   8. FIGURES SAVED TO PWD. savefig uses a relative filename; the
%      output should go to opt.analysis/plots/neural_dynamics/.
%   9. NO RETURN VALUE. The output variable neuralDynamics is never
%      assigned. NGL02 saves whatever this returns to neuralDynamics.mat
%      — which today is "Undefined function or variable" depending on
%      MATLAB version.
%
% SEE ALSO:
%   calculate_fireRate_general, calculate_fireRate_byBlock
%
% Last modified 29.05.2026 (Jesus) - audit + docstring pass (#15)
%
% ---- Original behaviour preserved below; refactor pending. ----

%% Default
if ~isfield(opt,'neurDyn'), opt.neurDyn = struct('binSize',0.1,'method',"tSNE",'synth',true); end
if ~isfield(opt.neurDyn,'binSize'), opt.neurDyn.binSize = 0.1;      end
if ~isfield(opt.neurDyn,'method'),  opt.neurDyn.method  = "tSNE";   end
if ~isfield(opt.neurDyn,'synth'),   opt.neurDyn.synth   = true;     end

%% Retrieve actual data and generate synth is required
numNeurons    = size(neurons.itiOn,1); % Number of neurons
numTrials     = size(neurons.itiOn{1,1},1); % Number of trials
trialDuration = median((trialdef{2,1}(:,2)-trialdef{2,1}(:,1))/1000); % Trial duration, (s)
firingRate.actual = fireRate.sps; % Mean firing rate (sp/s)
timeBins      = (0:opt.neurDyn.binSize:trialDuration+opt.neurDyn.binSize);
numTimeBins   = length(timeBins);

if opt.neurDyn.synth
    maxFR = max(max(cell2mat(fireRate.sps)));
    minFR = min(min(cell2mat(fireRate.sps)));

    firingRate.synth = minFR + (maxFR-minFR).*rand(numNeurons,1,"single"); % Mean firing rate (Hz)

    spTimes.synth = cell(numNeurons, numTrials);
    for neuron = 1:numNeurons
        for trial = 1:numTrials
            % firingRate is a struct (.actual / .synth); index into .synth.
            nSpikes = poissrnd(firingRate.synth(neuron) * trialDuration);
            spTimes.synth {neuron, trial} = sort(rand(nSpikes, 1) * trialDuration);
        end
    end
else
    spTimes = neurons;
end

%% Bin spike data
binnedSpikes = zeros(numNeurons, numTimeBins, numTrials);
for neuron = 1:numNeurons
    for trial = 1:numTrials
        spikes = spTimes{neuron, trial};
        binnedSpikes(neuron, :, trial) = histcounts(spikes, [timeBins, trialDuration]);
    end
end

% Step 3: Reshape Data for Dimensionality Reduction
% Combine data across trials: (neurons x time bins) x trials
reshapedData = reshape(binnedSpikes, numNeurons * numTimeBins, numTrials)';

%% User choice: PCA, t-SNE, or UMAP
switch lower(opt.neurDyn.method)
    case 'pca' % PCA
        dataMean = mean(reshapedData, 1);
        centeredData = reshapedData - dataMean;
        [coeff, score, ~] = pca(centeredData);
        reducedData = score(:, 1:3);
        titleStr = 'Neural Dynamics: PCA';

    case 'tsne' % t-SNE
        rng('default');
        [reducedData, loss] = tsne(reshapedData, 'NumDimensions', 3, 'Perplexity', 10, "Standardize", true);
        titleStr = 'Neural Dynamics: t-SNE';

    case 'umap'  % UMAP (requires UMAP toolbox: https://www.mathworks.com/matlabcentral/fileexchange/71902)
        addpath('umap'); % Adjust path as needed
        [reducedData, ~, ~] = run_umap(reshapedData, 'n_components', 3);
        titleStr = 'Neural Dynamics: UMAP';

    otherwise
        error('Invalid method. Choose "PCA", "tSNE", or "UMAP".');
end

%% Visualization
%
figure;
scatter3(reducedData(:,1), reducedData(:,2), reducedData(:,3), 30, turbo(numTrials), 'filled');
xlabel('Dim 1'); ylabel('Dim 2'); zlabel('Dim 3');
title(titleStr);
grid on; view(3); colorbar;

%
cmap = turbo(numTrials);
figure;
for trial = 1:numTrials
    plot3(reducedData(trial,1), ...
          reducedData(trial,2), ...
          reducedData(trial,3), ...
          'Color', cmap(trial,:), 'LineWidth', 1);
    hold on;
end
xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
title('Neural Manifold Trajectories (each trial colored differently)');

%% Save figure
savefig([opt.neurDyn.method '_neural_dynamics.fig']);

end