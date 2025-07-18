function neuralDynamics = calculate_neural_dynamics(neurons, fireRate, trialdef, opt)
% Neural Dynamics Analysis using PCA, t-SNE, or UMAP
% This script generates synthetic spike data and allows the user to analyze
% neural population dynamics using either PCA, t-SNE, or UMAP.

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
            nSpikes = poissrnd(firingRate(neuron) * trialDuration);
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