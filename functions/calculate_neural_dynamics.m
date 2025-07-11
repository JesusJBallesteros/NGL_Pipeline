function [projectedData, covMatrix, eigen] = calculate_neural_dynamics(spikeTimes, opt, input)
% Neural Manifold Construction from Spike Data

% Step 1: Set parameters and get data
% Initialize spike data structure
% spikeTimes = []; % input into function
numNeurons = 50;        % get from data
trialDuration = 10;     % get from data
numTrials = 50;        % get from data
binSize = 0.1;         % tell from ops

timeBins = 0:binSize:trialDuration-binSize; % Time bins constuction

% Step 2: Bin Spike Data
% Initialize binned spike matrix: neurons x time bins x trials
numTimeBins = length(timeBins);
binnedSpikes = zeros(numNeurons, numTimeBins, numTrials);

for neuron = 1:numNeurons
    for trial = 1:numTrials
        % Get spike times for this neuron and trial
        spikes = spikeTimes{neuron, trial};
        % Bin spikes into time bins
        binnedSpikes(neuron, :, trial) = histcounts(spikes, [timeBins, trialDuration]);
    end
end

% Step 3: Reshape Data for Dimensionality Reduction
% Combine data across trials: (neurons x time bins) x trials
reshapedData = reshape(binnedSpikes, numNeurons * numTimeBins, numTrials)';

% Step 4: Perform Principal Component Analysis (PCA)
% Subtract mean across trials
meanSubtractedData = reshapedData - mean(reshapedData, 1);
% Compute covariance matrix
covMatrix = cov(meanSubtractedData);

% Perform eigen decomposition
[eigen.Vectors, eigen.Values] = eig(covMatrix);
% Sort eigenvalues and eigenvectors in descending order
[eigen.Values, sortIdx] = sort(diag(eigen.Values), 'descend');
eigen.Vectors = eigen.Vectors(:, sortIdx);
% Project data onto the first three principal components
projectedData = meanSubtractedData * eigen.Vectors(:, 1:3);



end