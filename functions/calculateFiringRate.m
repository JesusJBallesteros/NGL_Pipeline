function [firing_rate, time_vector] = calculateFiringRate(timestamps, sample_rate, sigma)
    % calculateFiringRate - Calculate firing rate from timestamps with Gaussian smoothing
    %
    % Inputs:
    %   timestamps   - Array of spike times in seconds
    %   sample_rate  - Sampling rate in Hz (samples per second)
    %   sigma        - Standard deviation of the Gaussian kernel in seconds
    %
    % Outputs:
    %   time_vector  - Array of time points corresponding to the firing rate
    %   firing_rate  - Array of firing rates (spikes per second) at each time point
        
    % Define the time vector
    t_min = min(timestamps);        % Start time (in seconds)
    t_max = max(timestamps);        % End time (in seconds)
    time_vector = t_min:1/sample_rate:t_max; % Time bins based on the sample rate

    % Create the binary spike train
    spike_train = zeros(size(time_vector)); % Initialize the spike train
    
    % Map timestamps to the nearest time bins
    spike_indices = round(timestamps * sample_rate); % Convert timestamps to bin indices
    spike_indices = spike_indices(spike_indices > t_min*sample_rate & spike_indices <= time_vector(end)*sample_rate); % Bounds check
    spike_indices = spike_indices - spike_indices(1)+1;
    spike_train(spike_indices) = 1; % Mark spikes in the spike train    
    
    % Generate Gaussian kernel
    kernel_width = sigma * sample_rate; % Convert sigma from seconds to bins
    x = -3*kernel_width:3*kernel_width; % Kernel range (±3 sigma)
    gaussian_kernel = exp(-0.5 * (x / kernel_width).^2); % Gaussian formula
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel); % Normalize
    
    % Convolve spike train with Gaussian kernel
    smoothed_spike_train = conv(spike_train, gaussian_kernel, 'same');
    
    % Convert to firing rate
    firing_rate = smoothed_spike_train * sample_rate; % Scale to spikes per second
end
