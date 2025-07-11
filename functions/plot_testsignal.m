function plot_testsignal(FT_data,ch)
    if ~exist('ch','var')
        ch = input.test_ch;
    end
    addpath(genpath('toolboxes\multitaper_prerau'))
    
    % Parameters for MultiTapered FFT
    Fs              = FT_data.fsample; % double - sampling frequency (Hz)
    frequency_range = [0 100]; % [<min frequency>, <max frequency>]
    taper_params    = [2 3];   % [<TW>, <N>]
                               % Time-half bandwidth product (TW) is N*BW/2 where
                               %  N is window length (s), BW is main lobe bandwidth
                               % optimal number of tapers (N) is 2*TW-1
    window_params   = [0.5 .1];  % [window length, step size] (s)
    min_NFFT        = 0;       % double - minimum allowable NFFT size, adds zero 
                               %  padding for interpolation (closest 2^x)
    detrend_opt     = 'linear';% string - detrend data window ('linear' (def.), 'constant', 'off')
    weighting       = 'unity'; % string - weighting of tapers ('unity' (def.), 'eigen', 'adapt')

    plot_on         = false;   % boolean - plot results
    verbose         = false;   % boolean - display spectrogram properties

%     printsettings   = '800x300_hdpi_spectrum'; % string, name of required 'style'
                               %  to apply to the plots. See 'sdf' function for more.
    
    % Check if channel requirements are possible
    if ch(end) > size(FT_data.trial{1,1},1)
        ch = 1:size(FT_data.trial{1,1},1);
        disp('The channel list to plot has been modified because the requested list was not possible')
    end
    
    for i=ch
        data = FT_data.trial{1,1}(i,:)'; % samples x 1 vector - time series data
           
        % COMPUTE MULTITAPER SPECTROGRAM
        [spectral.data, spectral.t, spectral.f] = ...
            multitaper_spectrogram_mex(data, Fs, frequency_range, ...
                                       taper_params, window_params, min_NFFT, ...
                                       detrend_opt, weighting, plot_on, verbose);
        
        % Obtain max/min values (in dB)
        cmax = max(max(mag2db(spectral.data)));
        cmin = min(min(mag2db(spectral.data)));
        
        % PLOT MTSpectrogram, in dB
        figure(1), pcolor(spectral.t, spectral.f, mag2db(spectral.data)); hold on;
            shading interp; 
            colormap("turbo"); % improved 'jet'
            set(gca, 'clim', [cmin*.25 cmax*.95]);  % Adjust dB scale
            c = colorbar('location','eastoutside'); % set scale location
            c.Label.String = 'dB'; % scale label
            
            % plot voltage trace, scaled to fit on the low frequency range
            plot(FT_data.time{1}, ((data)/50)+35, 'Color', 'w', 'LineWidth', 1)
            
            ylabel('Frequency (Hz) and Voltage (uV/50)', 'fontsize', 12); % graph label
            xlabel('sec', 'fontsize', 12); % graph label
            xlim([15 45]); ylim([0 50]);   % graph limits

        title(sprintf('Spectrogram for Ch: %s',FT_data.label{i,1}));
%         sdf(1, printsettings);
        box('off'); 
        saveas(gca, sprintf('testsignal_ch%s',FT_data.label{i,1}), 'png')
        close all
    end
end