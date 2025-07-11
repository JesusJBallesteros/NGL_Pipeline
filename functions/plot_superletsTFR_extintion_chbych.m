function plot_superletsTFR_extintion_chbych(TFR, param, opt)
% TFR plotting function, optimized for extintion project showing data in a
% channel by channel fashion.
%
% It uses the TFR data obtained from FieldTrip functionalities, preferably
% via 'superlet' method, and plots average power across channels and trials.
% The power data (either absolute, relative or substraction) will be
% baselined to a givent amount of time and represented as a Z-value
% relative scale.
%
% Jesus 15.04.2025

%% Figure size, name and layouts
if strcmpi('adaptive', param.size)
    param.screen.size   = get(0, 'ScreenSize');  
    param.screen.width  = param.screen.size(3);
    param.screen.height = param.screen.size(4);
else
    param.screen.width  = param.size(1); 
    param.screen.height = param.size(2);
end
param.testname = 'ASL_CbC_';

%% min-max frequency
minfreq = min(cellfun(@min, opt.freqInterest));
maxfreq = max(cellfun(@max, opt.freqInterest));
logticks = round(logspace(log10(minfreq), log10(maxfreq), 5));

%% Check number of channels
chanMap = [];
nchans = opt.numChannels;
opt.chanMap = fullfile(opt.analysisCode, opt.KSchanMapFile);
chanMap = load(opt.chanMap);

% Add personalized labels
chanMap.label = num2cell(chanMap.chanMap);
chanMap.label = cellfun(@num2str, chanMap.label, 'UniformOutput', false);
for ch=1:nchans
    chanMap.label{ch} = sprintf('%03s', chanMap.label{ch});
end

% adaptive position of panels
cols = max(chanMap.kcoords);
while mod(nchans, cols) ~= 0
    cols = cols - 1;
end
rows = ceil(nchans/cols);

%% Plot Configuration
cfgp = [];
 cfgp.interactive   = 'no';
 cfgp.colormap      = turbo;
 cfgp.xlim          = [-3 3];
 cfgp.ylim          = [minfreq maxfreq+1];
 cfgp.clim          = [0 3];
 cfgp.interactive   = 'no';
 cfgp.fontsize      = 14;
 cfgp.title         = ' ';
 cfgp.baseline      = [-3 0]; % 'no'; %
 cfgp.baselinetype  = 'zscore';
 % cfgp.maskstyle     = 'saturation';
 if ischar(cfgp.baseline)
    cfgp.colorbartext  = 'NS-FS Power (no-baselined, abs)';
 else
    cfgp.colorbartext  = ['NS-FS Power (baselined, ', cfgp.baselinetype, ')'];
 end  
% No colorbar, use 'graphic current axes' to plot
cfgp.colorbar  = 'no';
cfgp.figure    = 'gca';

% titles for each panel (each half of each block) 
tlt = {'Acquisition Early' 'Acquisition Late' 'Extintion Early' 'Extintion Late'};
ntlt = 0;

%% Plot Baselined TFRs.
for b = 1:opt.blocks % Each block
    for i = 1:2 % early vs late
        ntlt = ntlt+1;
        colt = 0:rows:nchans-1;

        % File name
        figtitl = [param.testname, opt.SavFileName, '_StimOn2_NS-FS_zscored_',tlt{ntlt}];

        fig = figure('Visible', param.visible); % create figure, assign visibility
        set(fig, 'Position', [0, 0, param.screen.height, param.screen.width]); % Set fig size
        
        % Inside panel, a vertical tiled layout with no spacing
        t = tiledlayout(rows, cols, "TileSpacing", "tight", 'Padding', 'loose');
        
        % Title and subtitle with N info
        title(t, tlt{ntlt});

        for p = 1:nchans
            iseven = rem(p,2) == 0;

            % find where it shoulb be located
            if ~iseven % odds, for two shank
                colt(1) = colt(1)+1;
                cfgp.channel   = chanMap.label{colt(1)};
            else  % even, for two shank
                colt(2) = colt(2)+1;
                cfgp.channel = chanMap.label{colt(2)};
            end

            ax = nexttile(t); % Generate next tile

            % Plot using FT capabilities
            for f = 3:-1:1
                ft_singleplotTFR(cfgp, TFR{f,b}.rem{i}); hold on
            end
            
            % ticks out
            t.Children(1).TickDir = "out";
            
            % lines for approx. itiOn and exact StimOn times
            xline(0, '--w', 'LineWidth', 2);
            xline(-2, ':w', 'LineWidth', 2);
            
            % Adjust color scale limits
            if ~ischar(cfgp.clim)
                clim(cfgp.clim);
            end
            
            set(gca, 'YScale', 'log')

            % shared X-axis adjustment
            xticks(cfgp.xlim(1):1:cfgp.xlim(2));
            xticklabels(cfgp.xlim(1):1:cfgp.xlim(2));
            xlabel(t, 'Time since Choice Stim. appear (s)'); 
            
            % shared y-axis label
            ylabel(t, 'Frequency (Hz)');
            yticks(logticks);
            yticklabels(logticks);
            
            % Hide y axis if not in first column
            if mod(p,cols)==0
                ax.YColor = 'none';
                ax.YTick = [];  % opcional
            end

            % Hide X axis if not in last row
            if p <= nchans-cols
                ax.XColor = 'none';
                ax.XTick = [];
            end
            
            % Individual title
            title(cfgp.channel, "FontSize", 6);
        end
    
    %% Save the figure as .png
    if ~isfolder(fullfile(opt.analysis,'plots', 'TFR')), mkdir(fullfile(opt.analysis,'plots', 'TFR')); end
    exportgraphics(fig, fullfile('D:\S3_ExperimentArena\data\analysis\TFR', [figtitl, '.png']));

    close all hidden
    end
end

end