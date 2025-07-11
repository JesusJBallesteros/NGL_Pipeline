function plot_superletsTFR_extintion_tbt(TFR, param, opt)
% TFR plotting function, optimized for extintion project showing data in a
% trial by trial fashion.
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
param.testname = 'ASL_TbT_';

% min-max frequency
minfreq = min(cellfun(@min, opt.freqInterest));
maxfreq = max(cellfun(@max, opt.freqInterest));
logticks = round(logspace(log10(minfreq), log10(maxfreq), 10));

%% Plot Configuration
cfgp = [];
 cfgp.interactive   = 'no';
 cfgp.colormap      = turbo;
 cfgp.xlim          = [-3 3];
 cfgp.ylim          = [minfreq maxfreq+1];
 cfgp.clim          = [0 10];
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

% titles for each panel (each half of each block) 
tlt = {'Acquisition Early' 'Acquisition Late' 'Extintion Early' 'Extintion Late'};
ntlt = 0;

%% Plot Baselined TFRs.
for b = 1:opt.blocks % Each block
    for i = 1:2 % early vs late
        ntlt = ntlt+1;

        % File name
        figtitl = [param.testname, opt.SavFileName, '_StimOn2_NS-FS_zscored_',tlt{ntlt}];

        % Check number of trials
        ntrials = size(TFR{1,b}.rem{i}.powspctrm, 1);

        % adaptive position of panels
        rows = floor(sqrt(ntrials));
        while mod(ntrials, rows) ~= 0
            rows = rows - 1;
        end
        if rows == 1 && ntrials > 5, rows = 3; end
        columns = ceil(ntrials/rows);

        fig = figure('Visible', param.visible); % create figure, assign visibility
        set(fig, 'Position', [0, 0, param.screen.width, param.screen.height]); % Set fig size
        
        % Inside panel, a vertical tiled layout with no spacing
        t = tiledlayout(rows, columns, "TileSpacing", "compact", 'Padding', 'compact');
        
        % Title and subtitle with N info
        title(t, tlt{ntlt});

        for p = 1:ntrials 
            ax = nexttile(t); % Generate next tile
            col = mod(p-1, columns) + 1;
            row = ceil(p / columns);

            % No colorbar, use 'graphic current axes' to plot, set ylim
            cfgp.colorbar  = 'no';
            cfgp.figure    = 'gca';
            cfgp.trials    = p;

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
            if col ~= 1
                ax.YColor = 'none';
                ax.YTick = [];  % opcional
            end

            % Hide X axis if not in last row
            if row ~= rows
                ax.XColor = 'none';
                ax.XTick = [];
            end
            
            % Individual title
            title(['trial ', num2str(p)]);  
        end
    
    %% Save the figure as .png
    if ~isfolder(fullfile(opt.analysis,'plots', 'TFR')), mkdir(fullfile(opt.analysis,'plots', 'TFR')); end
    exportgraphics(fig, fullfile('D:\S3_ExperimentArena\data\analysis\TFR', [figtitl, '.png']));

    close all hidden
    end
end

end