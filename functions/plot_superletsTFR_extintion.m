function plot_superletsTFR_extintion(TFR, param, opt)
% 'Basic' TFR plotting function, although optimized for extintion project.
% It used the TFR data obtained from FieldTrip functionalities, preferably
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

% File name
figtitl = [param.testname, opt.SavFileName, '_TFR_StimOn2_NS-FS_zscored_',opt.TFRmethod];
% titles for each panel (each half of each block) 
tlt = {'Acquisition Early' 'Acquisition Late' 'Extintion Early' 'Extintion Late'};
% position of such panels
posit = [0 .5 .5 .5; .5 .5 .5 .5; 0 0 .5 .5; .5 0 .5 .5];

%% Plot Configuration
cfgp = [];
 cfgp.figure    = 'gca';
 cfgp.interactive   = 'no';
 cfgp.colormap      = turbo;
 cfgp.xlim          = [-3 3];
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
 
%% Plot Baselined TFRs.
fig = figure('Visible', param.visible); % create figure, assign visibility
set(fig, 'Position', [0, 0, param.screen.width, param.screen.height]); % Set fig size
p = 0; % reset panel count

for b = 1:opt.blocks % Each block
    for i = 1:2 % early vs late
        % Create panel as described in 'posit'
        p = p+1; % Next panel
        pn = uipanel(fig, 'Position', posit(p,:));

        % Inside panel, a vertical tiled layout with no spacing
        t = tiledlayout(pn, (size(opt.freqInterest,1)), 1, "TileSpacing", "none");
        
        % Loop over the frequency ranges
        for f = (size(opt.freqInterest,1)):-1:1 
            nexttile(t); % Generate next tile

            % No colorbar, use 'graphic current axes' to plot, set ylim
            cfgp.colorbar = 'no';
            cfgp.ylim      = [ceil(min(TFR{f,b}.rem{i}.freq)) ...
                              floor(max(TFR{f,b}.rem{i}.freq))];

            % Plot using FT capabilities
            try ft_singleplotTFR(cfgp, TFR{f,b}.rem{i}); hold on
            catch, continue, end % Next iteration if error
            
            % ticks out
            t.Children(1).TickDir = "out";
            
            % lines for approx. itiOn and exact StimOn times
            xline(0, '--w', 'LineWidth', 2);
            xline(-2, ':w', 'LineWidth', 2);
            
            % Adjust color scale limits
            if ~ischar(cfgp.clim)
                clim(cfgp.clim);
            end

            % superlet params annotation
            txt = {['c=', num2str(TFR{f,b}.rem{i}.cfg.width)], ...
                   ['o:', num2str(TFR{f,b}.rem{i}.cfg.order(1)) '-' num2str(TFR{f,b}.rem{i}.cfg.order(end))]};
            text(-2.9, cfgp.ylim(2)*0.9 , txt, 'FontSize', 8, 'Color', [1 1 1],'FontWeight','bold')
        end

        % shared X-axis adjustment
        xticks(cfgp.xlim(1):1:cfgp.xlim(2));
        xticklabels(cfgp.xlim(1):1:cfgp.xlim(2));
        xlabel(t, 'Time since Choice Stim. appear (s)'); 
        
        % shared y-axis label
        ylabel(t, 'Frequency (Hz)');
 
        % Title and subtitle with N info
        title(t, tlt{p}, sprintf('Ntrials = %d', size(TFR{f,b}.rem{i}.powspctrm,1)));

        % Colorbar adjustment
        if f == 1 % this makes the bar appear only once per plot
            cb = colorbar;
            cb.Layout.Tile = 'east'; % and this makes it be shared by all
            cb.Label.String = cfgp.colorbartext;
            cb.Label.FontSize = 12;
            cb.FontSize = 12;
            if ~ischar(cfgp.clim)
                cb.Ticks = cfgp.clim(1):1:cfgp.clim(2);
            end
        end
    end
end
    
%% Save the figure as .png
if ~isfolder(fullfile(opt.analysis,'plots', 'TFR')), mkdir(fullfile(opt.analysis,'plots', 'TFR')); end
saveas(fig, fullfile('D:\S3_ExperimentArena\data\analysis\TFR', figtitl), 'png');
% exportapp(fig, fullfile(opt.analysis,'plots', 'TFR', [figtitl '.png']));
close all hidden

end