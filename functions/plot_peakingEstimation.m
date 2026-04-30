function f = plot_peakingEstimation(accelData_trial, magnitude_vector_trial, accAxes_wf, valid_peaks_wf, valid_peaks, accAxes_avg, opt, info)

red_rate = 1000;
f = figure;
    fig = tiledlayout(f,'flow','TileSpacing','tight','Padding','tight');
        
    % sample filtered data
    r = randperm(length(accelData_trial),1);
    t0 = nexttile(fig, 1, [1 6]);
        plot(accelData_trial{r});
        hold on
        plot(magnitude_vector_trial{r}/red_rate, 'LineWidth', 2 , Color=[.7 .7 .7]); 
        str1 = sprintf('thr. (%i SD)', opt.SDs);
        yline(info.threshold/red_rate,'--', str1, FontSize=10);
        title('Downsampled, bias corrected, 1 KHz, filtered');
        ylabel('m/s^3'), xlabel('session time (ms)');
        box off
        set(t0, FontSize=14, LineWidth=1.5, TickDir='out');
    
    % Plot combined examples of 3-axis signal 
    r = randperm(size(accAxes_wf,3), min([size(accAxes_wf,3), 6]));
    for i = 1:length(r)
        nexttile(fig, i+6);
        plot(accAxes_wf(:,:,r(i))); 
        hold on
        plot(valid_peaks_wf(r(i),:)/red_rate, LineWidth=2, Color=[.7 .7 .7]);
        xlabel('ms'), xlim([0 100]), xticks(25:25:75), xticklabels(-25:25:25);
        ylim([-40 80]);
        if i==1, ylabel('m/s^2'), yticks(-20:20:80);
        else,   ylabel([]), yticks([]); end
        if i==3, title('examples X/Y/Z traces', FontSize=14, FontWeight='bold'); end
        yline(info.threshold/red_rate,'--');
        set(gca, FontSize=14, LineWidth=1.5, TickDir='in'); box off
    end
    clear i r
    
    % Plot individual axes as imagesc, overlay average waveform
    tx = nexttile(fig, 13, [1 2]);
        imagesc(squeeze(accAxes_wf(:,1,:))');
            colormap cool
        xlabel('ms'), xlim([20 80]), xticks(25:25:75), xticklabels(-25:25:75);
        ylabel('Trace #');
        box off
                
    ty = nexttile(fig, 14, [1 2]);
        imagesc(squeeze(accAxes_wf(:,2,:))'); hold on
            colormap cool
        xlabel('ms'), xlim([0 100]), xticks(25:25:75), xticklabels(-25:25:75);
        ylabel([]), yticks([]), yticklabels([]);
        title('All XYZ accel. traces at detection');        
        box off

    tz = nexttile(fig, 15, [1 2]);
        imagesc(squeeze(accAxes_wf(:,3,:))'); hold on
            colormap cool
        xlabel('ms'), xlim([20 80]), xticks(25:25:75), xticklabels(-25:25:75);
        ylabel([]), yticks([]), yticklabels([]);
        box off

        c = colorbar(tz,"eastoutside","Box","off");
            c.Label.String = 'm/s^2';

    set(tx, YDir='normal', FontSize=12, LineWidth=1.5, TickDir='in');
    set(ty, YDir='normal', FontSize=12, LineWidth=1.5, TickDir='in');
    set(tz, YDir='normal', FontSize=12, LineWidth=1.5, TickDir='in');

    % Plot random trial 3-axis raw data, with overlaid Acc Magnitude, around Event
    r = randperm(length(accelData_trial),10);
    t1 = nexttile(fig, 16, [2 6]);
        add_y = 1;
        for i = 1:length(r)
            plot(add_y + accelData_trial{r(i)}); hold on
            plot(add_y + magnitude_vector_trial{r(i)}/red_rate, 'LineWidth', .2 , Color=[.7 .7 .7]);
            % an asterisk on top of each labeled peak
            for w = 1:length(valid_peaks{r(i)}(:,2))
                peakt = valid_peaks{r(i)}(w,2)*opt.fs + opt.ev_minus*opt.fs;
                text(peakt, add_y + 10, "*", FontSize=20)
            end
            % % a red asterisk on the LAST peak before event 
            % peakt = last_validp{r(i)}*opt.fs + opt.ev_minus*opt.fs;
            % if ~isempty(peakt)
            %     text(peakt(1), add_y + 10, "*", FontSize=20, Color='r')      
            % end
            add_y = add_y + 30;
        end
        xline(opt.ev_minus*opt.fs, '-.', LineWidth=2, FontSize=12);
        plotops = struct('xlabel', {''}, 'ylabel', {'Accel and jerk @ trial#'}, ...
                     'xticks',  0:opt.fs:(opt.ev_minus+opt.ev_plus)*opt.fs,  'yticks', 1:30:331, ...
                     'xticklabels', -opt.ev_minus:1:opt.ev_plus, 'yticklabels', r);
        prettify(plotops);
        ylim([-15 301]);
        title('Accel. & jerk traces'); box off

        set(t1, FontSize=14, LineWidth=1.5, TickDir='out');
        xlim(t1, [0 (opt.ev_minus+opt.ev_plus)*opt.fs]);
    
    % Plot valid pecks raster around event of interest
    str1 = sprintf('%s events', opt.limitEvents{1});
    t2 = nexttile(fig, 28, [1 6]);
        % plot all events, even empty ones
        toplot = cellfun(@(x) x(:,2), valid_peaks, 'UniformOutput', false);
        plotRaster(toplot, 1);
            plotops = struct('xlabel', {''}, 'ylabel', {str1}, ...
                             'xticks', -opt.ev_minus:1:opt.ev_plus,  'yticks', 0:max(1, floor(size(toplot,1)/2)):size(toplot,1), ...
                             'xticklabels', -opt.ev_minus:1:opt.ev_plus, 'yticklabels', 0:max(1, floor(size(toplot,1)/2)):size(toplot,1));
            prettify(plotops);
        % clear empty events
        empty = cellfun(@(x) any(isnan(x)), toplot, 'UniformOutput', false); empty = cell2mat(empty);
        toplot(empty) = [];
        % Count detected pecks
        npeaks = cellfun(@(x) numel(x), toplot, 'UniformOutput', false);
        sumpeaks = int2str(sum([npeaks{:}]));
        xline(t2, 0, '-.', LineWidth=2, FontSize=12);
        title(sprintf('Estimated pecks: %s. Raster and rate', sumpeaks))

        set(t2, FontSize=14, LineWidth=1.5, TickDir='out');
        xlim(t2, [-opt.ev_minus opt.ev_plus]); 
    
    % Plot perievent peck rate
    t3 = nexttile(fig, 34, [1 6]);
        toplot = cellfun(@(x) x*opt.fs, toplot, 'UniformOutput', false);
        plotPSTH(toplot, opt.stpSize, opt.binSize, [-opt.ev_minus*opt.fs opt.ev_plus*opt.fs], opt.fs);
            plotops = struct('xlabel', {'time (s)'}, 'ylabel', {'Pecks/s'}, ...
                             'xticks', 0:(opt.binSize/opt.stpSize)*(opt.fs/opt.binSize):(opt.ev_minus+opt.ev_plus)*(opt.binSize/opt.stpSize)*(opt.fs/opt.binSize), ...
                             'yticks', 0:1:5, ...
                             'xticklabels', -opt.ev_minus:1:opt.ev_plus, ...
                             'yticklabels', 0:1:5);
            prettify(plotops);
            xline(t3, opt.ev_minus*(opt.binSize/opt.stpSize)*(opt.fs/opt.binSize),'-.', opt.limitEvents{1}, LineWidth=2, FontSize=16);

        set(t3, FontSize=14, LineWidth=1.5, TickDir='out');
        xlim(t3, [0 (opt.ev_minus+opt.ev_plus)*(opt.binSize/opt.stpSize)*(opt.fs/opt.binSize)]);
        ylim(t3, [0 10]); 
        % hide yaxes
        yticks(t3, []);

    % Insets in XYZ traces pcolor plots
    pause(4) % fig needs to render completely first
    locat = [0.77 0.6 0.2 0.38]; % left bottom width height, within the pcolor plots
        
    % X axis inset
    inset = get_insetPos(tx, locat);
    txi = axes('Parent', f, 'Position', inset);
        plot(txi, accAxes_avg(:,1), 'LineWidth', 2, 'Color','k');
        txi.Color = 'none';
        xlabel(txi,'ms'); xlim(txi,[0 100]); xticks(txi,25:25:75)
            xticklabels(txi,-25:25:75)
        ylabel(txi,'m/s^2')
        box(txi,'off')
    
    % Y axis inset
    inset = get_insetPos(ty, locat);
    tyi = axes('Parent', f, 'Position', inset);
        plot(tyi, accAxes_avg(:,2), 'LineWidth', 2, 'Color','k');
        tyi.Color = 'none';
        xlabel(tyi,'ms'); xlim(tyi,[0 100]); xticks(tyi,25:25:75)
            xticklabels(tyi,-25:25:75)
        ylabel(tyi,'m/s^2')
        box(tyi,'off')
    
    % Z axis inset
    inset = get_insetPos(tz, locat);
    tzi = axes('Parent', f, 'Position', inset);
        plot(tzi, accAxes_avg(:,3), 'LineWidth', 2, 'Color','k');
        tzi.Color = 'none';
        xlabel(tzi,'ms'); xlim(tzi,[0 100]); xticks(tzi,25:25:75)
            xticklabels(tzi,-25:25:75)
        ylabel(tzi,'m/s^2')
        box(tzi,'off')   
end

% Helper
function inset = get_insetPos(tile_ax, location)
    % position of the tile axes in figure coordinates
    pos = tile_ax.Position;
    
    % inset relative to tile
    inset = [pos(1)+location(1)*pos(3), ...
             pos(2)+location(2)*pos(4), ...
             location(3)*pos(3), ...
             location(4)*pos(4)];
end