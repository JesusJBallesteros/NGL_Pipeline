function drawChanMap(input, mapfile, varargin)
% This function allows to visualize the channel map file used during
% this recording. It allows to change the status of individual channels
% between 'connected' and 'disconnected'.
%
% Jesus 21.12.2023

% If left empty, return to invoker function
if strcmp(mapfile,'')
    disp('No map file especified.')  
    return
end

% Check if not valid channels are given
if nargin>2
    selectCh = varargin{1};
else
    % Empty selectCh
    selectCh = [];
end

% If valid, load the channel map file
if isfile(fullfile(input.analysisCode,mapfile))
    load(fullfile(input.analysisCode,mapfile));
    
    % Default colors for connected/disconnected channels
    bckgnd.On  = [0.00  1.00  1.00]; % cyan
    bckgnd.Off = [0.90  0.60  0.60]; % red
    
    if ~isempty(selectCh)
        connected(chanMap0ind==selectCh) = 0;
    end

    % Figure at position x,y with size w,h: [x y w h]
    figure("Position", [150 150 400 800]),
    
    % Plot layout of connected channels
    scatter(xcoords(connected), ycoords(connected), 250, ...
        'filled', 'square', 'MarkerFaceColor', bckgnd.On, 'MarkerEdgeColor', 'black')
    hold on

    % Plot layout of disconnected channels
    scatter(xcoords(~connected), ycoords(~connected), 250, ...
        'filled', 'square', 'MarkerFaceColor', bckgnd.Off, 'MarkerEdgeColor', 'black')
    hold on

    % Set axes limits, remove ticks and box
    xlim([min(xcoords)-50 max(xcoords)+50]);
    ylim([min(ycoords)-50 max(ycoords)+50]);
    xticks([]);
    yticks([]);
    box off

    % Label contacts with 0-indexed channel numbers
    for i=1:length(chanMap)
        txt = int2str(chanMap0ind(i));
        text(xcoords(i),ycoords(i)+1,txt,"FontWeight","bold", ...
            "VerticalAlignment","middle","HorizontalAlignment","center")
    end

else
    disp('Map file does not exist.')  
    return    
end


end
