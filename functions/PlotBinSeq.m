function [fig]= PlotBinSeq(seqID, seqData, varargin)
%% Takes all the data for each sequence and plots an histogram seen from above
% inputs: 
% seqID   - is a numerical value, sequence label -> relevant for the title
% seqData - camera tracking data -> matrix of position coordinates
%           (x,y,frame number) - data collected for all trials of seqID
%           OR, blob tracking (x,y) position
% varargin - optional argument containing parameters to modify the
%           histogram2 properties:
%             param.title     = 'text to display as title '; 
%             param.subtitle  = 'extra text' % here, used to further name saved files 
%             param.NumBins   = [76 76];  % number of desired bins
%             param.BinWidth  = [10 10];  % size of bins
%             param.shape     = [0 256; 256 256; ...; 256 530; 0 530]; shape defined as vertices (x,y) connected by lines
%             param.visible   = 'off';    % show or not the figure
%             param.AXvisible = 'off';    % show or not the grid and the axes
%             param.colmap    = 'hot';    % colormap of choice
%             param.flipY     = 0;        % if true, reverses Y axis to set 0 on top.
%             param.norm      = 'count';
%             param.barunit   = title the colorbar
%             param.savepath  = specify saving destiny
%             param.posspk    = special input to generate a histogram using hist2 data
% Outputs:
% hist2   - The 2 Dimensional distribution data created for the figure, to
%           use in further calculations
%
%   TODO: perhaps include seqID as varargin, or make it go as subtitle?
%
% Sara: version 1: 
%      completed: 27.04.2022
%      some comments: 08.01.2025
% Jesus, added: 24.01.2025
%       possibility for optional argument input 'param', 
%       checks for 'seqID' and 'seqData', for generalization,
%       flipY parameter, to rotate the axes so origin is top-left
%       'f' is not manually defined, but based on 'v',
%       X/Y BinLimits depend now on the shape given,
%       option to change the normalization of data ('param.norm')
%        and to title the color bar based on it
%       saving happens inside function
%       figure visibility is optional
%       colorbar maximum depends on histogram max value
%       Subtitle to further name individual files (case of neurons/sequences per animal)
%       colormap and colorbar lomits ('CLim') are choosen depending on the
%         data, for now. * In some case it might be good to be fixed, as if normalizations are done... ? 

%% Default param,
if nargin < 3
    % if not given as an additional input, it should fit the travelling 
    % pigeon paradigms representation of sequences (TO BE TESTED).
    %  Until Sara's approval at least.
    param.title     = 'Sequence ';  %
    param.subtitle  = '';           % Except if interesting for saving, empty
    param.NumBins   = [76 76];      % number of desired bins
    param.BinWidth  = [10 10];      % size of bins
    param.shape     = [0 256; 256 256; 256 0; 530 0; 530 256; 768 256; 768 530; 530 530; 530 768; 256 768; 256 530; 0 530];
    param.visible   = 'on';         % see the figure
    param.AXvisible   = 'off';      % hide axes
    param.colmap    = 'hot';        % colormap
    param.flipy     = 0;        % Do not Reverse Y axis
    param.norm      = 'count';
    param.barunit   = 'occurrences in bin';
    param.savepath  = pwd;

else % if given just take them all
    param = varargin{1};
    if strcmpi(param.norm, 'probability') || strcmpi(param.norm, 'pdf')
        param.barunit   = 'prob of occurrence';
    elseif strcmpi(param.norm, 'count')
        param.barunit   = 'occurrences in bin';
    end        
end

%% Check obligatory input arguments
if ~isempty(seqID) % If sequence is explicitly left empty
    param.title = [param.title, num2str(seqID)]; % use sequence in title
end

if size(seqData,2) < 3 % if position data has only (x,y) values
    JoinAllPlotsXSorted = seqData(:,1); % set as x
    JoinAllPlotsYSorted = seqData(:,2); % set as y
else % if data has 3 dimensions
    [~,SortOrder] = sort(seqData(:,3));         %
    JoinAllPlotsXSorted = seqData(SortOrder,1); %
    JoinAllPlotsYSorted = seqData(SortOrder,2); %
end
    
%% Create the histogram2 plot
fig = figure('visible', param.visible); % visibility; 'painters' allow for transparency
    
    PlotBin = histogram2(JoinAllPlotsXSorted, JoinAllPlotsYSorted, ... % X and Y data goes here
                         'DisplayStyle', 'tile', 'ShowEmptyBins', 'on', 'EdgeColor','none'); % Properties to leave fixed
    
    % Customizable properties for histogram2
    PlotBin.NumBins     = param.NumBins;
    PlotBin.XBinLimits  = [0 max(max(param.shape))]*1.1; % 10% bigger than shape border
    PlotBin.YBinLimits  = [0 max(max(param.shape))]*1.1;
    PlotBin.BinWidth    = param.BinWidth;
    PlotBin.Normalization = param.norm;
    axis([PlotBin.XBinLimits, PlotBin.YBinLimits])
    colormap(param.colmap)

    % Get CLim
    if strcmpi(param.norm, 'probability') || strcmpi(param.norm, 'pdf')
        param.limClim(2) = round(max(max(PlotBin.Values/2)),2); % obtain Max value/2 for colorbar
    elseif strcmpi(param.norm, 'count')
        param.limClim(2) = round(max(max(PlotBin.Values)/4)); % obtain Max value/4 for colorbar
    end
    title(param.title);
    
    % Clim determination
    c = gcf;
    c.Children.CLim(2) = param.limClim(2);

    % Main axes
    ax = gca;
    ax.XAxis.Visible    = param.AXvisible;
    ax.YAxis.Visible    = param.AXvisible;
    if param.flipY % Set (0,0) on top left
        set(ax,'XDir','normal')
        set(ax,'YDir','reverse')
    end

    % Color bar creation
    c = colorbar('Ticks', param.limClim, 'TickLabels', {param.limClim(1), param.limClim(2)}, 'TickLength', 0, 'FontSize', 10);
        c.Label.String = param.barunit;
        c.Label.Position = [1.25 param.limClim(2)/2 0];
    hold on

    % this creates the shape of the arena 
    v = param.shape;
    f = 1:size(param.shape,1);
    patch('Faces', f , 'Vertices' , v, ...
          'EdgeColor', 'white', 'FaceColor', 'none', 'LineWidth', 1);
    
%% each bin is 1 pixel 
%     mat = zeros(768,768);
%     for i=1:length(JoinAllPlotsXSorted)
%         mat(round(JoinAllPlotsYSorted(i),0),round(JoinAllPlotsXSorted(i),0)) = ...
%             mat(round(JoinAllPlotsYSorted(i),0),round(JoinAllPlotsXSorted(i),0))+1;
%     end
%     figure
%     imshow(mat)
%     PlotBin = image(mat/max(mat(:)),'CDataMapping','scaled');
%     axis square
%     colormap('hot')
%     colorbar;
%     f = gcf;
%     f.Children(2).CLim = [0 0.05];
% 
%     hold on
%     v = [0 256; 256 256; 256 0; 530 0; 530 256; 768 256; 768 530; 530 530; 530 768; 256 768; 256 530; 0 530];
%     f = [1 2 3 4 5 6 7 8 9 10 11 12];
%     patch('Faces',f,'Vertices',v,...
%         'EdgeColor','white','FaceColor','none','LineWidth',1);
%     title(['Date', num2str(sessionDate),' ', 'trials:', num2str(trlNoID(find(trlNoID)))]);
%
%% Save
% Check/Create folder.
if ~exist(param.savepath,"dir"), mkdir(param.savepath), end
exportgraphics(fig, fullfile(param.savepath, [param.title, ' ', [param.subtitle{:}], '.png']), 'Resolution', 300);
close all
end 
