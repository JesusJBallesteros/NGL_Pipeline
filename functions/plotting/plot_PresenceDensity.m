function plot_PresenceDensity(csvFile, KeepLabels, shape)
    % Generates presence density maps for bodypart positions with an 
    % enclosing hexagon overlay. Reads a resultant file from DLC analisys,
    % processes it to remove 'likelihood' columns (perhaps worth to keep?),
    % to rearrange the headers of labels and facilitate the extraction.
    % INPUTS
    %   csvFile:    full path to the CSV file containing bodypart coordinates
    %               and likelihood. ie: 'C:\...\...\csvFile.csv'
    %   bodyLabels: a logical array (ie [1 0 1 1 0]) where each integer
    %               corresponds to a label/bodypard, ordered as in the csvFile headers,
    %               that indicates which bodyLabels to keep and which to discard.
    %               Needs human inspection of the csv file to determine the 
    %               order of the labels.
    %               TODO: Read and find desired labels
    %
    %   shape:      string 'Hex', 'Plus' ... To create a pre-determined shape in the plot. 
    %
    % Jesus. 13.03.2026

    %% USE INPUTS
    % Read CSV fully
    raw = readcell(csvFile);

    % Force logical array.
    KeepLabels = logical(KeepLabels);

    %% Labels and data Extraction
    % Extract bodyLabels (placed in file's 2nd row)
    bodyLabels = raw(2,:);
    % Extract single column headers ('x' 'y' 'likelihood')
    coords = raw(3,:);

    % TODO:CHECK IF OTHER SETTINGS GENERATE EXTRA COLUMNS

    % Data starts at row 4. TODO: CHECK IF ALWAYS TRUE
    data = cell2mat(raw(4:end,:));

    % Remove every 3rd column (likelihood) TODO: GOOD TO TRASH?
    keepIdx = ~strcmp(coords,'likelihood');
    data = data(:, keepIdx);
    bodyLabels = bodyLabels(keepIdx);

    % Rename bodyLabels to unique labels (bodyLabel_x, bodyLabel_y)
    % BC in csv file, they have the same header name!
    newNames = cell(size(bodyLabels)); % Allocate space
    seen = containers.Map; % Prepare key-map for already seen labels
    % loop over headers
    for i = 1:numel(bodyLabels)
        bL = bodyLabels{i}; % get label
        if ~isKey(seen, bL) % already seen?
            seen(bL) = 1; % if not, flag 'first'
            suffix = 'x'; % as first appearance, it is X
        else
            seen(bL) = seen(bL) + 1; % if seen, flag 'second'
            suffix = 'y'; % as second appearance, it is Y
        end
        newNames{i} = [bL '_' suffix]; % Make new header eg. 'label_x'
    end
    bodyLabels_coord = newNames;

    % Check how many labels we will need to plot
    bL_Unique = unique(erase(bodyLabels, {'_x','_y'}));
    bL_Unique = bL_Unique(KeepLabels);
    N = numel(bL_Unique);

    %% From data, check limits and calculate bins for 2D histogram
    allX = data(:, endsWith(bodyLabels_coord,'_x'));
    allY = data(:, endsWith(bodyLabels_coord,'_y'));
    xlimData = [min(allX(:)), max(allX(:))];
    ylimData = [min(allY(:)), max(allY(:))];

    % Set bin size and bin edges
    binSize = 50;
    xEdges = xlimData(1):binSize:xlimData(2);
    yEdges = ylimData(1):binSize:ylimData(2);

    %% Prepare plotting
    % Set the center of the x/y map (depends on the camera x/y image center)
    % HARDCODED FOR NOW (only HEXAGON, PLUS camera will have a different center)
    cx = 650; % GENERALIZE to:? mean(allX(:),'omitnan');
    cy = 540; % GENERALIZE to:? mean(allY(:),'omitnan');

    % Compute radius as maximum distance from centroid.
    % HARDCODED FOR NOW (HEXAGON)
    r = 630; % to test: max(sqrt((allX(:)-cx).^2 + (allY(:)-cy).^2));
    
    if strcmp(shape, 'Hex')
        % Generate 6 equidistant vertices in a cicle
        theta = (0:6) * pi/3;
        
        % Offset to center, vertices' x/y coordinates 
        xHex = cx + r*cos(theta); 
        yHex = cy + r*sin(theta);

    elseif strcmp(shape, 'Plus')
        % Arm half-width (w) and arm length from center (L)
        % Chosen so all 8 outer corners lie exactly on circle of radius r:
        %   L^2 + w^2 = r^2,  
        % with ratio L:w = 3:1
        w = r / sqrt(10);       % arm half-width
        L = 3 * r / sqrt(10);   % arm length from center

        px = [L,  w,  w, -w, -w, -L, -L, -w, -w,  w,  w,  L,  L];
        py = [w,  w,  L,  L,  w,  w, -w, -w, -L, -L, -w, -w,  w];

        % Offset to center
        xPlus = px + cx;
        yPlus = py + cy;
    end
   
    %%  Prepare tiled layout
    t = tiledlayout(N/2, 2,'TileSpacing','compact','Padding','compact');
    
    % Add a title row
    title(t,'Presence Density per bodyLabel','FontSize',14)

    for i = 1:N
        % Extract coordinates for this bodypart
        x = data(:, strcmp(bodyLabels_coord,[bL_Unique{i} '_x']));
        y = data(:, strcmp(bodyLabels_coord,[bL_Unique{i} '_y']));

        % 2D histogram counts
        counts = histcounts2(x, y, xEdges, yEdges);

        % Plot counts
        nexttile;
            imHandlesCounts(i) = imagesc(xEdges, yEdges, counts');
            axis equal tight;
            set(gca,'YDir','reverse');
            hold on;
            
            if strcmp(shape, 'Hex')
                plot(xHex, yHex,'w-','LineWidth',2); % HEX SHAPE
            elseif strcmp(shape, 'Plus')
                plot(xPlus, yPlus, 'w-', 'LineWidth', 2); % PLUS SHAPE
            end

            colorbar;
            colormap turbo

            ylabel('Y (pixels)');
            if i == N, xlabel('X (pixels)'); end
            title([bL_Unique{i} ' - Counts']);
    end

    % Print to A4
    set(t, 'PaperOrientation', 'portrait');
    set(t, 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);
end
