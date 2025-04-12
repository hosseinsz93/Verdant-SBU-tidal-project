function VelPlot_Verdant(velMag, velDir, seadepths, units, DMY, workfolder, deploymentName)
    %% VELPLOT_VERDANT - Creates velocity profile plots from ADCP data
    %
    % By: Jonathan Craig
    % Date: 26 February, 2025
    %
    % Last modified: Hossein Seyedzadeh
    % Date: 10 April, 2025
    %
    % Syntax:
    %   VelPlot_Verdant(velMag, velDir, seadepths, units, DMY, workfolder)
    %   VelPlot_Verdant(velMag, velDir, seadepths, units, DMY, workfolder, deploymentName)
    %
    % Inputs:
    %   velMag        - Velocity magnitude values (1D array)
    %   velDir        - Velocity direction values in degrees (1D array)
    %   seadepths     - Depth values for each measurement (1D array)
    %   units         - 2-element array with [depthUnit, velocityUnit]
    %   DMY           - Datetime value for this profile
    %   workfolder    - Output folder for saving plots
    %   deploymentName- (Optional) Name for deployment identification (default: "Default")
    %
    % Outputs:
    %   None - Function saves three plot types to output folders:
    %     1. Speed profiles (2D) - Speed vs. depth
    %     2. Direction profiles - Compass plot showing current direction
    %     3. 3D velocity profiles - 3D quiver plot showing U,V components
    %
    % Example:
    %   VelPlot_Verdant(velMag, velDir, depths, ["m","m/s"], datetime('now'), pwd, "Deployment1")
    
    %% Input validation
    arguments
        velMag (1,:) {mustBeNumeric}
        velDir (1,:) {mustBeNumeric}
        seadepths (1,:) {mustBeNumeric}
        units (1,2) {mustBeText}
        DMY datetime
        workfolder {mustBeFolder}
        deploymentName = "Default" % Parameter for deployment identification
    end
    
    %% Extract and prepare data
    % Get units for plotting
    unitLength = units(1,1);
    unitVel = units(1,2);
    
    % Find any NaN values in the data
    velMag_NaN = find(ismissing(velMag));
    velDir_NaN = find(ismissing(velDir));
    
    %% Calculate velocity components
    % Convert magnitude and direction to cartesian components
    U = velMag.*sind(velDir); % East component of velocity
    V = velMag.*cosd(velDir); % North component of velocity
    W = zeros(size(U));       % Vertical component (zero for ADCP measurements)
    
    % Calculate maximum velocity for plot scaling
    maxvel = max(velMag,[],'all');
    
    %% Prepare depth values for plotting
    seadepths = sort(seadepths,'ascend');
    depths_NaN = seadepths(velMag_NaN);
    
    % Determine sea surface and bottom based on depth values
    if  all(seadepths > 0, 'all')
        seatop = max(seadepths);
        seabottom = 0;
    elseif all(seadepths < 0, 'all')
        seatop = 0;
        seabottom = min(seadepths);
    else
        seatop = max(seadepths);
        seabottom = min(seadepths);
    end
    
    %% Define output folders for different plot types
    % Get the profile directories for this deployment
    speedFolder = fullfile(workfolder, "SpeedProfiles");
    dirFolder = fullfile(workfolder, "DirProfiles");
    velFolder = fullfile(workfolder, "VelProfiles");
    
    %% PLOT TYPE 1: Speed Profiles (2D)
    % Create 2D velocity magnitude vs. depth plots
    for i = 1:length(DMY)
        % Create figure with specified name, but don't display it
        fig = figure('Name', sprintf("%s: ADCP at Orient Point on %s", deploymentName, DMY(i)), 'visible', 'off');
        
        % Plot velocity magnitude vs. depth
        plot(velMag(i,:), seadepths);
        hold on
        axis([0 maxvel seabottom seatop]);
        
        % Add title and labels
        title(sprintf('2D Vertical Velocity Profile - %s - %s', deploymentName, DMY(i)));
        ylabel(sprintf("Depth [%s]", unitLength));
        xlabel(sprintf("Current speed [%s]", unitVel));
        
        % Save the plot to file
        filename = sprintf("%s_2D_VVP_%s.png", deploymentName, string(DMY(i)));
        filename = replace(filename, ':', '-');
        saveas(gcf, fullfile(speedFolder, filename));
        clf
        close
    end
    
    %% PLOT TYPE 2: Direction Profiles (Compass)
    % Calculate mean velocity components for the average vector
    Ucompass = velMag.*sind(velDir);
    Vcompass = velMag.*cosd(velDir);
    Uavg = mean(Ucompass,2);
    Vavg = mean(Vcompass,2);
    
    for i = 1:length(DMY)
        % Create figure with specified name, but don't display it
        fig = figure('Name', sprintf("%s: Direction at Orient Point on %s", deploymentName, DMY(i)), 'visible', 'off');
        
        % Create compass plot for each depth bin
        c1 = compass(Ucompass(i,:), Vcompass(i,:), 'r');
        hold on
        
        % Add the average direction as a thicker green vector
        c2 = compass(Uavg(i), Vavg(i));
        avgdir = c2(1);
        avgdir.LineWidth = 5;
        avgdir.Color = "green";
    
        % Remove only compass labels while preserving the title
        titleObj = get(gca, 'Title');
        titleText = get(titleObj, 'String');
    
        % Selectively remove only compass labels
        ax = gca;
        allTextObjs = findall(ax, 'Type', 'Text');
        for j = 1:length(allTextObjs)
            % Skip the title object
            if ~isequal(allTextObjs(j), titleObj)
                delete(allTextObjs(j));
            end
        end
    
        %% Add custom compass labels with cardinal directions
        compassRadius = maxvel * 1.2;
    
        % Loop through all directions in 30-degree increments
        for angle = 0:30:330
            % Calculate position on the circle
            x = compassRadius * sind(angle);
            y = compassRadius * cosd(angle);
            
            % Determine label text based on angle
            switch angle
                case 0
                    labelText = 'N (0° T)';
                    hAlign = 'center';
                    vAlign = 'bottom';
                case 90
                    labelText = 'E (90° T)';
                    hAlign = 'left';
                    vAlign = 'middle';
                case 180
                    labelText = 'S (180° T)';
                    hAlign = 'center';
                    vAlign = 'top';
                case 270
                    labelText = 'W (270° T)';
                    hAlign = 'right';
                    vAlign = 'middle';
                otherwise
                    labelText = sprintf('%d° T', angle);
                    
                    % Adjust alignment based on position in the circle
                    if angle < 90
                        hAlign = 'left';
                        vAlign = 'bottom';
                    elseif angle < 180
                        hAlign = 'left';
                        vAlign = 'top';
                    elseif angle < 270
                        hAlign = 'right';
                        vAlign = 'top';
                    else
                        hAlign = 'right';
                        vAlign = 'bottom';
                    end
            end
            
            % Add the text label
            text(x, y, labelText, 'HorizontalAlignment', hAlign, 'VerticalAlignment', vAlign, 'FontSize', 8);
        end
    
        % Add title and save the plot
        title(sprintf('2D Velocity Compass - %s - %s', deploymentName, DMY(i)));
        filename = sprintf("%s_2D_VC_%s.png", deploymentName, string(DMY(i)));
        filename = replace(filename, ':', '-');
        saveas(gcf, fullfile(dirFolder, filename));
        clf
        close
    end
    
    %% PLOT TYPE 3: 3D Velocity Profiles (Quiver)
    % Define axis limits for 3D plot
    umin = -1*round(max(velMag,[],'all'),2);
    umax = round(max(velMag,[],'all'),2);
    vmin = -1*round(max(velMag,[],'all'),2);
    vmax = round(max(velMag,[],'all'),2);
    zmin = 0;
    zmax = seatop;
    
    for i = 1:length(DMY)
        % Create larger figure for 3D plot
        fig = figure('Name', sprintf("%s: 3D Profile at Orient Point on %s", deploymentName, DMY(i)), 'visible', 'off', 'Position', [10 10 900 600]);
        
        % Create 3D quiver plot with velocity vectors
        scalfac = 0.25;  % Scale factor for arrows
        quiver3(W(i,:), W(i,:), seadepths, U(i,:), V(i,:), W(i,:), scalfac, 'r');
        hold on
        
        % Add station location information
        annotation('textbox', [0.75 0.05 0.2 0.08], 'String', "Latitude = 41° 152993' N Longitude = -72° 350476' W");
        
        % Set axis limits and adjust plot appearance
        axis([umin*2 umax*2 vmin*2 vmax*2 zmin zmax]);
        axscale = 0.8;
        os = [0.0 0.05];
        ax = gca;
        hax = gca;
        sc = 0.5*(1-axscale)*[1 1];
        set(hax, 'position', [sc+os 1-2*sc])
        
        % Add title and axis labels
        title(sprintf('3D Vertical Velocity Profile - %s - %s', deploymentName, DMY(i)));
        zlabel(sprintf("Depth [%s]", unitLength));
        ylabel(sprintf("V [%s]", unitVel));
        xlabel(sprintf("U [%s]", unitVel));
        
        % Save the 3D plot to file
        filename = sprintf("%s_3D_VVP_%s.png", deploymentName, string(DMY(i)));
        filename = replace(filename, ':', '-');
        saveas(gcf, fullfile(velFolder, filename));
        clf
        close
    end
end