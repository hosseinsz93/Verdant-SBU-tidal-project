function flowVisualized(mainOutputDir,DMY,velSigned_depthAvg,velDir_depthAvg,floodSign_depthAvg,ebbSign_depthAvg,Eas_depthAvg,Nor_depthAvg,floodVelMag_depthtimeAvg,ebbVelMag_depthtimeAvg,floodMeanDir,ebbMeanDir,siteID)

% Original author: Hossein Seyedzadeh
% Date modified: 18 May, 2026

% Check whether output folder exists
if ~exist(mainOutputDir, 'dir')
    error("The flow visualization function does not receive a valid main output directory.");
end

%% Figure 1: Summary of Info about Tidal Flow
    figure('Name',sprintf('%s (%s to %s)',siteID,DMY(1),DMY(end)),'Position',[50, 50, 1000, 450]);
    ax1 = axes('Position', [0.1 0.15 0.35 0.7]);
    bar([floodVelMag_depthtimeAvg, ebbVelMag_depthtimeAvg]);
    set(ax1, 'XTickLabel', {'Flood', 'Ebb'}, 'FontName', 'Times New Roman');
    ylabel('Average Magnitude (m/s)', 'FontName', 'Times New Roman');
    title(sprintf('%s (%s to %s)',siteID,DMY(1),DMY(end)), 'FontName', 'Times New Roman');
    grid on;  

    ax2 = polaraxes('Position', [0.55 0.1 0.4 0.8]);
    hold on;

    ax2.ThetaZeroLocation = 'top';
    ax2.ThetaDir = 'clockwise';
    ax2.FontName = 'Times New Roman';

    floodRad = deg2rad(floodMeanDir);
    ebbRad = deg2rad(ebbMeanDir);

    polarplot([0 floodRad], [0 floodVelMag_depthtimeAvg], 'LineWidth', 2, 'Color', [0.2 0.4 0.9]);
    polarplot([0 ebbRad], [0 ebbVelMag_depthtimeAvg], 'LineWidth', 2, 'Color', [0.9 0.2 0.2]);

    polarscatter(floodRad, floodVelMag_depthtimeAvg, 70, [0.2 0.4 0.9], 'filled');
    polarscatter(ebbRad, ebbVelMag_depthtimeAvg, 70, [0.9 0.2 0.2], 'filled');

    thetaticks(0:45:315);
    thetaticklabels({'N','NE','E','SE','S','SW','W','NW'});

    text(floodRad, floodVelMag_depthtimeAvg*1.15, sprintf('%.1f° T', floodMeanDir), 'Color', [0.2 0.4 0.9], 'FontWeight', 'bold');
    text(ebbRad, ebbVelMag_depthtimeAvg*1.15, sprintf('%.1f° T', ebbMeanDir), 'Color', [0.9 0.2 0.2], 'FontWeight', 'bold');

    legend({'Flood', 'Ebb', '', ''}, 'Location', 'southoutside', 'FontName', 'Times New Roman');
    title(sprintf('%s (%s to %s): Dominant Current Directions (True North)',siteID,DMY(1),DMY(end)), 'FontName', 'Times New Roman');
    
    saveas(gcf, fullfile(mainOutputDir, '_FlowSummary.png'));
    saveas(gcf, fullfile(mainOutputDir, '_FlowSummary.fig'));
    
    %% Figure 2: Time Series
    figure('Name',sprintf('%s (%s to %s): Time Series',siteID,DMY(1),DMY(end)), 'Position', [50, 50, 1200, 600]);
    plot(DMY, velSigned_depthAvg, 'k-', 'LineWidth', 1);
    hold on;
    plot(DMY(floodSign_depthAvg), velSigned_depthAvg(floodSign_depthAvg), 'b.', 'MarkerSize', 10);
    plot(DMY(ebbSign_depthAvg), velSigned_depthAvg(ebbSign_depthAvg), 'r.', 'MarkerSize', 10);
    ylabel('Velocity (m/s)');
    xlabel('Time');
    datetick('x', 'mm/dd HH:MM', 'keepticks');
    title(sprintf('%s (%s to %s): Tidal Flow Analysis: Flood (+) and Ebb (-) Tides',siteID,DMY(1),DMY(end)));
    legend('Velocity', 'Flood Tide', 'Ebb Tide');
    grid on;   
    saveas(gcf, fullfile(mainOutputDir, '_TimeSeries.png'));
    
    %% Figure 3: Velocity Components
    
    figure('Name', sprintf('%s (%s to %s): Velocity Components',siteID,DMY(1),DMY(end)), 'Position', [50, 50, 1200, 400]);
    plot(DMY, Eas_depthAvg, 'r-', 'LineWidth', 1.2);
    hold on;
    plot(DMY, Nor_depthAvg, 'g-', 'LineWidth', 1.2);
    legend('Eastward', 'Northward');
    title(sprintf('%s (%s to %s): Velocity Components',siteID,DMY(1),DMY(end)));
    datetick('x', 'mm/dd HH:MM', 'keepticks');
    grid on;
    saveas(gcf, fullfile(mainOutputDir, '_Components.png'));
    
    %% Figure 4: Statistical Analysis
    floodStats = velSigned_depthAvg(floodSign_depthAvg);
    ebbStats = velSigned_depthAvg(ebbSign_depthAvg);

    figure('Name', sprintf('%s (%s to %s): Statistical Analysis',siteID,DMY(1),DMY(end)), 'Position', [50, 50, 800, 600]);
    
    subplot(2, 1, 1);
    histogram(floodStats, 20, 'FaceColor', 'b', 'FaceAlpha', 0.7);
    hold on;
    histogram(ebbStats, 20, 'FaceColor', 'r', 'FaceAlpha', 0.7);
    title(sprintf('%s (%s to %s): Velocity Distribution by Tidal Phase',siteID,DMY(1),DMY(end)));
    xlabel('Velocity (m/s)');
    ylabel('Frequency');
    legend('Flood Tide', 'Ebb Tide');
    grid on;

    subplot(2, 1, 2);
    boxplot([floodStats; -ebbStats], [ones(size(floodStats)); 2*ones(size(ebbStats))], ...
        'Labels', {'Flood', 'Ebb'}, 'Whisker', 1.5);
    title(sprintf('%s (%s to %s): Statistical Comparison of Flood and Ebb Magnitudes',siteID,DMY(1),DMY(end)));
    ylabel('Velocity Magnitude (m/s)');
    grid on;
    
    saveas(gcf,  fullfile(mainOutputDir, '_Statistics.png'));
    
    %% Figure 5: Current Rose
    figure('Name', sprintf('%s (%s to %s): Tidal Current Rose',siteID,DMY(1),DMY(end)), 'Position', [50, 50, 700, 600]);

    directions_rad = deg2rad(velDir_depthAvg);
    magnitudes = velSigned_depthAvg;

    polarhistogram(directions_rad, 36, 'DisplayStyle', 'stairs', 'Normalization', 'count');
    hold on;

    edges = linspace(0, 2*pi, 37);
    magnitude_edges = linspace(0, max(magnitudes)*1.1, 5);

    colormap(jet);
    colors = jet(length(magnitude_edges)-1);

    for i = length(magnitude_edges)-1:-1:1
        band_idx = magnitudes >= magnitude_edges(i) & magnitudes < magnitude_edges(i+1);
        if any(band_idx)
            h = polarhistogram(directions_rad(band_idx), edges, 'DisplayStyle', 'bar', 'Normalization', 'count');
            h.FaceColor = colors(i,:);
            h.FaceAlpha = 0.7;
        end
    end

    ax = gca;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    thetaticks(0:45:315);
    thetaticklabels({'N','NE','E','SE','S','SW','W','NW'});
    ax.FontName = 'Times New Roman';

    title(sprintf('%s (%s to %s): Current Rose (True North)',siteID,DMY(1),DMY(end)), 'FontName', 'Times New Roman');
    c = colorbar;
    c.Ticks = linspace(0, 1, length(magnitude_edges)-1);
    c.TickLabels = arrayfun(@(x,y) sprintf('%.1f-%.1f', x, y), ...
        magnitude_edges(1:end-1), magnitude_edges(2:end), 'UniformOutput', false);
    c.Label.String = 'Current Speed (m/s)';
    c.FontName = 'Times New Roman';
    
    saveas(gcf, fullfile(mainOutputDir, '_CurrentRose.png'));
    
    %% Figure 6: Tidal Cycles
    figure('Name', sprintf('%s (%s to %s): Tidal Cycles',siteID,DMY(1),DMY(end)), 'Position', [50, 50, 1200, 400]);

    subplot(2,1,1);
    plot(DMY, velSigned_depthAvg, 'k-', 'LineWidth', 0.5);
    hold on;
    plot(DMY(floodSign_depthAvg), velSigned_depthAvg(floodSign_depthAvg), 'b.', 'MarkerSize', 6);
    plot(DMY(ebbSign_depthAvg), velSigned_depthAvg(ebbSign_depthAvg), 'r.', 'MarkerSize', 6);

    yline(0, 'k--');

    ylabel('Signed Velocity (m/s)');
    title(sprintf('%s (%s to %s): Tidal Cycles',siteID,DMY(1),DMY(end)));
    grid on;
    datetick('x', 'mm/dd', 'keeplimits');

    subplot(2,1,2);
    edges = linspace(-max(abs(velSigned_depthAvg)), max(abs(velSigned_depthAvg)), 30);
    histogram(velSigned_depthAvg(floodSign_depthAvg), edges, 'FaceColor', 'b', 'FaceAlpha', 0.7);
    hold on;
    histogram(velSigned_depthAvg(ebbSign_depthAvg), edges, 'FaceColor', 'r', 'FaceAlpha', 0.7);
    xlabel('Signed Velocity (m/s)');
    ylabel('Frequency');
    legend('Flood', 'Ebb');
    grid on;
    
    saveas(gcf, fullfile(mainOutputDir, '_TidalCycles.png'));
    
end