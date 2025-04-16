%% ADCP Tidal Flow Analysis
% This script extracts depth-averaged velocity data from multiple ADCP files
% and identifies flood and ebb tides based on current direction.
%
% Author: Hossein Seyedzadeh
% Date: April 16, 2025
%
% The script performs the following operations:
% 1. Extracts depth-averaged velocity data from ADCP Excel files
% 2. Identifies flood and ebb tide periods based on flow direction
% 3. Calculates statistics for tidal flow patterns
% 4. Creates multiple visualizations of tidal patterns
% 5. Outputs results to both MAT file and text file
%
% Required functions:
% - ExtractCvelAvgWithTimestamps.m
% - identifyTidalFlow.m

%% Initialize workspace
clear;
clc;
close all;

% Define working folder
workfolder = 'C:\Users\Student\Documents\ADCP\Matlab\other-codes\depth-avg-analysis';

%% Step 1: Extract ADCP data from files
fprintf('Extracting data from ADCP files...\n');

% Extract velocity data and get file data in one call to avoid duplicate processing
[CvelAvg, fileData] = ExtractCvelAvgWithTimestamps(workfolder);

% Extract timestamps from the file data
if iscell(fileData)
    C = fileData{1};  % Start with first file
else
    C = fileData;
end

% Find timestamp column
DT_EST = contains(C.Properties.VariableNames, "EST");
if ~any(DT_EST)
    DT_EST = contains(C.Properties.VariableNames, "Date"); % Try alternate column name
end

% Extract timestamps
timeStamps = [];
for i = find(DT_EST)
    timeStamps = table2array(C(:,i));
    break; % Just use the first timestamp column found
end

% Ensure timestamps match data rows
if length(timeStamps) > size(CvelAvg, 1)
    timeStamps = timeStamps(1:size(CvelAvg, 1));
end

%% Step 2: Identify tidal flows using direction-based classification
fprintf('Identifying flood and ebb tides...\n');
% Extract velocity magnitude and direction
velocityMagnitude = CvelAvg(:,3);
direction = CvelAvg(:,4);

% Call identifyTidalFlow function to classify tidal phases
[floodVel, ebbVel, isFlood, isEbb, signedVelocity] = identifyTidalFlow(velocityMagnitude, direction);

%% Step 3: Calculate tidal flow statistics
fprintf('Computing tidal flow statistics...\n');

% Calculate average magnitudes for flood and ebb
floodMagnitude = mean(CvelAvg(isFlood,3), 'omitnan');
ebbMagnitude = mean(CvelAvg(isEbb,3), 'omitnan');

% Calculate dominant directions using circular statistics
floodDirections = CvelAvg(isFlood,4);
ebbDirections = CvelAvg(isEbb,4);

% Convert to radians for circular calculations
floodDirRadians = deg2rad(floodDirections);
ebbDirRadians = deg2rad(ebbDirections);

% Calculate mean direction using vector approach (circular mean)
floodMeanDir = mod(rad2deg(atan2(mean(sin(floodDirRadians), 'omitnan'), mean(cos(floodDirRadians), 'omitnan'))), 360);
ebbMeanDir = mod(rad2deg(atan2(mean(sin(ebbDirRadians), 'omitnan'), mean(cos(ebbDirRadians), 'omitnan'))), 360);

% Display results to console
fprintf('\n==== TIDAL FLOW ANALYSIS RESULTS ====\n');
fprintf('Time period: %s to %s\n', ...
    datestr(min(timeStamps), 'mm/dd/yyyy HH:MM'), ...
    datestr(max(timeStamps), 'mm/dd/yyyy HH:MM'));
fprintf('Total data points: %d\n', length(timeStamps));

fprintf('\nFLOOD TIDE:\n');
fprintf('  Average Magnitude: %.2f m/s\n', floodMagnitude);
fprintf('  Dominant Direction: %.1f° (True North)\n', floodMeanDir);
fprintf('  Duration: %.1f%% of record\n', 100*sum(isFlood)/length(isFlood));

fprintf('\nEBB TIDE:\n');
fprintf('  Average Magnitude: %.2f m/s\n', ebbMagnitude);
fprintf('  Dominant Direction: %.1f° (True North)\n', ebbMeanDir);
fprintf('  Duration: %.1f%% of record\n', 100*sum(isEbb)/length(isEbb));

%% Step 4: Create visualizations

%% Figure 1: Tidal Flow Summary (Bar chart and direction rose)
figure('Name', 'Tidal Flow Summary', 'Position', [50, 50, 1000, 450]);

% Create first axes for bar chart with explicit position
ax1 = axes('Position', [0.1 0.15 0.35 0.7]);
bar([floodMagnitude, ebbMagnitude]);
set(ax1, 'XTickLabel', {'Flood', 'Ebb'});
ylabel('Average Magnitude (m/s)');
title('Tidal Current Magnitudes');
grid on;

% Create second axes for polar plot with explicit position
ax2 = axes('Position', [0.55 0.1 0.4 0.8]);
delete(ax2);  % Delete the regular axes
ax2 = polaraxes('Position', [0.55 0.1 0.4 0.8]);  % Create polar axes directly
hold on;

% Set correct orientation for true north
ax2.ThetaZeroLocation = 'top';    % 0 degrees at top (North)
ax2.ThetaDir = 'clockwise';       % Clockwise angle increase (like a compass)

% Calculate radians for polar plotting
floodRad = deg2rad(floodMeanDir);
ebbRad = deg2rad(ebbMeanDir);

% Plot the vectors
polarplot([0 floodRad], [0 floodMagnitude], 'LineWidth', 2, 'Color', [0.2 0.4 0.9]);
polarplot([0 ebbRad], [0 ebbMagnitude], 'LineWidth', 2, 'Color', [0.9 0.2 0.2]);

% Add markers at the tips
polarscatter(floodRad, floodMagnitude, 70, [0.2 0.4 0.9], 'filled');
polarscatter(ebbRad, ebbMagnitude, 70, [0.9 0.2 0.2], 'filled');

% Configure compass ticks and labels
thetaticks(0:45:315);  
thetaticklabels({'N','NE','E','SE','S','SW','W','NW'});

% Add text labels showing exact directions
text(floodRad, floodMagnitude*1.15, sprintf('%.1f° T', floodMeanDir), 'Color', [0.2 0.4 0.9], 'FontWeight', 'bold');
text(ebbRad, ebbMagnitude*1.15, sprintf('%.1f° T', ebbMeanDir), 'Color', [0.9 0.2 0.2], 'FontWeight', 'bold');

% Add legend
legend({'Flood', 'Ebb', '', ''}, 'Location', 'southoutside');
title('Dominant Current Directions (True North)');

%% Figure 2: Time Series of Signed Velocities
figure('Name', 'Tidal Flow Analysis', 'Position', [50, 50, 1200, 600]);
% Plot the full time series
plot(timeStamps, signedVelocity, 'k-', 'LineWidth', 1);
hold on;
% Highlight flood and ebb periods
plot(timeStamps(isFlood), signedVelocity(isFlood), 'b.', 'MarkerSize', 10);
plot(timeStamps(isEbb), signedVelocity(isEbb), 'r.', 'MarkerSize', 10);
ylabel('Velocity (m/s)');
xlabel('Time');
datetick('x', 'mm/dd HH:MM', 'keepticks');
title('Tidal Flow Analysis: Flood (+) and Ebb (-) Tides');
legend('Velocity', 'Flood Tide', 'Ebb Tide');
grid on;

%% Figure 3: Velocity Components (East and North)
figure('Name', 'Velocity Components', 'Position', [50, 50, 1200, 400]);
plot(timeStamps, CvelAvg(:,1), 'r-', 'LineWidth', 1.2); % East
hold on;
plot(timeStamps, CvelAvg(:,2), 'g-', 'LineWidth', 1.2); % North
legend('Eastward', 'Northward');
title('Velocity Components');
datetick('x', 'mm/dd HH:MM', 'keepticks');
grid on;

%% Figure 4: Statistical Analysis of Flood and Ebb Tides
% Extract flood and ebb statistics
floodStats = signedVelocity(isFlood);
ebbStats = signedVelocity(isEbb);

figure('Name', 'Tidal Statistics', 'Position', [50, 50, 800, 600]);
% Histogram of velocity distribution
subplot(2, 1, 1);
histogram(floodStats, 20, 'FaceColor', 'b', 'FaceAlpha', 0.7);
hold on;
histogram(ebbStats, 20, 'FaceColor', 'r', 'FaceAlpha', 0.7);
title('Velocity Distribution by Tidal Phase');
xlabel('Velocity (m/s)');
ylabel('Frequency');
legend('Flood Tide', 'Ebb Tide');
grid on;

% Boxplot for comparing magnitude distributions
subplot(2, 1, 2);
boxplot([floodStats; -ebbStats], [ones(size(floodStats)); 2*ones(size(ebbStats))], ...
    'Labels', {'Flood', 'Ebb'}, 'Whisker', 1.5);
title('Statistical Comparison of Flood and Ebb Magnitudes');
ylabel('Velocity Magnitude (m/s)');
grid on;

%% Figure 5: Current Rose Diagram
figure('Name', 'Current Rose', 'Position', [50, 50, 700, 600]);

% Convert directions to radians for polar plotting
directions_rad = deg2rad(CvelAvg(:,4));
magnitudes = CvelAvg(:,3);

% Create polar histogram outline
polarhistogram(directions_rad, 36, 'DisplayStyle', 'stairs', 'Normalization', 'count');
hold on;

% Create a detailed current rose with intensity bands
edges = linspace(0, 2*pi, 37); % 36 bins (every 10 degrees)
magnitude_edges = linspace(0, max(magnitudes)*1.1, 5); % 4 magnitude bands

% Loop through magnitude bands from highest to lowest
colormap(jet);
colors = jet(length(magnitude_edges)-1);

for i = length(magnitude_edges)-1:-1:1
    % Select data in this magnitude band
    band_idx = magnitudes >= magnitude_edges(i) & magnitudes < magnitude_edges(i+1);
    if any(band_idx)
        % Create polar histogram for this band
        h = polarhistogram(directions_rad(band_idx), edges, 'DisplayStyle', 'bar', 'Normalization', 'count');
        h.FaceColor = colors(i,:);
        h.FaceAlpha = 0.7;
    end
end

% Configure polar axis to show compass directions correctly
ax = gca;
ax.ThetaZeroLocation = 'top';    % Put 0 degrees at the top (North)
ax.ThetaDir = 'clockwise';       % Make angles increase clockwise (like a compass)
thetaticks(0:45:315);            % Show tick marks every 45 degrees
thetaticklabels({'N','NE','E','SE','S','SW','W','NW'}); % Label with compass directions

% Add title and colorbar
title('Current Rose Diagram (True North)');
c = colorbar;
c.Ticks = linspace(0, 1, length(magnitude_edges)-1);
c.TickLabels = arrayfun(@(x,y) sprintf('%.1f-%.1f', x, y), ...
    magnitude_edges(1:end-1), magnitude_edges(2:end), 'UniformOutput', false);
c.Label.String = 'Current Speed (m/s)';

fprintf('Current rose diagram created with correct True North orientation.\n');

%% Step 5: Save results to files

% Save MATLAB data for future reference
fprintf('Saving results to MAT file...\n');
matFileName = fullfile(workfolder, 'TidalAnalysisResults.mat');
save(matFileName, 'CvelAvg', 'floodVel', 'ebbVel', 'isFlood', 'isEbb', ...
     'signedVelocity', 'floodMagnitude', 'ebbMagnitude', 'floodMeanDir', 'ebbMeanDir', 'timeStamps');
fprintf('Results saved to %s\n', matFileName);

%% Write results to text output file
fprintf('Writing results to output file...\n');
outputFilename = fullfile(workfolder, 'TidalFlowAnalysis_Results.txt');
fid = fopen(outputFilename, 'w');

% Write header information
fprintf(fid, '===============================================\n');
fprintf(fid, '          ADCP TIDAL FLOW ANALYSIS RESULTS     \n');
fprintf(fid, '===============================================\n\n');
fprintf(fid, 'Analysis Date: %s\n\n', datestr(now, 'mmm dd, yyyy HH:MM:SS'));
fprintf(fid, 'Data Source: %s\n', workfolder);
fprintf(fid, 'Time period: %s to %s\n', ...
    datestr(min(timeStamps), 'mm/dd/yyyy HH:MM'), ...
    datestr(max(timeStamps), 'mm/dd/yyyy HH:MM'));
fprintf(fid, 'Total data points: %d\n\n', length(timeStamps));

% Write flood tide statistics
fprintf(fid, 'FLOOD TIDE STATISTICS:\n');
fprintf(fid, '-------------------------\n');
fprintf(fid, '  Average Magnitude: %.2f m/s\n', floodMagnitude);
fprintf(fid, '  Dominant Direction: %.1f° (True North)\n', floodMeanDir);
fprintf(fid, '  Duration: %.1f%% of record\n', 100*sum(isFlood)/length(isFlood));
fprintf(fid, '  Total Samples: %d\n\n', sum(isFlood));

% Write ebb tide statistics
fprintf(fid, 'EBB TIDE STATISTICS:\n');
fprintf(fid, '-------------------------\n');
fprintf(fid, '  Average Magnitude: %.2f m/s\n', ebbMagnitude);
fprintf(fid, '  Dominant Direction: %.1f° (True North)\n', ebbMeanDir);
fprintf(fid, '  Duration: %.1f%% of record\n', 100*sum(isEbb)/length(isEbb));
fprintf(fid, '  Total Samples: %d\n\n', sum(isEbb));

% Write additional statistics
fprintf(fid, 'ADDITIONAL STATISTICS:\n');
fprintf(fid, '-------------------------\n');
fprintf(fid, '  Max Flood Velocity: %.2f m/s\n', max(floodVel));
fprintf(fid, '  Max Ebb Velocity: %.2f m/s\n', max(ebbVel));
fprintf(fid, '  Flood/Ebb Magnitude Ratio: %.2f\n\n', floodMagnitude/ebbMagnitude);

% Write information about saved files
fprintf(fid, 'SAVED OUTPUT FILES:\n');
fprintf(fid, '-------------------------\n');
fprintf(fid, '  MATLAB Data: %s\n', matFileName);
fprintf(fid, '  Results Summary: %s\n', outputFilename);

% Close the file
fclose(fid);
fprintf('Results written to %s\n', outputFilename);
fprintf('Analysis complete!\n');