function TimeAveragedProfiles(workfolder)
    %% TIMEAVERAGEDPROFILES - Generates time-averaged velocity profiles from ADCP data
    %
    % Combines data from multiple ADCP deployment periods into a continuous timeline 
    % and calculates a time-averaged velocity profile.
    %
    % Creator: Hossein Seyedzadeh
    % Date: April 10, 2025
    % Last Modified: April 12, 2025
    %
    % Syntax:
    %   TimeAveragedProfiles(workfolder)
    %
    % Inputs:
    %   workfolder - Path to folder containing ADCP Excel files
    %
    % Outputs:
    %   None - Creates a figure showing the time-averaged velocity profile
    %
    % Example:
    %   TimeAveragedProfiles(pwd)
    %
    % Notes:
    %   - Processes Excel files from workfolder containing ADCP data
    %   - Can handle multiple deployments with different bin structures
    %   - Handles multiple sites (e.g., Station A, Station B)
    %   - Removes problematic bottom measurements
    %   - Calculates various statistics (mean velocity, duration, etc.)
    
    %% STEP 1: INITIALIZATION
    
    % List of sheet names to try for each deployment period
    % Site1 corresponds to Station B, Site2 corresponds to Station A
    sheetPairs = {{"Station B", "Station B-2"}, {"Station A", "Station A-2"}};
    
    % Get list of Excel files, excluding temporary Excel files (starting with ~$)
    workfolder = pwd;
    files = dir(fullfile(workfolder, '*.xlsx'));
    files = files(~startsWith({files.name}, '~$'));
    
    % Initialize data structure to store data by site
    allDataBySite = struct();
    
    %% STEP 2: DATA LOADING AND COMBINATION
    % Process each Excel file
    for fileIdx = 1:length(files)
        [~, fileName] = fileparts(files(fileIdx).name);
        fprintf('\n==== Processing file %d/%d: %s ====\n', fileIdx, length(files), files(fileIdx).name);
        
        % Try each possible site (Station A, Station B)
        for siteIdx = 1:length(sheetPairs)
            siteName = sprintf('Site%d', siteIdx); % Site1 = B, Site2 = A
            
            % Try each sheet name for this site
            for sheetIdx = 1:length(sheetPairs{siteIdx})
                currentSheet = sheetPairs{siteIdx}{sheetIdx};
                try
                    fprintf('  Trying sheet: %s\n', currentSheet);
                    
                    % Extract deployment data from this Excel file and sheet
                    [DMY, magData, dirData, depths, binDataNumeric, magBins, fileStats] = extractDeploymentData(workfolder, files(fileIdx).name, currentSheet);
                    
                    % If successful, report the number of time points read
                    fprintf('  Success! Read %d time points from sheet: %s\n', length(DMY), currentSheet);
                    
                    % Initialize site data structure if this is the first file for this site
                    if ~isfield(allDataBySite, siteName)
                        allDataBySite.(siteName) = struct('DMY', [], 'magData', [], 'depths', [], 'filesProcessed', 0);
                    end
                    
                    % Get current site data for updating
                    siteData = allDataBySite.(siteName);
                    
                    %% STEP 3: DATA MERGING LOGIC
                    % Handle first data for this site
                    if isempty(siteData.depths)
                        % First data for this site - use as is
                        siteData.depths = depths;
                        siteData.DMY = DMY;
                        siteData.magData = magData;
                        siteData.dirData = dirData;
                    else
                        % We already have data for this site - need to merge
                        % Check if depths match between files (within small tolerance)
                        if length(depths) == length(siteData.depths) && all(abs(depths - siteData.depths) < 0.01)
                            % Same bin structure, just append time series data
                            siteData.DMY = [siteData.DMY; DMY];
                            siteData.magData = [siteData.magData; magData];
                            siteData.dirData = [siteData.dirData; dirData];
                        else
                            % Bins don't match - need to interpolate new data to match existing depth grid
                            fprintf('  Warning: Depth bins in this file don''t match previous file, interpolating...\n');
                            
                            % Interpolate each profile to match site depths
                            magDataInterp = zeros(length(DMY), length(siteData.depths));
                            dirDataInterp = zeros(length(DMY), length(siteData.depths));
                            for i = 1:length(DMY)
                                % Linear interpolation with extrapolation for out-of-range values
                                magDataInterp(i,:) = interp1(depths, magData(i,:), siteData.depths, 'linear', 'extrap');
                                dirDataInterp(i,:) = interp1(depths, dirData(i,:), siteData.depths);
                            end
                            
                            % Append interpolated data to existing time series
                            siteData.DMY = [siteData.DMY; DMY];
                            siteData.magData = [siteData.magData; magDataInterp];
                            siteData.dirData = [siteData.dirData; dirDataInterp];
                        end
                    end
                    
                    % Update processed file count for this site
                    siteData.filesProcessed = siteData.filesProcessed + 1;
                    
                    % Store updated data back in the structure
                    allDataBySite.(siteName) = siteData;
                    
                    % Break the sheet loop - we found a working sheet for this file
                    break;
                catch e
                    % Report failure with specific error message
                    fprintf('  Failed with sheet %s: %s\n', currentSheet, e.message);
                end
            end
            
            % If we processed this file with a sheet from this site, don't try other sites
            % This prevents processing the same file twice if it contains data for multiple sites
            if isfield(allDataBySite, siteName) && allDataBySite.(siteName).filesProcessed > 0
                if allDataBySite.(siteName).filesProcessed == fileIdx
                    break;
                end
            end
        end
    end
    
    %% STEP 4: SELECT BEST SITE FOR ANALYSIS
    % Process the site with the most complete data
    siteNames = fieldnames(allDataBySite);
    if isempty(siteNames)
        error('No data could be processed from any file');
    end
    
    % Find the site with the most files processed
    maxFiles = 0;
    bestSite = '';
    for i = 1:length(siteNames)
        if allDataBySite.(siteNames{i}).filesProcessed > maxFiles
            maxFiles = allDataBySite.(siteNames{i}).filesProcessed;
            bestSite = siteNames{i};
        end
    end
    
    fprintf('\nUsing data from %s with %d files processed\n', bestSite, maxFiles);
    
    % Get data from the best site
    siteData = allDataBySite.(bestSite);
    allDMY = siteData.DMY;
    allMagData = siteData.magData;
    allDirData = siteData.dirData;
    combinedDepths = siteData.depths;
    
    %% STEP 4B: SEPARATION INTO EBB AND FLOOD
    % Binary conditional statement with loop assigns signs to denote
    % west as positive (0 < dir < 180) and east as negative (0 < dir  <
    % 270)
    
%     tsAvg = mean(diff(allDMY)); % For given break in datasets' timeframes, the code takes mean
%     tsTide = hours(12); % Tide is semidiurnal period of 12 hours and 25 minutes; however, slightly shorter period accommodates more peaks
%     tsMinPeak = 1;
%     while tsMinPeak*tsAvg <= tsTide
%         tsMinPeak = tsMinPeak+1;
%     end
    
    [allDMY, sortIdx] = sort(allDMY);
    allDirData = allDirData(sortIdx, :);
    avgDir = mean(allDirData, 2, 'omitnan');
    binaryDir = zeros(size(avgDir));
    for i = 1:length(allDirData)
        if (0 <= avgDir(i)) && (avgDir(i) < 180)
            binaryDir(i) = -1;
        elseif (180 <= avgDir(i)) && (avgDir(i) < 360)
            binaryDir(i) = 1;
        elseif isnan(allDirData(i))
            warning("There is a NaN value at index %i", i);
        else
            warning("There is no value at index %i", i);
        end
    end
    
    allMagData = allMagData(sortIdx, :);
    allVelData = allMagData.*binaryDir;
%     velAvgMin = min(allVelData,[],2);
%     velAvgMax = max(allVelData,[],2);
    velAvgMean = mean(allVelData,2,"omitnan");
    velAvgstd = std(velAvgMean,0,"omitnan");
    
    %[tides,tideLocs,tideTol] = findpeaks(velMagAvg,'MinPeakProminence',velAvgstd,'MinPeakDistance',timestepHalfTide);
    [floods,floodLocs] = findpeaks(velAvgMean,'MinPeakProminence',velAvgstd);
    [ebbs,ebbLocs] = findpeaks(-velAvgMean,'MinPeakProminence',velAvgstd);
    
    floodVel = zeros(size(floods));
    for i = 1:length(allMagData)
        for j = 1:length(floodLocs)
            if i == floodLocs(j)
                for k = 1:size(allMagData,2)
                    floodVel(j,k) = allMagData(i,k);
                end
            else
                continue
            end
        end
    end
    
    ebbVel = zeros(size(ebbs));
    for i = 1:length(allMagData)
        for j = 1:length(ebbLocs)
            if i == ebbLocs(j)
                for k = 1:size(allMagData,2)
                    ebbVel(j,k) = allMagData(i,k);
                end
            else
                continue
            end
        end
    end

    
    %% STEP 5: TIME-AVERAGING AND PROFILE CREATION
    % Sort the combined data chronologically (only if we have data)
    
    if ~isempty(allDMY) && ~isempty(allMagData)
        % Sort data chronologically to ensure proper time ordering
        [allDMY, sortIdx] = sort(allDMY);
        allMagData = allMagData(sortIdx, :);
        
        % THE TIME AVERAGING STEP
        % Calculate time-averaged velocity profile from the combined data
        % Average along dimension 1 (rows/time) and handle NaN values
        avgMag = mean(allMagData, 1, 'omitnan');
        avgFlood = mean(floodVel, 1, 'omitnan');
        avgEbb = mean(ebbVel, 1, 'omitnan');
        %avgDir = mean(allDirData, 1, 'omitnan');
        
        % Convert from mm/s to m/s (if not already done)
        avgMag = avgMag / 1000;
        avgFlood = avgFlood / 1000;
        avgEbb = avgEbb / 1000;

        %% STEP 6: SURFACE DATA CORRECTION
        % Handle the shallowest 4 points by the surface, which are often
        % contaminated with NaN % WE CAN MODIFY THIS SUBSECTION TO IDENTIFY
        % THE CONTAMINATED DATA ACCORDING TO THE PRESENCE OF NAN VALUES
        numBottom = 5;  % Number of bottom points to fix
        if length(avgMag) > numBottom + 2
            fprintf('Fixing the %d shallowest velocity measurements (near seafloor)\n', numBottom);
            
            % Option 1: Exclude bottom points completely
            avgMag = avgMag(1:end-numBottom);
            avgFlood = avgFlood(1:end-numBottom);
            avgEbb = avgEbb(1:end-numBottom);
            combinedDepths = combinedDepths(1:end-numBottom);
            
            fprintf('Bottom points removed from profile\n');
         
        end
        
        %% STEP 7: VISUALIZATION AND REPORTING
        % Plot the time-averaged magnitude profile
        
        % Create a figure for the combined analysis, with formatted dates
        startDate = min(allDMY);
        endDate = max(allDMY);
        figure('Name', sprintf('Time-Averaged Velocity Profiles\nCombined Deployment Periods: %s to %s', ...
            datestr(startDate, 'mm/dd/yyyy'), datestr(endDate, 'mm/dd/yyyy'), 'Position', [50, 50, 900, 700]));
        subplot(2,2,1);
        plot(avgMag, combinedDepths, 'b-o', 'LineWidth', 2);
        
        % Add data point labels at regular intervals for readability
        labelStep = max(1, round(length(combinedDepths)/10));
        for i = 1:labelStep:length(combinedDepths)
            text(avgMag(i) + 0.01, combinedDepths(i), sprintf('%.3f', avgMag(i)), 'FontSize', 8);
        end
        
        % Set plot properties
        grid on;
        xlabel('Velocity Magnitude (m/s)', 'FontSize', 12);
        ylabel('Distance from bed = 0 (m)', 'FontSize', 12);

        % Format date range for title
        title(sprintf('Time-Averaged Magnitude Profile\nCombined Deployment Periods: %s to %s', ...
            datestr(startDate, 'mm/dd/yyyy'), datestr(endDate, 'mm/dd/yyyy')), 'FontSize', 10);
        
        % Plot the time-averaged flood profile
        subplot(2,2,2);
        plot(avgFlood, combinedDepths, 'r-o', 'LineWidth', 2);
        
        % Add data point labels at regular intervals for readability
        labelStep = max(1, round(length(combinedDepths)/10));
        for i = 1:labelStep:length(combinedDepths)
            text(avgFlood(i) + 0.01, combinedDepths(i), sprintf('%.3f', avgFlood(i)), 'FontSize', 8);
        end
        
        % Set plot properties
        grid on;
        xlabel('Flood Magnitude (m/s)', 'FontSize', 12);
        ylabel('Distance from bed = 0 (m)', 'FontSize', 12);

        % Format date range for title
        
        title(sprintf('Time-Averaged Flood Profile\nCombined Deployment Periods: %s to %s', ...
            datestr(startDate, 'mm/dd/yyyy'), datestr(endDate, 'mm/dd/yyyy')), 'FontSize', 10);
        
                % Plot the time-averaged magnitude profile
        subplot(2,2,3);
        plot(avgEbb, combinedDepths, 'g-o', 'LineWidth', 2);
        
        % Add data point labels at regular intervals for readability
        labelStep = max(1, round(length(combinedDepths)/10));
        for i = 1:labelStep:length(combinedDepths)
            text(avgEbb(i) + 0.01, combinedDepths(i), sprintf('%.3f', avgEbb(i)), 'FontSize', 8);
        end
        
        % Set plot properties
        grid on;
        xlabel('Ebb Magnitude (m/s)', 'FontSize', 12);
        ylabel('Distance from bed = 0 (m)', 'FontSize', 12);

        % Format date range for title
        startDate = min(allDMY);
        endDate = max(allDMY);
        title(sprintf('Time-Averaged Ebb Profile\nCombined Deployment Periods: %s to %s', ...
            datestr(startDate, 'mm/dd/yyyy'), datestr(endDate, 'mm/dd/yyyy')), 'FontSize', 10);
        
        % Display comprehensive statistics about the combined dataset
        fprintf('\n==== COMBINED DEPLOYMENT ANALYSIS ====\n');
        fprintf('Date range: %s to %s\n', ...
            datestr(startDate, 'mm/dd/yyyy HH:MM'), ...
            datestr(endDate, 'mm/dd/yyyy HH:MM'));
        fprintf('  Number of profiles averaged: %d\n', size(allMagData, 1));
        fprintf('  Overall average magnitude: %.3f m/s\n', mean(avgMag, 'omitnan'));
        fprintf('  Near-Surface magnitude: %.3f m/s\n', avgMag(end));
        fprintf('  Near-Bottom magnitude: %.3f m/s\n', avgMag(1));
        fprintf('  Number of ebbs: %d\n', size(ebbs, 1));
        fprintf('  Overall average ebb: %.3f m/s\n', mean(avgEbb, 'omitnan'));
        fprintf('  Near-Surface ebb: %.3f m/s\n', avgEbb(end));
        fprintf('  Near-Bottom ebb: %.3f m/s\n', avgEbb(1));
        fprintf('  Number of floods: %d\n', size(floods, 1));
        fprintf('  Overall average flood: %.3f m/s\n', mean(avgEbb, 'omitnan'));
        fprintf('  Near-Surface flood: %.3f m/s\n', avgEbb(end));
        fprintf('  Near-Bottom flood: %.3f m/s\n', avgEbb(1));

        % After the visualization section (around line 350), add:

%% STEP 8: FIT PROFILES TO ANALYTICAL MODELS
% Fit profiles to log law and power law for boundary condition formulation
fprintf('\nFitting velocity profiles to analytical models...\n');

% Create a new figure for model fits
figure('Name', sprintf('Velocity Profile Models\nCombined Deployment Periods: %s to %s', ...
    datestr(startDate, 'mm/dd/yyyy'), datestr(endDate, 'mm/dd/yyyy')), 'Position', [100, 100, 1000, 800]);

% 1. Polynomial fit
subplot(2,3,1);
[ebbPolyCoeffs, ebbPolyFormula] = createPolynomialFit(avgEbb, combinedDepths, 'Ebb', 5);
title('Polynomial Fit - Ebb');

subplot(2,3,2);
[floodPolyCoeffs, floodPolyFormula] = createPolynomialFit(avgFlood, combinedDepths, 'Flood', 5);
title('Polynomial Fit - Flood');

% 2. Power law fit
subplot(2,3,3);
[ebbPowerFormula, ebbAlpha, ebbSurface] = createPowerLawProfile(avgEbb, combinedDepths, 'Ebb');
title('Power Law Fit - Ebb');

subplot(2,3,4);
[floodPowerFormula, floodAlpha, floodSurface] = createPowerLawProfile(avgFlood, combinedDepths, 'Flood');
title('Power Law Fit - Flood');

% 3. Log law fit
subplot(2,3,5);
[ebbLogFormula, ebbUstar, ebbZ0] = createLogLawProfile(avgEbb, combinedDepths, 'Ebb');
title('Log Law Fit - Ebb');

subplot(2,3,6);
[floodLogFormula, floodUstar, floodZ0] = createLogLawProfile(avgFlood, combinedDepths, 'Flood');
title('Log Law Fit - Flood');

% Print formulas for boundary conditions
fprintf('\n==== BOUNDARY CONDITION FORMULAS ====\n');
fprintf('EBB VELOCITY PROFILE FORMULAS:\n');
fprintf('1. Polynomial (degree 5): v(z) = %.6f*z^5 + %.6f*z^4 + %.6f*z^3 + %.6f*z^2 + %.6f*z + %.6f\n', ...
    ebbPolyCoeffs(1), ebbPolyCoeffs(2), ebbPolyCoeffs(3), ebbPolyCoeffs(4), ebbPolyCoeffs(5), ebbPolyCoeffs(6));
fprintf('2. Power law: v(z) = %.4f * (z/%.2f)^%.4f\n', ebbSurface, combinedDepths(end), ebbAlpha);
fprintf('3. Log law: v(z) = %.4f/0.41 * ln(z/%.6f)\n', ebbUstar, ebbZ0);
fprintf('   where u_* = %.4f m/s, z_0 = %.6f m\n', ebbUstar, ebbZ0);

fprintf('\nFLOOD VELOCITY PROFILE FORMULAS:\n');
fprintf('1. Polynomial (degree 5): v(z) = %.6f*z^5 + %.6f*z^4 + %.6f*z^3 + %.6f*z^2 + %.6f*z + %.6f\n', ...
    floodPolyCoeffs(1), floodPolyCoeffs(2), floodPolyCoeffs(3), floodPolyCoeffs(4), floodPolyCoeffs(5), floodPolyCoeffs(6));
fprintf('2. Power law: v(z) = %.4f * (z/%.2f)^%.4f\n', floodSurface, combinedDepths(end), floodAlpha);
fprintf('3. Log law: v(z) = %.4f/0.41 * ln(z/%.6f)\n', floodUstar, floodZ0);
fprintf('   where u_* = %.4f m/s, z_0 = %.6f m\n', floodUstar, floodZ0);

% Save formulas to a text file
formulaFile = fullfile(workfolder, 'boundary_condition_formulas.txt');
fid = fopen(formulaFile, 'w');
fprintf(fid, 'BOUNDARY CONDITION FORMULAS\n');
fprintf(fid, 'Generated from ADCP data: %s to %s\n\n', ...
    datestr(startDate, 'mm/dd/yyyy'), datestr(endDate, 'mm/dd/yyyy'));
fprintf(fid, 'EBB VELOCITY PROFILE FORMULAS:\n');
fprintf(fid, '1. Polynomial: v(z) = %.6f*z^5 + %.6f*z^4 + %.6f*z^3 + %.6f*z^2 + %.6f*z + %.6f\n', ...
    ebbPolyCoeffs(1), ebbPolyCoeffs(2), ebbPolyCoeffs(3), ebbPolyCoeffs(4), ebbPolyCoeffs(5), ebbPolyCoeffs(6));
fprintf(fid, '2. Power law: v(z) = %.4f * (z/%.2f)^%.4f\n', ebbSurface, combinedDepths(end), ebbAlpha);
fprintf(fid, '3. Log law: v(z) = %.4f/0.41 * ln(z/%.6f)\n', ebbUstar, ebbZ0);
fprintf(fid, 'FLOOD VELOCITY PROFILE FORMULAS:\n');
fprintf(fid, '1. Polynomial: v(z) = %.6f*z^5 + %.6f*z^4 + %.6f*z^3 + %.6f*z^2 + %.6f*z + %.6f\n', ...
    floodPolyCoeffs(1), floodPolyCoeffs(2), floodPolyCoeffs(3), floodPolyCoeffs(4), floodPolyCoeffs(5), floodPolyCoeffs(6));
fprintf(fid, '2. Power law: v(z) = %.4f * (z/%.2f)^%.4f\n', floodSurface, combinedDepths(end), floodAlpha);
fprintf(fid, '3. Log law: v(z) = %.4f/0.41 * ln(z/%.6f)\n', floodUstar, floodZ0);
fclose(fid);
fprintf('Formulas saved to: %s\n', formulaFile);
        
        % Calculate temporal statistics
        totalDuration = endDate - startDate;
        daysElapsed = days(totalDuration);
        fprintf('  Total monitoring period: %.1f days\n', daysElapsed);
        
        % Add annotation with key information
        annotation('textbox', [0.15, 0.01, 0.7, 0.05], 'String', ...
            sprintf('Combined analysis of %d profiles across %.1f days', size(allMagData, 1), daysElapsed), ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    else
        title('No valid data could be processed', 'FontSize', 14);
        fprintf('ERROR: No valid data could be processed from any file\n');
    end
end

%% HELPER FUNCTIONS
function [DMY, magData, dirData, depths, binDataNumeric, magBins, stats] = extractDeploymentData(workfolder, fileName, sheetName)
    %% EXTRACTDEPLOYMENTDATA - Extracts ADCP data from a specific deployment file
    %
    % This function reads and processes ADCP data from a specific Excel file and sheet.
    % It extracts timestamps, velocity magnitudes, depths, and other relevant information.
    %
    % Inputs:
    %   workfolder - Path to the folder containing the Excel file
    %   fileName   - Name of the Excel file to process
    %   sheetName  - Name of the sheet to read from the Excel file
    %
    % Outputs:
    %   DMY            - Array of datetime values (timestamps)
    %   magData        - Matrix of velocity magnitude data (rows=time, cols=depths)
    %   depths         - Array of depth values for each measurement bin
    %   binDataNumeric - Full matrix of numeric bin data including all measurements
    %   magBins        - Logical array identifying magnitude columns
    %   stats          - Structure containing statistics about this deployment
    
    % Build filepath for this specific Excel file
    filepath = fullfile(workfolder, fileName);
    
    % Define keywords for data extraction
    keyword = "nor";     % For velocity component identification
    keyword2 = "distance"; % For depth identification
    keywordbin = "bin";    % For bin number identification
    
    % Read the Excel file using specialized helper function
    [C, DistT] = readSpecificFile(filepath, keyword, keyword2, keywordbin, sheetName);
    
    %% Extract timestamps
    % Find date/time columns using header names
    DT = contains(C.Properties.VariableNames, "Date");
    DT_EST = contains(C.Properties.VariableNames, "EST");
    dateCols = DT | DT_EST;
    lastDateCol = find(dateCols, 1, 'last');
    
    % Extract datetime values
    for i = find(DT_EST)
        DMY = table2array(C(:,i));
    end
    
    % Remove missing timestamp values
    if any(ismissing(DMY),'all')
        DMY = rmmissing(DMY,1);
    end
    
    %% Extract velocity data
    % Get velocity data columns (after date columns)
    Cmod = C(:, lastDateCol+1:end);
    
    % Remove empty rows and columns
    absentCols = all(ismissing(Cmod), 1);
    absentRows = all(ismissing(Cmod), 2);
    Cmod(:, absentCols) = [];
    Cmod(absentRows, :) = [];
    
    % Trim timestamps to match data rows if needed
    if length(DMY) > size(Cmod, 1)
        DMY = DMY(1:size(Cmod, 1));
    end
    
    % Analyze variable names to identify types of measurements
    CmodVarNames = Cmod.Properties.VariableNames;
    velVars = contains(CmodVarNames, ["Eas", "Nor", "Mag", "Dir"]); % Velocity components
    velNames = CmodVarNames(velVars);
    
    % Count unique measurement types (e.g., East, North, Magnitude, Direction)
    velNums = unique(extractBefore(velNames, "_")); 
    velGages = numel(velNums);
    
    % Separate depth-averaged and bin-specific data
    % First velGages columns are typically depth-averaged values
    velAvg = Cmod(:, 1:velGages);
    binData = Cmod(:, velGages+1:end);  % Remaining are bin-specific data
    
    % Get bin data variable names
    binVarNames = binData.Properties.VariableNames;
    
    %% Extract bin numbers from variable names
    binNumbers = zeros(1, length(binVarNames));
    for i = 1:length(binVarNames)
        varName = binVarNames{i};
        % Try to extract bin number from end of variable name
        matches = regexp(varName, '_(\d+)$', 'tokens');
        if ~isempty(matches) && ~isempty(matches{1})
            binNumbers(i) = str2double(matches{1}{1});
        end
    end
    
    % Identify specific column types
    magBins = contains(binVarNames, "Mag"); % Magnitude columns
    dirBins = contains(binVarNames, "Dir"); % Direction columns
    
    % Extract unique bin numbers for magnitude columns
    magBinNumbers = binNumbers(magBins);
    uniqueBins = unique(magBinNumbers);
    numUniqueBins = length(uniqueBins);
    
    fprintf('Found %d unique depth bins in file %s\n', numUniqueBins, fileName);
    
    %% Extract depths from the distance table
    DistBinVars = string(DistT.Properties.RowNames);
    distIdx = find(any(strncmpi(DistBinVars, "distance", 8), 2), 1);
    DistT_ind = DistT(distIdx,:);
    DistA = str2double(table2array(DistT_ind));
    
    % Get valid depths (remove NaNs)
    allDepths = DistA(~isnan(DistA));
    
    %% Match depths to bin numbers
    % ADCP data often has depths repeated for each measurement type
    if length(allDepths) >= numUniqueBins*4
        % If depths are repeated for each measurement type (E,N,M,D), take every 4th
        depthStep = 4;  
        depths = allDepths(1:depthStep:end);
        
        if length(depths) < numUniqueBins
            % If we don't have enough depths, revert to using all available depths
            depths = allDepths(1:numUniqueBins);
            fprintf('Warning: Depth count doesn''t match unique bins. Using first %d depths.\n', numUniqueBins);
        end
    else
        % Just use available depths
        depths = allDepths;
        if length(depths) < numUniqueBins
            fprintf('Warning: Not enough depths (%d) for all bins (%d).\n', length(depths), numUniqueBins);
            % Fill with artificially spaced depths if needed
            depths = [depths, linspace(depths(end)+1, depths(end)+numUniqueBins-length(depths), numUniqueBins-length(depths))];
        end
    end
    
    % Ensure we have correct number of depths
    depths = depths(1:numUniqueBins);
    fprintf('Using %d unique depths for velocity profile\n', length(depths));
    
    % Convert bin data to numeric values
    binDataNumeric = str2double(table2array(binData));
    
    %% Extract time window of interest
    % Use entire time series by default
    timeWindows = {[min(DMY), max(DMY)]};
    
    % Get time window limits
    startTime = timeWindows{1}(1);
    endTime = timeWindows{1}(2);
    
    % Find indices within this time window
    timeIdx = DMY >= startTime & DMY <= endTime;
    
    % Extract magnitude and direction data for this time window
    magData = binDataNumeric(timeIdx, magBins);
    dirData = binDataNumeric(timeIdx, dirBins);
    DMY = DMY(timeIdx);
    
    % Verify dimensions match
    if size(magData, 2) ~= length(depths)
        error('Dimension mismatch: %d magnitude columns but %d depths', size(magData, 2), length(depths));
    end
    
    % Build statistics structure for this deployment
    stats = struct();
    stats.startTime = startTime;
    stats.endTime = endTime;
    stats.numProfiles = sum(timeIdx);
end

function [T, DistT] = readSpecificFile(filepath, keyword, keyword2, keywordbin, sheetName)
    %% READSPECIFICFILE - Reads ADCP data from a specific Excel file and sheet
    %
    % This function is a simplified version of FolderReadCSV for reading a single file.
    % It extracts ADCP data by finding specific keywords in the headers.
    %
    % Inputs:
    %   filepath   - Full path to the Excel file
    %   keyword    - Keyword to identify velocity column headers (e.g., "nor")
    %   keyword2   - Keyword to identify distance/depth data (e.g., "distance")
    %   keywordbin - Keyword to identify bin numbers (e.g., "bin")
    %   sheetName  - Name or index of the sheet to read
    %
    % Outputs:
    %   T     - Table containing ADCP data with processed headers
    %   DistT - Table containing depth/bin information
    
    try
        %% Read Excel file
        % Read the Excel file into a table using the specified sheet
        opts = detectImportOptions(filepath, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
        opts = setvaropts(opts, opts.VariableNames, 'Type', 'string'); % Read everything as text
        T = readtable(filepath, opts); % Read the table
    
        %% Find keyword rows for header detection
        keylength = strlength(keyword);
        keylength2 = strlength(keyword2);
        keylengthbin = strlength(keywordbin);
        Tcell = string(table2cell(T));
        
        % Find rows containing our keywords
        rowIdx = find(any(strncmpi(Tcell, keyword, keylength), 2), 1);
        [rowIdxdist,colIdydist] = find(any(strncmpi(Tcell, keyword2, keylength2), 2));
        [rowIdxbin,colIdybin] = find(any(strncmpi(Tcell, keywordbin, keylengthbin), 2));
        
        % Verify that all required keywords were found
        if isempty(rowIdx) || isempty(rowIdxdist) || isempty(rowIdxbin)
            error('Missing required keywords in file');
        end
        
        %% Extract variable names from headers
        rawNames1 = string(Tcell(rowIdx, :));       % First header row
        rawNames2 = string(Tcell(rowIdx + 1, :));   % Second header row (often units)
        
        % Get distance and bin information
        dist1 = string(Tcell(rowIdxdist, :));       % Distance values
        binNumbers = string(Tcell(rowIdxbin,:));    % Bin numbers
        
        %% Process depth and bin data
        if colIdydist == colIdybin
            % Distance and bin data are aligned - create distance table
            DistBin = [binNumbers; dist1];
            
            % Remove any missing values
            if any(ismissing(DistBin),'all')
                DistBin = rmmissing(DistBin,2);
            end
            
            % Get bin variable names
            DistBinVars = DistBin(:,1);
            [uniquebinNames, idx_bin] = unique(DistBinVars, 'stable');
            duplicatebinIdx = setdiff(1:numel(DistBinVars), idx_bin);
            
            % Create distance table and handle bin name duplicates
            if isempty(duplicatebinIdx)
                DistT = table(DistBin);
                DistBinVars = matlab.lang.makeValidName(DistBinVars);
                DistT.Properties.RowNames = DistBinVars;
            else
                warning('Issue with bin names - continuing anyway');
                DistT = table(DistBin);
                DistBinVars = matlab.lang.makeValidName(DistBinVars);
                DistT.Properties.RowNames = DistBinVars;
            end
        else
            warning('Distances and bins not aligned - continuing anyway');
            DistT = table();
        end
        
        %% Create combined variable names for the data table
        % Combine header rows for compound variable names
        newVarNames = rawNames1;
        for i = 1:length(rawNames2)
            if i <= length(newVarNames)
                if newVarNames(i) == ""
                    newVarNames(i) = rawNames2(i);
                elseif i <= length(rawNames2) && rawNames2(i) ~= ""
                    newVarNames(i) = newVarNames(i) + "_" + rawNames2(i);
                end
            end
        end
        
        % Replace missing names with numeric placeholders
        for i = 1:length(newVarNames)
            if ismissing(newVarNames(i))
                newVarNames(i) = string(i);
            end
        end

        % Make variable names MATLAB-compatible
        newVarNames = matlab.lang.makeValidName(newVarNames);
        
        % Handle duplicate variable names by appending bin numbers
        [uniqueNames, idx] = unique(newVarNames, 'stable');
        duplicateIdx = setdiff(1:numel(newVarNames), idx);

        for i = duplicateIdx
            if i <= numel(binNumbers) && ~ismissing(binNumbers(i))
                newVarNames(i) = newVarNames(i) + "_" + binNumbers(i);
            else
                newVarNames(i) = newVarNames(i) + "_" + i;
            end
        end

        %% Process data rows and apply variable names
        % Remove header rows from the data
        T(1:rowIdx, :) = [];
        
        % Assign variable names to the table
        if length(newVarNames) >= width(T)
            T.Properties.VariableNames = newVarNames(1:width(T));
        else
            % Extend newVarNames if needed
            while length(newVarNames) < width(T)
                newVarNames = [newVarNames, "Column_" + (length(newVarNames) + 1)];
            end
            T.Properties.VariableNames = newVarNames;
        end

        %% Convert date/time columns
        % Identify and convert date columns to datetime objects
        dateVars = contains(T.Properties.VariableNames, "Date");
        for i = find(dateVars)
            try
                T.(T.Properties.VariableNames{i}) = datetime((T.(T.Properties.VariableNames{i})), 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
            catch
                % If conversion fails, keep as string
                fprintf('Warning: Could not convert date column %s to datetime\n', T.Properties.VariableNames{i});
            end
        end
        
        %% Clean up data table
        % Remove entirely empty rows and columns
        if any(ismissing(T), 'all')
            T = rmmissing(T, 2, 'MinNumMissing', height(T));
            T = rmmissing(T, 1, 'MinNumMissing', width(T));
        end
        
    catch e
        % If any error occurs, propagate it with useful message
        error('Error processing file: %s', e.message);
    end
end

function [logFormula, uStar, z0] = createLogLawProfile(velocityProfile, depths, labelText)
    % LOG LAW: u(z) = (u*/κ) * ln(z/z0) where κ = 0.41 (von Karman constant)
    
    % Convert to column vectors if needed
    velocityProfile = velocityProfile(:);
    depths = depths(:);
    
    % Remove any zero or negative depths (can't take log of these)
    validIdx = depths > 0;
    depthsValid = depths(validIdx);
    velValid = velocityProfile(validIdx);
    
    % Only fit to bottom portion (lower 80% of water column)
    maxDepthToUse = max(depthsValid) * 0.8;
    fitIdx = depthsValid <= maxDepthToUse;
    depthsFit = depthsValid(fitIdx);
    velFit = velValid(fitIdx);
    
    % Set up log law as linear equation: u = a * ln(z) + b
    % where a = u*/κ and b = -u*/κ * ln(z0)
    X = [ones(length(depthsFit), 1), log(depthsFit)];
    coeffs = X \ velFit; % Linear least squares
    b = coeffs(1);
    a = coeffs(2);
    
    % Extract parameters
    kappa = 0.41; % von Karman constant
    uStar = a * kappa;
    z0 = exp(-b/a);
    
    % Create function handle for the log law
    logFormula = @(z) (uStar/kappa) * log(z/z0);
    
    % Plot the original data and fit
    plot(velocityProfile, depths, 'o', 'DisplayName', 'Data');
    hold on;
    
    % Plot the fit on a finer grid
    zFine = logspace(log10(max(0.001, min(depthsValid))), log10(max(depthsValid)), 100);
    vFine = logFormula(zFine);
    plot(vFine, zFine, '-', 'LineWidth', 2, 'DisplayName', 'Log law fit');
    
    % Highlight fitted region
    plot(velFit, depthsFit, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', 'Fitted region');
    
    % Set plot properties
    grid on;
    xlabel('Velocity (m/s)', 'FontSize', 12);
    ylabel('Distance from bed (m)', 'FontSize', 12);
    title(sprintf('%s Profile - Log Law Fit', labelText), 'FontSize', 12);
    legend('Location', 'best');
    
    % Add text with fit parameters
    text(0.05, 0.15, sprintf('u_* = %.4f m/s\nz_0 = %.6f m', uStar, z0), ...
         'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', [1 1 1 0.7]);
end

function [polyCoeffs, velocityFormula] = createPolynomialFit(velocityProfile, depths, labelText, polynomialOrder)
    % Fit polynomial to velocity profile
    
    % Convert to column vectors if needed
    velocityProfile = velocityProfile(:);
    depths = depths(:);
    
    % Fit polynomial
    polyCoeffs = polyfit(depths, velocityProfile, polynomialOrder);
    
    % Create function handle for simulation
    velocityFormula = @(z) polyval(polyCoeffs, z);
    
    % Plot original data and fit
    plot(velocityProfile, depths, 'o', 'DisplayName', 'Data');
    hold on;
    
    % Plot the fit on a finer grid
    zFine = linspace(min(depths), max(depths), 100);
    plot(velocityFormula(zFine), zFine, '-', 'LineWidth', 2, 'DisplayName', 'Polynomial fit');
    
    % Set plot properties
    grid on;
    xlabel('Velocity (m/s)', 'FontSize', 12);
    ylabel('Distance from bed (m)', 'FontSize', 12);
    title(sprintf('%s Profile - Polynomial Fit (Order %d)', labelText, polynomialOrder), 'FontSize', 12);
    legend('Location', 'best');
end

function [velocityFormula, alpha, vSurface] = createPowerLawProfile(velocityProfile, depths, labelText)
    % Power law: v(z) = vSurface * (z/zSurface)^alpha
    
    % Convert to column vectors if needed
    velocityProfile = velocityProfile(:);
    depths = depths(:);
    
    % Extract near-surface velocity as reference
    vSurface = velocityProfile(end);
    zSurface = depths(end);
    
    % Remove any zero depths (can't use in power law)
    validIdx = depths > 0;
    depthsValid = depths(validIdx);
    velValid = velocityProfile(validIdx);
    
    % Find alpha (power law exponent) using least squares
    alpha = fminsearch(@(a) sum((velValid - vSurface*(depthsValid/zSurface).^a).^2), 1/7);
    
    % Create function handle for simulation
    velocityFormula = @(z) vSurface * (z/zSurface).^alpha;
    
    % Plot original data and fit
    plot(velocityProfile, depths, 'o', 'DisplayName', 'Data');
    hold on;
    
    % Plot the fit on a finer grid
    zFine = linspace(max(0.001, min(depths)), max(depths), 100);
    plot(velocityFormula(zFine), zFine, '-', 'LineWidth', 2, 'DisplayName', 'Power law fit');
    
    % Set plot properties
    grid on;
    xlabel('Velocity (m/s)', 'FontSize', 12);
    ylabel('Distance from bed (m)', 'FontSize', 12);
    title(sprintf('%s Profile - Power Law Fit', labelText), 'FontSize', 12);
    legend('Location', 'best');
    
    % Add text with fit parameters
    text(0.05, 0.15, sprintf('v_surface = %.4f m/s\nalpha = %.4f', vSurface, alpha), ...
         'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', [1 1 1 0.7]);
end