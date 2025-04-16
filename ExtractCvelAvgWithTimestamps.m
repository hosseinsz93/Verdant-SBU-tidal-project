function [CvelAvg, fileData] = ExtractCvelAvgWithTimestamps(workfolder, sheetName)
    % EXTRACTCVELAVGWITHTIMESTAMPS - Extracts velocity data and timestamps from ADCP files
    %
    % This function reads ADCP Excel files once and extracts both velocity data
    % and original file data containing timestamps. It avoids duplicate file 
    % processing by combining file reading and data extraction operations.
    %
    % Syntax:
    %   [CvelAvg, fileData] = ExtractCvelAvgWithTimestamps(workfolder)
    %   [CvelAvg, fileData] = ExtractCvelAvgWithTimestamps(workfolder, sheetName)
    %
    % Inputs:
    %   workfolder - Path to folder containing ADCP Excel files
    %   sheetName  - (Optional) Excel sheet name or number (default: 1)
    %
    % Outputs:
    %   CvelAvg  - Depth-averaged velocity data (m/s)
    %               Column 1: East velocity (m/s)
    %               Column 2: North velocity (m/s)
    %               Column 3: Velocity magnitude (m/s)
    %               Column 4: Direction (degrees from true North)
    %   fileData - Original file data containing timestamps and other information
    %
    % Notes:
    %   - This function automatically converts velocity units from mm/s to m/s
    %   - Only processes the first file if multiple files are found
    %
    % Author: Hossein Seyedzadeh
    % Date: April 16, 2025

    % Default sheet if not specified
    if nargin < 2
        sheetName = 1;
    end
    
    % Set keywords for data identification
    keyword = "nor";      % For velocity components
    keyword2 = "distance"; % For distance/depth data
    keywordbin = "bin";    % For bin information
    
    % Read the Excel file once
    [fileData, ~] = FolderReadCSV(workfolder, keyword, keyword2, keywordbin, sheetName);
    
    % Process the data directly instead of calling ExtractCvelAvg again
    % This avoids duplicate file reading
    
    % Check if fileData is a cell array (multiple files processed)
    if iscell(fileData)
        fprintf('Processing first file only (out of %d files)\n', numel(fileData));
        C = fileData{1}; % Use only the first file's data
    else
        C = fileData;
    end
    
    % Identify date/time columns
    DT = contains(C.Properties.VariableNames, "Date");
    DT_EST = contains(C.Properties.VariableNames, "EST");
    dateCols = DT | DT_EST;
    lastDateCol = find(dateCols, 1, 'last');
    
    % Extract data after date columns
    Cmod = C(:, lastDateCol+1:end);
    
    % Remove empty rows and columns
    absentCols = all(ismissing(Cmod), 1);
    absentRows = all(ismissing(Cmod), 2);
    Cmod(:, absentCols) = [];
    Cmod(absentRows, :) = [];
    
    % Identify velocity-related variables
    CmodVarNames = Cmod.Properties.VariableNames;
    velVars = contains(CmodVarNames, ["Eas", "Nor", "Mag", "Dir"]);
    velNames = CmodVarNames(velVars);
    
    % Count unique measurement types
    velNums = unique(extractBefore(velNames, "_")); 
    velGages = numel(velNums);
    
    % Extract just the depth-averaged velocity data
    velAvg = Cmod(:, 1:velGages);
    
    % Convert to numeric array
    CvelAvg = str2double(table2array(velAvg));
    
    % Handle missing values
    ncol_velAvg = size(CvelAvg, 2);
    if any(ismissing(CvelAvg), 'all')
        CvelAvg = rmmissing(CvelAvg, 1, 'MinNumMissing', ncol_velAvg);
    end
    
    % Convert velocity values from mm/s to m/s
    % The first three columns are East, North, and Magnitude in mm/s
    % The fourth column is Direction in degrees (no conversion needed)
    if size(CvelAvg, 2) >= 3
        % Convert East, North, and Magnitude from mm/s to m/s
        CvelAvg(:, 1:3) = CvelAvg(:, 1:3) / 1000;
        fprintf('Converted velocity units from mm/s to m/s\n');
    else
        fprintf('Warning: Expected at least 3 velocity columns but found %d\n', size(CvelAvg, 2));
    end
end