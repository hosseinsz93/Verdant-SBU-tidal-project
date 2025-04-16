function CvelAvg = ExtractCvelAvg(workfolder, sheetName)
    % EXTRACTCVELAVG - Extracts depth-averaged velocity data from ADCP Excel files
    %
    % This function reads ADCP Excel files and extracts depth-averaged velocity 
    % components (east, north, magnitude, direction). It automatically converts
    % velocity units from mm/s to m/s.
    %
    % Syntax:
    %   CvelAvg = ExtractCvelAvg(workfolder)
    %   CvelAvg = ExtractCvelAvg(workfolder, sheetName)
    %
    % Inputs:
    %   workfolder - Path to folder containing ADCP Excel files
    %   sheetName  - (Optional) Excel sheet name or number (default: 1)
    %
    % Outputs:
    %   CvelAvg - Numeric array of depth-averaged velocity measurements
    %             Column 1: East velocity (m/s)
    %             Column 2: North velocity (m/s)
    %             Column 3: Velocity magnitude (m/s)
    %             Column 4: Direction (degrees from true North)
    %
    % Notes:
    %   - If multiple files are found, only the first file is processed
    %   - Velocity values are automatically converted from mm/s to m/s
    %
    % Author: Hossein Seyedzadeh
    % Date: April 16, 2025
    
    % Default sheet if not specified
    if nargin < 2
        sheetName = 1;
    end
    
    % Set keywords for data identification
    keyword = "nor";       % For velocity components
    keyword2 = "distance"; % For distance/depth data
    keywordbin = "bin";    % For bin information
    
    % Read the Excel file using FolderReadCSV
    [C, ~] = FolderReadCSV(workfolder, keyword, keyword2, keywordbin, sheetName);
    
    % Check if C is a cell array (multiple files processed)
    if iscell(C)
        fprintf('Processing first file only (out of %d files)\n', numel(C));
        C = C{1}; % Use only the first file's data
    end
    
    % Now C is guaranteed to be a table
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