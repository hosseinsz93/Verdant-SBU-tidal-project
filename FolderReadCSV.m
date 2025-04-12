function [allData, allDistT] = FolderReadCSV(folderpath, keyword, keyword2, keywordbin, sheetName)
    %% FOLDERREADCSV - Extracts ADCP data from Excel files in a folder
    % By: Jonathan Craig
    % Date: 26 February, 2025
    %
    % Last modified: Hossein Seyedzadeh
    % Date: 10 April, 2025
    %
    % Syntax:
    %   [allData, allDistT] = FolderReadCSV(folderpath, keyword, keyword2, keywordbin)
    %   [allData, allDistT] = FolderReadCSV(folderpath, keyword, keyword2, keywordbin, sheetName)
    %
    % Inputs:
    %   folderpath - Path to folder containing Excel files with ADCP data
    %   keyword    - Keyword to identify velocity column headers (e.g., "nor")
    %   keyword2   - Keyword to identify distance/depth data (e.g., "distance")
    %   keywordbin - Keyword to identify bin numbers (e.g., "bin")
    %   sheetName  - (Optional) Sheet name or index to read (default: 1)
    %
    % Outputs:
    %   allData    - Table or cell array of tables containing ADCP velocity data
    %   allDistT   - Table or cell array of tables containing depth/bin information
    %
    % Example:
    %   [data, distTable] = FolderReadCSV(pwd, "nor", "distance", "bin", "Station B-2");
    %
    % Notes:
    %   - The function searches Excel files for specific keywords in headers
    %   - Handles multiple Excel files and returns a cell array if >1 file processed
    %   - Automatically converts date/time columns to MATLAB datetime objects
    %   - Skips temporary Excel files (starting with ~$)
    
    arguments
        folderpath {mustBeText}
        keyword {mustBeText}
        keyword2 {mustBeText}
        keywordbin {mustBeText}
        sheetName = 1  % Default to first sheet if not specified
    end

    %% Input validation
    % Validate sheet name manually since arguments block can't handle OR logic
    if ~(isstring(sheetName) || ischar(sheetName) || isnumeric(sheetName))
        error('Sheet specifier must be a string, character array, or numeric value');
    end

    %% Find Excel files in the specified folder
    % Filter out temporary Excel files (files starting with ~$)
    files = dir(fullfile(folderpath, '*.xlsx'));
    files = files(~startsWith({files.name}, '~$'));
    
    if isempty(files)
        error('No valid Excel files found in the folder');
    end
    
    %% Initialize data structures
    % Create cell arrays to store data from each file
    allData = cell(numel(files), 1);
    allDistT = cell(numel(files), 1);
    fileNames = cell(numel(files), 1);
    
    fprintf('Found %d Excel files to process\n', numel(files));
    
    %% Process each Excel file
    for k = 1:numel(files)
        % Record the file name without extension
        [~, fileNames{k}] = fileparts(files(k).name);
        fprintf('Processing file %d/%d: %s\n', k, numel(files), files(k).name);
        
        % Construct full file path
        filepath = fullfile(folderpath, files(k).name);

        try
            %% Read Excel file data
            % Read the Excel file into a table, using the specified sheet
            opts = detectImportOptions(filepath, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
            opts = setvaropts(opts, opts.VariableNames, 'Type', 'string'); % Read everything as text
            T = readtable(filepath, opts); % Read the table
        
            %% Find keyword rows for header detection
            % Search for rows containing velocity, distance, and bin keywords
            keylength = strlength(keyword);
            keylength2 = strlength(keyword2);
            keylengthbin = strlength(keywordbin);
            Tcell = string(table2cell(T));
            rowIdx = find(any(strncmpi(Tcell, keyword, keylength), 2), 1);
            [rowIdxdist,colIdydist] = find(any(strncmpi(Tcell, keyword2, keylength2), 2));
            [rowIdxbin,colIdybin] = find(any(strncmpi(Tcell, keywordbin, keylengthbin), 2));
            
            % Skip this file if required keywords are missing
            if isempty(rowIdx) || isempty(rowIdxdist) || isempty(rowIdxbin)
                warning('Skipping file %s: Missing required keywords', files(k).name);
                continue;
            else
                fprintf('All keywords are valid for file: %s\n', files(k).name);
            end
            
            %% Extract variable names and distance data
            % Get variable names from header rows
            rawNames1 = string(Tcell(rowIdx, :));
            rawNames2 = string(Tcell(rowIdx + 1, :));
            
            % Get distance and bin information
            dist1 = string(Tcell(rowIdxdist, :));
            binNumbers = string(Tcell(rowIdxbin,:));
            
            %% Process depth and bin data
            if colIdydist == colIdybin
                % Distance and bin data are aligned - create distance table
                DistBin = [binNumbers; dist1];
                if any(ismissing(DistBin),'all')
                    DistBin = rmmissing(DistBin,2);
                end
                DistBinVars = DistBin(:,1);
                [uniquebinNames, idx_bin] = unique(DistBinVars, 'stable');
                duplicatebinIdx = setdiff(1:numel(DistBinVars), idx_bin);
                
                % Create distance table and handle bin name duplicates
                if isempty(duplicatebinIdx)
                    DistT = table(DistBin);
                    DistBinVars = matlab.lang.makeValidName(DistBinVars);
                    DistT.Properties.RowNames = DistBinVars;
                else
                    warning('Issue with bin names in file %s - continuing anyway', files(k).name);
                    DistT = table(DistBin);
                    DistBinVars = matlab.lang.makeValidName(DistBinVars);
                    DistT.Properties.RowNames = DistBinVars;
                end
            else
                warning('Distances and bins not aligned in file %s - continuing anyway', files(k).name);
                DistT = table();
            end
            
            %% Create variable names for the data table
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
            
            %% Store processed data in output arrays
            allData{k} = T;
            allDistT{k} = DistT;
            
            % Add metadata to the table
            allData{k}.Properties.UserData.FileName = fileNames{k};
            allData{k}.Properties.UserData.DeploymentPeriod = k;
            allData{k}.Properties.UserData.FilePath = filepath;
            
        catch e
            warning('Error processing file %s: %s', files(k).name, e.message);
            % Continue with next file
        end
    end
    
    %% Prepare return values
    % Remove empty entries if any files were skipped
    validIdx = ~cellfun(@isempty, allData);
    allData = allData(validIdx);
    allDistT = allDistT(validIdx);
    fileNames = fileNames(validIdx);
    
    % If only one file was processed, return single tables instead of cell arrays
    if length(allData) == 1
        allData = allData{1};
        allDistT = allDistT{1};
        fprintf('Processed 1 file: %s\n', fileNames{1});
    else
        fprintf('Successfully processed %d files\n', length(allData));
    end
end