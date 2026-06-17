function [files_num,VarNames,C,DMY,dateInterval,convertFactor,siteID,seadepths,depth_units] = FolderReadCSV(folderpath)
    
    arguments
        folderpath {mustBeText}
    end

    %% Preprocessing 1
    % Get list of all CSV files in the folder
    files = dir(fullfile(folderpath, '*.csv'));
    files_num = numel(files);

    dataCells = cell(numel(files), 1);
    dateCells = cell(numel(files), 1);
    filenames = strings(files_num,1);
    
    prev = struct();
    
    % This loop reads the info of files' names, and then their contents of
    % datetime, speed, and direction. These data populate the appropriate
    % arrays. Also, it reads the units automatically whereas it takes the 
    % user's input to concatenate datetime values for datasets longer than
    % one month
    for k = 1:files_num
        % Construct full file path
        filepath = fullfile(folderpath, files(k).name); % D = dir('NAME') returns the results in an M-by-1 structure with the fields: name        -- Filename
        filenames(k) = files(k).name;
        
        % Read the CSV file into a table
        
        T = readCSVFile(filepath);
        
        % Extract metadata
        meta = extractMetadata(T,files(k).name);
        
        if k > 1 && exist("convertFactor",'var')
            valVarNames(meta.VarNames,prev.VarNames,files(k-1).name,files(k).name);
        else
            convertFactor = unitConverter(meta.fileUnits);
        end
        
        [dateData,numericData] = extractData(T);

        if k > 1 && isfield(prev,'dateData') % k > 1 && exist('dateInterval','var')
            [dateOut, dateOut_unique, dataOut, uniquedatetime_Length] = processTimeSeries(prev, dateData, numericData, meta.fileDepth, files(k), files(k-1), dateInterval);
            if ~isequal(dateOut_unique,prev.dateOut_unique)
                warning("There is mismatch between unique dates.");
                dateOut_uniqueSwitch = input("Write '1' for the new set of dates, or '0' to keep the old set of dates: ");
                switch (dateOut_uniqueSwitch)
                    case 1
                        prev.dateOut_unique = dateOut_unique;
                    case 0
                        dateOut_unique = prev.dateOut_unique;
                    otherwise
                        error("Invalid input. Input must be numeric '1' or '0'.");
                end
            end
        else
            [dateOut, dataOut, dateInterval, uniquedatetime_Length, dateOut_unique] = initializeTimeSeries(dateData, numericData);
            prev.dateOut_unique = dateOut_unique;
        end
    
        % Temporary data storage
        dateCells{k} = dateOut;
        dataCells{k} = dataOut;
        
        % Update of previous structure
        prev = updatePrev(prev,meta,dateData,numericData,uniquedatetime_Length);
    
    end
    
    %% Preprocessing 2
    
    VarNames = prev.VarNames;
    % Extract site ID from file
    siteID = extractSiteID(filenames);
    % Gather water depths into array
    [seadepths,depth_units] = findDepth(filenames,files_num);
    % Assignment of concatenated datetime values
    DMY = dateOut_unique;
    
    dataCells = alignDataLengths(dataCells, length(DMY));
    C = horzcat(dataCells{:});

end

function T = readCSVFile(filepath)

    T = readtable(filepath,'VariableNamingRule','preserve');
    
end

function meta = extractMetadata(T,filename)

        VarNames = T.Properties.VariableNames;
        
        if ~any(contains(VarNames,"Speed"))
            error("%s does not contain 'Speed' column.",filename);
        end
        
        for i = 1:numel(VarNames)
            if contains(VarNames(i),"Speed")
                % This line extracts the units from the "Speed" column
                % VarKey(i) = VarNames(i);
                fileUnits = extract(string(extractAfter(VarNames(i),"Speed")),lettersPattern);
                if contains(fileUnits,["cm","s"])
                    fileUnits = join(fileUnits,"/");
                end
                break
            end
        end
        
        fileDepth = string(extractBetween(filename,"_","-",'Boundaries','exclusive'));

        
        meta.VarNames = VarNames;
        meta.fileUnits = fileUnits;
        meta.fileDepth = fileDepth;

end

function valVarNames(current,previous,filePrev,fileCurr)

    if ~isequal(current, previous)
        error("Column mismatch between %s and %s.", filePrev, fileCurr);
    end
    
end

function [dateData, numericData] = extractData(T)

    dateData = datetime(T{:,1});
    numericData = T{:,2:end};
    
end

function [dateOut, dataOut, dateInterval, uniquedatetime_Length, dateOut_unique] = initializeTimeSeries(dateData, numericData)
    
    dateDiff = diff(dateData);

    if all(dateDiff == dateDiff(1))
        dateInterval = dateDiff(1);
    else
        warning("The is uneven date interval within current dataset, and so first interval is chosen by default.");
        dateInterval = dateDiff(1);
    end

    dateOut = dateData;
    dataOut = numericData;
    dateOut_unique = dateData; % Preliminary value for processTimeSeries
    uniquedatetime_Length = length(dateData);
    
end

function [dateOut, dateOut_unique, dataOut, uniquedatetime_Length] = processTimeSeries(prev, dateData, numericData, fileDepth, fileCurr, filePrev, dateInterval)
    
    if dateData(1)>prev.dateData(end) && ...
            abs(dateInterval - (dateData(1)-prev.dateData(end))) < seconds(1) && ...
            isequal(fileDepth, prev.fileDepth) % For datasets within timespans longer than 1 month
        % Condition for concatenation of datarange longer than 1 month at same water depth
        
        uniquedatetime = unique(cat(1,prev.dateData,dateData));
        uniquedatetime_Length = length(uniquedatetime);
        
        if uniquedatetime_Length ~= prev.uniquedatetime_Length
            timeDiff = abs(uniquedatetime_Length-prev.uniquedatetime_Length);
            timeDiff_datetime = string(timeDiff*dateInterval);
            warning("There is difference of %d date entries between files %s and %s.\nThe previous file ranges from %s to %s.\nThe current file ranges from %s to %s.\nThis number of entries correspond to a time difference of %s.",...
                timeDiff,fileCurr.name,filePrev.name,string(prev.dateData(1)),string(prev.dateData(end)),string(dateData(1)),string(dateData(end)),timeDiff_datetime);
            user_concatenate = input("Do you approve the concatenation of these datasets? Please write '1' to approve or '0' to disapprove: ");
            if user_concatenate
                dateOut_unique = uniquedatetime(1:uniquedatetime_Length);
            else
                dateOut_unique = uniquedatetime(1:prev.uniquedatetime_Length);
            end
        end

        dateOut = uniquedatetime;
        dataOut = [prev.numericData; numericData];

    elseif isequal(dateData,prev.dateData) % && (dateData(1)==prev.dateData(1)) % For datasets within timespans of 1 month

        dateOut = prev.dateData;
        dataOut = numericData;
    
    elseif ~isequal(dateData,prev.dateData) && isequal(fileDepth,prev.fileDepth)
        
        warning("Date mismatch between files %s and %s.", filePrev.name, fileCurr.name);
        dateOut = dateData;
        dataOut = numericData;

    else
        
        dateOut = dateData;
        dataOut = numericData;
        
    end
    
    if ~exist('uniquedatetime_Length','var')
        uniquedatetime_Length = prev.uniquedatetime_Length;
    end
    if ~exist('dateOut_unique','var')
        dateOut_unique = prev.dateOut_unique;
    end
    
end

function prev = updatePrev(prev, meta, dateData, numericData, uniquedatetime_Length)

    prev.VarNames = meta.VarNames;
    prev.fileDepth = meta.fileDepth;
    prev.dateData = dateData;
    prev.numericData = numericData;
    prev.uniquedatetime_Length = uniquedatetime_Length;

end

function [waterDepth,depth_units] = findDepth(filenames,files_num)
	
    arguments
        filenames {mustBeText}
        files_num {mustBeNumeric,mustBePositive}
    end
    
    fileDepths = extractBetween(filenames,"_","-",'Boundaries','exclusive');
    depthUnits = extract(fileDepths,lettersPattern);
    if all(depthUnits == depthUnits(1,:))
        depth_units = depthUnits(1,1); % To be used for plotting of depth either as meters (m) or as feet (m)
    else
        error("waterDepth:incorrectFormat","Error in file format. \nFile format must list same units after site ID. \nAcceptable formats are the following: \nLIS1001_05m76cm \nLIS1001_18ft09df")
    end
    
    depthDigits = str2double(extract(fileDepths,digitsPattern));
    waterDepth = depthDigits(:,1)+depthDigits(:,end)./100; % Same as m + cm or ft + 100ths of ft
    
    % For long-running dataset with repeated water depths
    waterDepth = unique(waterDepth);
    
end

function convertFactor = unitConverter(fileUnits)

    arguments
        fileUnits {mustBeText}
    end
    
    sprintf("Setup of conversion factor from extracted units to meters per second (m/s) for all files' tidal analysis");
    if contains(fileUnits,"knots")
        convertFactor = 0.51444448824222; % knots to m/s
    elseif  contains(fileUnits,"cm/sec")
        convertFactor = 0.01; % cm/s to m/s
    else
        % error("\n  In CO-OPS CSV file, 'Speed' column must be defined in 'knots' or 'cm/s'. Please check format.\n");
        warning("File defines water speed neither in knots or cm/s. Conversion factor of 1 is given by default. Please check file's units.");
        convertFactor = 1;
    end
    
end

function siteID = extractSiteID(filenames)

    siteName = extract(filenames,lettersPattern);
    siteNum = extract(filenames,digitsPattern);
    
    if all(siteName(1,1) == siteName(:,1)) && all(siteNum(1,1) == siteNum(:,1))
        siteID = siteName(1)+siteNum(1);
    else
        warning("There is mismatch of site ID within provided data.");
        siteID = inputID("User-defined ID requested for plotting: ",'s');
    end
    
end

function siteID = inputID()

    user_decision = input('Do you want to continue by defining the site ID? (Y/N): ', 's');
    if isequal(upper(user_decision), 'N')
        disp('Exiting from FolderReadCSV. Recommendation to revise data in files.');
        return; % Ends the function immediately
    else
        siteID = input("\n Please define the site ID to appear in plots (e.g., LIS1001): ", 's');
    end
    disp('Continuing with user');
    
end

function dataCells = alignDataLengths(dataCells, targetLength)

    for k = 1:numel(dataCells)
        
        if isempty(dataCells{k})
            continue
        end
        
        dataLength = size(dataCells{k},1);
        
        if dataLength < targetLength
            dataCells{k} = []; % data(k) = [] removes the element outright, but this makes an error in the loop's indices
        elseif dataLength >= targetLength
            dataCells{k} = dataCells{k}(1:targetLength,:);
        else
            continue
        end
    end
    
    dataCells(cellfun(@isempty, dataCells)) = [];
    
end

% % Original function for testing
% function M = FolderReadCSV(folderpath)
% arguments
%     folderpath {mustBeText}
% %     range {mustBeText}
% end
% 
% folder = folderpath; %'C:\Users\jonny\OneDrive\Pictures\Documents\MATLAB\Test\'
% d = dir(folder);
% e = {d.name};
% f=e(~cellfun(@isempty,regexp(e,'.+(?=\.xlsx)','match')));
% 
% for k=1:numel(f)-1
% %     if ~exist('range','var')
% %         data{k,1}=xlsread(f{k}); %'A1:C20'
% %     else
%     data{k,1}=xlsread(f{k},range);
% end
% M=cell2mat(data);
% xlswrite('new_file',M);
% end