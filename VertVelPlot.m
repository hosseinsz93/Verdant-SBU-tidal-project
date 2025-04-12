%==========================================================================
%% VERTVELPLOT - Generate velocity profile plots from ADCP data
%
% Goal: Plots of ADCP Info at Orient Point from 10/23/24 to 12/03/24
% By: Jonathan Craig
% Date: 26 February, 2025
%
% Last modified: Hossein Seyedzadeh
% Date: 10 April 2025
%
% This script:
% 1. Imports ADCP data from Excel files
% 2. Processes velocity magnitude and direction data
% 3. Generates vertical velocity profile plots for each time point
% 4. Creates speed, direction, and 3D vector visualizations
%==========================================================================
%%
clc;
clear;
close all;

% Start timing the simulation
tic;

%% STEP 1: DATA IMPORT
% Define search terms for data extraction
workfolder = pwd;
keyword = "nor";        % For velocity component identification
keyword2 = "distance";  % For depth/distance identification
keywordbin = "bin";     % For bin number identification

% Option 1: To read from the first sheet (default):
% [C, DistT] = FolderReadCSV(workfolder, keyword, keyword2, keywordbin);

% Option 2: To read from a specific sheet by name:
% Read the excel file with specific sheet name
[C, DistT] = FolderReadCSV(workfolder, keyword, keyword2, keywordbin, "Station B-2");

% Check if any data was returned
if isempty(C) || ~istable(C)
    error('No valid data found in Excel files. Check that files exist and contain required keywords.');
end

% Option 3: To read from a specific sheet by number:
% [C, DistT] = FolderReadCSV(workfolder, keyword, keyword2, keywordbin, 3);

% Get table dimensions
[nrow,ncol] = size(C);

%% STEP 2: EXTRACT TIMESTAMPS
% Find date/time columns
DT = contains(C.Properties.VariableNames, "Date");
DT_EST = contains(C.Properties.VariableNames, "EST");

% Combine logical arrays to find any date columns
dateCols = DT | DT_EST; 
lastDateCol = find(dateCols,1,'last');

% Extract datetime values from EST column
for i = find(DT_EST)
    DMY = table2array(C(:,i));
end

% Clean up datetime values
if any(ismissing(DMY),'all')
    DMY = rmmissing(DMY,1);
end

%% STEP 3: EXTRACT VELOCITY DATA
% Separate vertically averaged measurements from bin measurements
Cmod = C(:, lastDateCol+1:end);

% Remove empty rows and columns
absentCols = all(ismissing(Cmod),1);
absentRows = all(ismissing(Cmod),2);
Cmod(:,absentCols) = [];
Cmod(absentRows,:) = [];

% Identify velocity measurements by type
CmodVarNames = Cmod.Properties.VariableNames;
velVars = contains(CmodVarNames, ["Eas","Nor","Mag","Dir"]);
velNames = CmodVarNames(velVars);

% Count unique measurement types
velNums = unique(extractBefore(velNames, "_")); 
velGages = numel(velNums); 

% Separate depth-averaged and bin-specific data
velAvg = Cmod(:,1:velGages);
Cmod = Cmod(:,velGages+1:end);

% Update variable names after removing averaged values
CmodVarNames = Cmod.Properties.VariableNames; 
velVars = contains(CmodVarNames, ["Eas","Nor","Mag","Dir"]);
velNames = CmodVarNames(velVars);

% Extract units from variable names
unitLengthGages = string(unique(extractBetween(CmodVarNames,'__','_','Boundaries','exclusive')));
unitTimeGages = string(unique(extractBetween(CmodVarNames,'_','_','Boundaries','exclusive')));

%% STEP 4: EXTRACT DEPTH DATA
% Get row names from the distance table
DistBinVars = string(DistT.Properties.RowNames);

% Find distance/depth information
distIdx = find(any(strncmpi(DistBinVars, keyword2, strlength(keyword2)), 2), 1);
unitDist = extractBetween(DistBinVars(distIdx),'_','_','Boundaries','exclusive'); 

% Extract numeric depth values
DistT_ind = DistT(distIdx,:);
DistA = str2double(table2array(DistT_ind));

% Clean up depth values
if any(ismissing(DistA),'all')
    DistA = rmmissing(DistA,2);
end

%% STEP 5: CALCULATE VELOCITY COMPONENTS
% Identify magnitude and direction columns
mag = contains(CmodVarNames, ["Mag"]);
dir = contains(CmodVarNames, ["Dir"]);
magdir = cat(1,mag,dir);

% Convert string data to numeric values
Cvel = str2double(table2array(Cmod));
ncol_vel = size(Cvel,2);

% Remove rows with missing data
if any(ismissing(Cvel),'all')
    Cvel = rmmissing(Cvel,1,'MinNumMissing',ncol_vel);
end
nrow_vel = size(Cvel,1);

%% STEP 6: CONVERT UNITS
% Define unit conversion options
lengthUnits = ["m","mm","cm","knot","ft","in"];
timeUnits = ["s","min","hour","day"];

% Identify units in variable names
varLengthUnits = contains(CmodVarNames, lengthUnits);
varTimeUnits = contains(CmodVarNames, timeUnits);
lengthInd = contains(unitLengthGages, lengthUnits);
timeInd = contains(unitTimeGages, timeUnits);
unitLength = unitLengthGages(lengthInd);
unitTime = unitTimeGages(timeInd);
unitVel = join([unitLength,unitTime],'/');

% Convert all velocity measurements to standard units (m/s)
for i = 1:nrow_vel
   for j = 1:ncol_vel
       if varLengthUnits(j)
           % Convert length units to meters
           Cvel(i,j) = convert_m(Cvel(i,j),unitLength);
           % Convert time units to seconds
           Cvel(i,j) = convert_s(Cvel(i,j),unitTime);
       end
   end
end

%% STEP 7: GENERATE VELOCITY PROFILE PLOTS
% Extract unique depths and identify column indices
depths = unique(DistA);
magInd = find(mag);
dirInd = find(dir);
magFreq = length(magInd);
dirFreq = length(dirInd);
units = [unitDist,unitVel];

% Verify magnitude and direction columns match
if magFreq ~= dirFreq
    error("Number of columns for magnitude does not agree with number of columns for direction");
end

% Initialize velocity arrays
[nrow_dmy,~] = size(DMY);
velMag = zeros(magFreq,2);
velDir = zeros(dirFreq,2);

%% STEP 8: PROCESS EACH TIME POINT WITH PROGRESS BAR
% Create progress bar
h = waitbar(0, 'Processing velocity plots...', 'Name', 'ADCP Data Processing');

% Create output folders if they don't exist
speedFolder = fullfile(workfolder, "SpeedProfiles");
dirFolder = fullfile(workfolder, "DirProfiles");
velFolder = fullfile(workfolder, "VelProfiles");
if ~exist(speedFolder, 'dir'), mkdir(speedFolder); end
if ~exist(dirFolder, 'dir'), mkdir(dirFolder); end
if ~exist(velFolder, 'dir'), mkdir(velFolder); end

% Check dimensions of timestamps and velocity data
if nrow_dmy == nrow_vel
    % Process each time point
    for i = 1:nrow_dmy
        % Update progress bar (less frequently to reduce flickering)
        if mod(i, 10) == 0 || i == nrow_dmy
            waitbar(i/nrow_dmy, h, sprintf('Processing velocity plots: %d/%d (%.1f%%)', i, nrow_dmy, 100*i/nrow_dmy));
        end
        
        % Extract velocity data for this time point
        velMag = Cvel(i,magInd);
        velDir = Cvel(i,dirInd);
        
        % Generate plots if data is valid
        if all(velMag) && all(velDir)
            VelPlot_Verdant(velMag, velDir, depths, units, DMY(i), workfolder, "Deployment1");
        end
    end
    
    % Close progress bar when done
    close(h);
else
    % Close progress bar if there's an error
    close(h);
    error("Number of date times (%d) does not agree with number of velocity measurements (%d). Please check data or program for mismatch.", nrow_dmy, nrow_vel);
end

%% UTILITY FUNCTIONS FOR UNIT CONVERSION
% Convert distance to meters
function distance_m = convertDistance_m(distance,unitdist)
    switch lower(unitdist)
        case "m"
            distance_m = distance;
        case "mm"
            distance_m = distance*1000;
        case "cm"
            distance_m = distance*100;
        case "km"
            distance_m = distance/1000;
        case "ft"
            distance_m = distance*0.3048; % 1 ft = 0.3048 m
        case "in"
            distance_m = distance*0.0254; % 1 in = 0.0254 m
        case "mi"
            distance_m = distance*0.0006214; % 1 mi = 0.000621371 m
        otherwise
            error("Unrecognized unit: %s. Accepted units include the following: m, cm, mm, km, mi, ft, and in");
    end
end

% Convert velocity length unit to meters
function velocity_m = convert_m(velocity,unitlength)
    switch lower(unitlength)
        case "m"
            velocity_m = velocity;
        case "mm"
            velocity_m = velocity/1000;
        case "cm"
            velocity_m = velocity/100;
        case "knots"
            velocity_m = velocity*0.514444; % 1 knot = 0.514444 m/s
        case "ft"
            velocity_m = velocity*0.3048; % 1 ft/s = 0.3048 m/s
        case "in"
            velocity_m = velocity*0.0254; % 1 in = 0.0254 m
        otherwise
            error("Unrecognized unit: %s. Accepted units include the following: m, cm, mm, knots, ft, and in");
    end
end

% Convert velocity time unit to seconds
function velocity_s = convert_s(velocity,unittime)
    switch lower(unittime)
        case "s"
            velocity_s = velocity;
        case "min"
            velocity_s = velocity/60;
        case "hour"
            velocity_s = velocity/3600;
        case "day"
            velocity_s = velocity/86400; % 1 day = 86400 s
        otherwise
            error("Unrecognized unit: %s. Accepted units include the following: s, min, hour, day");
    end
end

% Display total execution time
totalTime = toc;
fprintf('Processing completed in %.2f seconds (%.2f minutes)\n', totalTime, totalTime/60);