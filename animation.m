%% ADCP ANIMATION - Creates animation from sequential velocity profile plots
%
% This script creates an animation by sequentially displaying PNG images of ADCP
% velocity profiles that are sorted chronologically based on timestamps embedded
% in the filenames.
%
% Purpose:
%   - Automatically finds PNG files from a specified directory
%   - Extracts date/time information from filenames
%   - Sorts images chronologically
%   - Creates a video animation showing how flow conditions change over time
%
% Creator: Hossein Seyedzadeh
% Date: April 12, 2025
% Last Modified: April 12, 2025
%==========================================================================
%% STEP 1: INITIALIZATION AND CONFIGURATION
% Clean up workspace and command window
clc;
clear;
close all;

% Animation settings
frameRate = 10;                    % Frames per second (adjust for faster/slower playback)
outputFileName = 'animation.avi';  % Name of output video file
inputDirectory = 'C:\Users\Student\Documents\ADCP\Matlab\DirProfiles';  % Source directory

% Get current directory for reference
current_dir = pwd;
disp(['Looking for PNG files in: ', inputDirectory]);

%% STEP 2: FILE DISCOVERY AND VALIDATION
% Get all PNG files (case-insensitive) from the input directory
files = dir(fullfile(inputDirectory, '*.png'));
filenames = {files.name};

% Verify that files were found and display examples
disp('First few filenames found:');
if ~isempty(filenames)
    % Display a sample of the filenames for verification
    for i = 1:min(5, length(filenames))
        disp(filenames{i});
    end
else
    % No files found - provide troubleshooting tips
    disp('No PNG files found in the specified directory!');
    disp('Possible solutions:');
    disp('  1. Ensure the directory path is correct');
    disp('  2. Check that PNG files exist in that directory');
    disp('  3. Verify file permissions allow reading from that location');
    return;
end

%% STEP 3: EXTRACT DATETIME INFORMATION FROM FILENAMES
% Extract dates and times and create a sortable array
validFiles = true(size(filenames));    % Keep track of which files have valid timestamps
dateNums = zeros(size(filenames));     % Store numeric date values for sorting

% Define the expected pattern for date/time extraction
% For filenames like "Deployment1_2D_VC_01-Apr-2025 09-30-00.png"
pattern = '(\d{2}-\w{3}-\d{4}) (\d{2}-\d{2}-\d{2})';

for i = 1:length(filenames)
    % Extract date/time substring from filename using regular expression
    tokens = regexp(filenames{i}, pattern, 'tokens');
    
    if ~isempty(tokens)
        % Successfully extracted date and time components
        date_str = tokens{1}{1};    % e.g., "01-Apr-2025"
        time_str = tokens{1}{2};    % e.g., "09-30-00"
        
        % Convert to MATLAB datetime and then to serial date for numeric sorting
        % Replace hyphens in time with colons for proper datetime formatting
        dt = datetime([date_str, ' ', strrep(time_str, '-', ':')], ...
                      'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
        dateNums(i) = datenum(dt);
    else
        % If pattern doesn't match, mark file as invalid
        validFiles(i) = false;
        warning('Pattern did not match for file: %s', filenames{i});
    end
end

%% STEP 4: SORT FILES CHRONOLOGICALLY
% Filter out invalid files (those that didn't match the datetime pattern)
filenames = filenames(validFiles);
dateNums = dateNums(validFiles);

% Sort the filenames based on date numbers
[~, idx] = sort(dateNums);
sorted_filenames = filenames(idx);

% Report sorting results
disp(['Successfully sorted ', num2str(length(sorted_filenames)), ' PNG files chronologically']);

%% STEP 5: CREATE ANIMATION
% Create a figure for displaying the animation frames
figure('Position', [100, 100, 800, 600], 'Name', 'ADCP Profile Animation');

% Initialize the video writer with specified frame rate
v = VideoWriter(outputFileName);
v.FrameRate = frameRate;
open(v);

% Display progress bar information
fprintf('Creating animation: [');
progressPoints = round(linspace(1, length(sorted_filenames), 20));
nextPoint = 1;

%% STEP 6: PROCESS IMAGES AND BUILD ANIMATION
% Process frames one by one and write directly to video file
for i = 1:length(sorted_filenames)
    % Read the current image
    img = imread(fullfile(inputDirectory, sorted_filenames{i}));
    
    % Display the image in the figure
    imshow(img);
    
    % Add title showing the current filename (useful for debugging)
    title(sorted_filenames{i}, 'Interpreter', 'none');
    drawnow;
    
    % Capture the current figure as a video frame
    frame = getframe(gcf);
    
    % Write the frame to the video file
    writeVideo(v, frame);
    
    % Show text-based progress indicator
    if i >= progressPoints(nextPoint)
        fprintf('=');
        nextPoint = nextPoint + 1;
    end
    
    % For longer animations, also print numeric progress at intervals
    if mod(i, 50) == 0 || i == length(sorted_filenames)
        fprintf(' %d%% ', round(100*i/length(sorted_filenames)));
    end
end

% Complete the progress bar
fprintf('] Done!\n');

%% STEP 7: FINALIZE AND CLEANUP
% Close the video file
close(v);

% Display completion message with output information
disp(['Animation created successfully as ''', outputFileName, '''']);
disp(['Video contains ', num2str(length(sorted_filenames)), ' frames at ', ...
      num2str(frameRate), ' fps (duration: ', ...
      num2str(length(sorted_filenames)/frameRate, '%.1f'), ' seconds)']);

% Add information about the first and last timestamps
if ~isempty(sorted_filenames)
    % Extract dates from first and last filenames
    firstTokens = regexp(sorted_filenames{1}, pattern, 'tokens');
    lastTokens = regexp(sorted_filenames{end}, pattern, 'tokens');
    
    if ~isempty(firstTokens) && ~isempty(lastTokens)
        firstDateTime = [firstTokens{1}{1}, ' ', strrep(firstTokens{1}{2}, '-', ':')];
        lastDateTime = [lastTokens{1}{1}, ' ', strrep(lastTokens{1}{2}, '-', ':')];
        disp(['Time period: ', firstDateTime, ' to ', lastDateTime]);
    end
end