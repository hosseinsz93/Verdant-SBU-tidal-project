# Verdant ADCP Data Analysis

This repository contains scripts and data for the Verdant-SBU tidal project, focusing on the analysis and visualization of ADCP (Acoustic Doppler Current Profiler) data and tidal flow patterns.

## Repository Contents

### Data Files
- `ADCP-Data_Deployment-Period1.xlsx` - ADCP measurements from first deployment period
- `ADCP-Data_Deployment-Period2.xlsx` - ADCP measurements from second deployment period

### MATLAB Scripts
- [`VertVelPlot.m`](VertVelPlot.m) - Main script that processes ADCP data and generates velocity profile plots
- [`VelPlot_Verdant.m`](VelPlot_Verdant.m) - Function that creates various velocity profile visualizations
- [`FolderReadCSV.m`](FolderReadCSV.m) - Function that extracts ADCP data from Excel files
- [`animation.m`](animation.m) - Creates animations from sequential velocity profile plots
- [`AnalyzeTidalFlow.m`](AnalyzeTidalFlow.m) - Analyzes tidal flow patterns from ADCP data
- [`ExtractCvelAvg.m`](ExtractCvelAvg.m) - Extracts depth-averaged velocity data
- [`ExtractCvelAvgWithTimestamps.m`](ExtractCvelAvgWithTimestamps.m) - Extracts velocity data with timestamps
- [`identifyTidalFlow.m`](identifyTidalFlow.m) - Classifies ADCP measurements into flood and ebb tide periods
- [`TimeAveragedProfiles.m`](TimeAveragedProfiles.m) - Generates time-averaged velocity profiles

### Utilities
- [`DMStoDD.ipynb`](DMStoDD.ipynb) - Python notebook for converting coordinates from DMS to decimal degrees

## Features

The scripts in this repository perform the following operations:

### Velocity Profile Analysis
- Import ADCP data from Excel files
- Process velocity magnitude and direction data
- Generate multiple visualization types:
  - 2D vertical velocity profiles
  - Directional compass plots
  - 3D vector visualizations
- Create animations from time-series plots

### Tidal Flow Analysis
- Identify flood and ebb tides based on current direction
- Calculate tidal flow statistics (magnitude, direction, duration)
- Visualize tidal patterns with:
  - Current rose diagrams
  - Time series of signed velocities
  - Statistical comparisons of flood vs. ebb tides
- Generate analytical models of velocity profiles using:
  - Power law approximations
  - Logarithmic profiles
  - Polynomial fits

## Usage

### Basic Velocity Profile Processing
1. Place ADCP data Excel files in the working directory
2. Run `VertVelPlot.m` to process data and generate profile plots
3. Output visualizations will be saved to folders:
   - `SpeedProfiles/` - 2D velocity magnitude vs. depth plots
   - `DirProfiles/` - Compass plots showing current direction
   - `VelProfiles/` - 3D quiver plots showing U,V components

### Tidal Flow Analysis
1. Configure the working folder in `AnalyzeTidalFlow.m`
2. Run the script to perform comprehensive tidal analysis
3. Examine the generated plots and statistical outputs
4. Results are saved to:
   - MAT file with all calculated parameters
   - Text file with detailed statistics
   - Multiple visualization figures

### Animation Creation
1. Run `animation.m` to create animations from sequential PNG files
2. Adjust frame rate and output settings in the script as needed

## Contributors
- Hossein Seyedzadeh
- Jonathan Craig

## Last Updated
April 2025