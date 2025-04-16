# Verdant ADCP Data Analysis

This repository contains scripts and data for the Verdant-SBU tidal project, focusing on the analysis and visualization of ADCP (Acoustic Doppler Current Profiler) data.

## Repository Contents

### Data Files
- `ADCP-Data_Deployment-Period1.xlsx` - ADCP measurements from first deployment period
- `ADCP-Data_Deployment-Period2.xlsx` - ADCP measurements from second deployment period

### MATLAB Scripts
- [`VertVelPlot.m`](VertVelPlot.m) - Main script that processes ADCP data and generates velocity profile plots
- [`VelPlot_Verdant.m`](VelPlot_Verdant.m) - Function that creates various velocity profile visualizations
- [`FolderReadCSV.m`](FolderReadCSV.m) - Function that extracts ADCP data from Excel files
- [`animation.m`](animation.m) - Creates animations from sequential velocity profile plots

### Utilities
- [`DMStoDD.ipynb`](DMStoDD.ipynb) - Python notebook for converting coordinates from DMS to decimal degrees

## Features

The scripts in this repository perform the following operations:
1. Import ADCP data from Excel files
2. Process velocity magnitude and direction data
3. Generate multiple visualization types:
   - 2D vertical velocity profiles
   - Directional compass plots
   - 3D vector visualizations
4. Create animations from time-series plots

## Usage

### Data Processing and Visualization
1. Place ADCP data Excel files in the working directory
2. Run `VertVelPlot.m` to process data and generate profile plots
3. Output visualizations will be saved to folders:
   - `SpeedProfiles/` - 2D velocity magnitude vs. depth plots
   - `DirProfiles/` - Compass plots showing current direction
   - `VelProfiles/` - 3D quiver plots showing U,V components

### Animation Creation
1. Run `animation.m` to create animations from sequential PNG files
2. Adjust frame rate and output settings in the script as needed

## Contributors
- Hossein Seyedzadeh
- Jonathan Craig

## Last Updated
April 2025