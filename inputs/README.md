# Inputs Folder

This folder contains the formatted files from CO-OPS for reading by FolderReadDF.ipynb

## Format

Each CSV file has the following name: {Site ID}+{Depth bin}+{Year-month-day}+.csv
- The CO-OPS website gives the site ID: e.g., LIS1001
- The depth bin is given in units of meters and centimeters, with two digit places for each unit. CO-OPS lists depths in both feet and meters
- Dates are written as year, month, and day of the deployment
- Timeframes greater than 1 month entail files with the same site ID and depth bin, but with distinct deployment dates
- Underscores (_) separate the site ID from the depth bin, and dashes (-) separate the depth bin and the date as well as the date units from each other
