import numpy as np
from matio import load_from_mat
import pandas as pd
import warnings


mat_file_path = "VertVel_Cornfield_2425_flowVisualized.mat"


print(f"Loading data via mat-io from {mat_file_path}...")

warnings.filterwarnings("ignore", message=".*mat_to_datetime: Ignoring 'fmt' property.*")
mat_data = load_from_mat(mat_file_path, raw_data=False)


# 1. Extract SiteID
siteID = mat_data['siteID'].item().strip()
print(f"siteID: {siteID}\n")

# 2. Extract arrays
velSigned_depthAvg = np.array(mat_data['velSigned_depthAvg']).squeeze()
print(f"velSigned_depthAvg: {velSigned_depthAvg}\n")

velDir_depthAvg = np.array(mat_data['velDir_depthAvg']).squeeze()
print(f"velDir_depthAvg: {velDir_depthAvg}\n")

Eas_depthAvg = np.array(mat_data['Eas_depthAvg']).squeeze()
print(f"Eas_depthAvg: {Eas_depthAvg}\n")

Nor_depthAvg = np.array(mat_data['Nor_depthAvg']).squeeze()
print(f"Nor_depthAvg: {Nor_depthAvg}\n")

# 3. Extract scalars
floodVelMag_depthtimeAvg = float(np.array(mat_data['floodVelMag_depthtimeAvg']).item())
print(f"floodVelMag_depthtimeAvg: {floodVelMag_depthtimeAvg}\n")

ebbVelMag_depthtimeAvg = float(np.array(mat_data['ebbVelMag_depthtimeAvg']).item())
print(f"ebbVelMag_depthtimeAvg: {ebbVelMag_depthtimeAvg}\n")

floodMeanDir = float(np.array(mat_data['floodMeanDir']).item())
print(f"floodMeanDir: {floodMeanDir}\n")

ebbMeanDir = float(np.array(mat_data['ebbMeanDir']).item())
print(f"ebbMeanDir: {ebbMeanDir}\n")

# 4. Extract boolean arrays
floodSign_depthAvg = np.array(mat_data['floodSign_depthAvg']).squeeze().astype(bool)
print(f"floodSign_depthAvg: {floodSign_depthAvg}\n")

ebbSign_depthAvg = np.array(mat_data['ebbSign_depthAvg']).squeeze().astype(bool)
print(f"ebbSign_depthAvg: {ebbSign_depthAvg}\n")


# 5. Extract DMY
raw_dmy = mat_data['DMY'].squeeze()
pandas_dates = pd.to_datetime(raw_dmy)
DMY = pandas_dates.strftime('%m/%d/%Y %H:%M')

print(f"raw_dmy: {raw_dmy}\n")
print(f"DMY: {DMY}\n")
