#!/usr/bin/env python3
# ebbFloodClassifier.py
# 
# Classifies ADCP tidal flow data into Ebb/Flood/Slack phases using along-channel velocity.
# 
# ===========================================================================================
# LONG ISLAND SOUND TIDAL FLOW CONVENTION (EAST END):
# ===========================================================================================
# At the east end of Long Island Sound:
#   - FLOOD TIDE: Water flows INTO the Sound from Atlantic → WESTWARD flow
#   - EBB TIDE:   Water flows OUT of the Sound to Atlantic → EASTWARD flow
#
# ===========================================================================================
# QUICK START - DIRECT RUN (uses defaults for Verdant Power ADCP data):
# ===========================================================================================
#   python ebbFloodClassifier.py
#
# This will:
#   1. Load: depthAvg_ADCPdata.csv (corrected ADCP data with orientation adjustments)
#   2. Process: Rotate velocities to along-channel/across-channel coordinates
#   3. Classify: Assign Ebb/Flood/Slack phases based on smoothed along-channel velocity
#   4. Output: depthAvg_ADCPdata_labeled.csv with 6 columns:
#      - Date & Time.2 (timestamp)
#      - Eas (East velocity in m/s)
#      - Nor (North velocity in m/s)
#      - u_along (along-channel velocity in m/s, + = westward/flood, - = eastward/ebb)
#      - u_across (across-channel velocity in m/s)
#      - Phase (classification: "EBB", "FLOOD", or "SLACK")
#
# ===========================================================================================
# OPTIONAL VISUALIZATION:
# ===========================================================================================
#   python ebbFloodClassifier.py --plot
#   (Shows time series plot of along-channel velocity with phase thresholds)
#
# ===========================================================================================
# CUSTOM USAGE (override defaults):
# ===========================================================================================
#   python ebbFloodClassifier.py mydata.csv --time-col "DateTime" --east-col "U_east" --north-col "U_north"
#   python ebbFloodClassifier.py data.csv --floodward-bearing 88  # Custom flood direction
#   python ebbFloodClassifier.py data.csv --estimate-axis  # Estimate principal axis from data
#   python ebbFloodClassifier.py data.csv --ebbward-bearing 268  # If you only know ebb direction
#   python ebbFloodClassifier.py data.csv --thr 0.1 --smooth-n 5  # Adjust slack threshold & smoothing
#
# ===========================================================================================
# DEFAULT SETTINGS:
# ===========================================================================================
#   Input CSV:          depthAvg_ADCPdata.csv
#   Time column:        "Date & Time.2"
#   East/North columns: "Eas", "Nor" (auto-converted from mm/s to m/s)
#   Floodward bearing:  252.38°T (westward, into Sound) = 72.38° + 180°
#   Slack threshold:    ±0.05 m/s (velocities between -0.05 and +0.05 m/s are "SLACK")
#   Smoothing window:   3 samples (moving average to reduce noise)
#
# Note: The original PCA analysis gave 72.38°T which pointed EASTWARD (seaward, ebb direction).
#       We add 180° to get the WESTWARD (landward, flood) direction for Long Island Sound.
#
# ===========================================================================================
# CLASSIFICATION METHOD:
# ===========================================================================================
#   1. Rotate (East, North) velocities to (along-channel, across-channel) using flood bearing
#   2. Apply 3-sample centered moving average to along-channel velocity
#   3. Classify based on smoothed velocity (Long Island Sound east end convention):
#      - FLOOD: u_along_smooth > +0.05 m/s  (WESTWARD flow INTO Sound from Atlantic)
#      - EBB:   u_along_smooth < -0.05 m/s  (EASTWARD flow OUT of Sound to Atlantic)
#      - SLACK: |u_along_smooth| ≤ 0.05 m/s (low/transitional flow)

import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def to_mps(series, units_hint):
    """Convert to m/s if the values look like mm/s."""
    if units_hint.lower().endswith("(mm/s)") or units_hint.lower().endswith("mmps"):
        return series / 1000.0
    # If column header didn't have units, guess from magnitude (heuristic).
    if series.abs().median() > 5:  # values like 100~500 look like mm/s
        return series / 1000.0
    return series

def principal_axis_bearing(Eas, Nor):
    """
    Estimate principal flow axis from the data (largest variance direction).
    Return bearing in degrees from North, clockwise.
    """

    # Step 1: Create velocity matrix
    U = np.vstack([Eas, Nor]).T # Shape: [n_samples, 2]

    # Step 2: Center the data (remove mean)
    U = U - U.mean(axis=0, keepdims=True)

    # Step 3: Calculate covariance matrix
    C = (U.T @ U) / (len(U) - 1)

    # Step 4: Find eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eig(C)

    # Step 5: Get the principal direction (largest eigenvalue)
    imax = np.argmax(eigvals)
    v = eigvecs[:, imax]  # [vx, vy] in (E,N)

    # Step 6: Convert to bearing (degrees from North, clockwise)
    # Bearing φ satisfies unit vector (sin φ, cos φ) = (vx, vy)
    phi_deg = math.degrees(math.atan2(v[0], v[1])) % 360.0
    return phi_deg

def rotate_along_across(Eas, Nor, bearing_deg):
    """Rotate (E,N) to (along, across) for bearing_deg (degT) pointing flood-ward."""
    phi = math.radians(bearing_deg)
    u_along  = Eas*np.sin(phi) + Nor*np.cos(phi)
    u_across = Eas*np.cos(phi) - Nor*np.sin(phi)
    return u_along, u_across

def label_phase(u_along_smooth, Uthr=0.05):
    if u_along_smooth >  Uthr: return "FLOOD"
    if u_along_smooth < -Uthr: return "EBB"
    return "SLACK"

def main():
    ap = argparse.ArgumentParser(description="Classify ADCP currents into Ebb/Flood/Slack via axis rotation.")
    ap.add_argument("csv", nargs='?', default="outputs/depthAvg_ADCPdata.csv", help="Input CSV path (default: ./outputs/depthAvg_ADCPdata.csv)")
    ap.add_argument("--time-col", default="Date & Time.2", help="Timestamp column name (default: 'Date & Time.2')")
    ap.add_argument("--east-col", default="Eas", help="Eastward velocity column (default: 'Eas')")
    ap.add_argument("--north-col", default="Nor", help="Northward velocity column (default: 'Nor')")
    ap.add_argument("--floodward-bearing", type=float, default=252.38,
                    help="Channel bearing (deg True) pointing FLOOD-ward (westward into LI Sound) (default: 252.38)")
    ap.add_argument("--ebbward-bearing", type=float, default=None,
                    help="Channel bearing (deg True) pointing EBB-ward (seaward)")
    ap.add_argument("--estimate-axis", action="store_true",
                    help="Estimate axis from data (principal axis); assumed FLOOD-ward")
    ap.add_argument("--thr", type=float, default=0.05, help="Slack threshold in m/s (default 0.05)")
    ap.add_argument("--smooth-n", type=int, default=3, help="Moving-average window (samples; default 3)")
    ap.add_argument("--time-format", default=None,
                    help="Optional strptime format for timestamps (e.g., '%%m/%%d/%%Y %%H:%%M')")
    ap.add_argument("--out-csv", default=None, help="Output CSV path (default: input basename + '_labeled.csv')")
    ap.add_argument("--plot", action="store_true", help="Show a matplotlib plot of along-channel velocity")
    args = ap.parse_args()

    # --- Load data
    df = pd.read_csv(args.csv)
    if args.time_col not in df.columns:
        raise ValueError(f"Time column '{args.time_col}' not found. Available: {list(df.columns)}")
    if args.east_col not in df.columns or args.north_col not in df.columns:
        raise ValueError(f"East/North columns not found. Available: {list(df.columns)}")

    # Parse time
    #args.time_format=True # DEBUGGING on 04, August
    if args.time_format:
        df["Time"] = pd.to_datetime(df[args.time_col], format=args.time_format)
    else:
        df["Time"] = pd.to_datetime(df[args.time_col], errors="coerce")
    
    # Check for parsing errors and provide helpful info
    if df["Time"].isna().any():
        print(f"Warning: {df['Time'].isna().sum()} timestamps could not be parsed.")
        print("Sample problematic timestamps:")
        bad_times = df[df["Time"].isna()][args.time_col].head(5)
        for i, bad_time in enumerate(bad_times):
            print(f"  Row {bad_times.index[i]}: '{bad_time}'")
        
        # If most timestamps are bad, raise error; otherwise continue with valid ones
        if df["Time"].isna().sum() > len(df) * 0.5:  # More than 50% bad
            raise ValueError("More than 50% of timestamps could not be parsed. Provide --time-format if needed.")
        else:
            print(f"Continuing with {(~df['Time'].isna()).sum()} valid timestamps...")
            df = df[~df["Time"].isna()].reset_index(drop=True)

    # Convert to m/s if necessary
    df["Eas"] = to_mps(df[args.east_col].astype(float), args.east_col)
    df["Nor"] = to_mps(df[args.north_col].astype(float), args.north_col)

    # Determine bearing
    bearing_flood = None
    if args.estimate_axis:
        estimated_axis = principal_axis_bearing(df["Eas"].values, df["Nor"].values)
        # For Long Island Sound east end:
        # - Principal axis likely points EASTWARD (seaward, dominant ebb direction ~72°T)
        # - Add 180° to get WESTWARD (landward, flood direction ~252°T)
        bearing_flood = (estimated_axis + 180.0) % 360.0
        bearing_source = f"principal-axis estimate ({estimated_axis:.1f}°T) + 180° (eastward→westward for LI Sound)"
    elif args.ebbward_bearing is not None:
        # User provided ebb direction (eastward), add 180° to get flood (westward)
        bearing_flood = (args.ebbward_bearing + 180.0) % 360.0
        bearing_source = "user (ebb-ward + 180° → flood-ward)"
    else:
        # Use floodward bearing (either user-specified or default 252.38°T westward)
        bearing_flood = args.floodward_bearing % 360.0
        bearing_source = "default (flood-ward westward)" if args.floodward_bearing == 252.38 else "user (flood-ward)"

    # Rotate
    u_along, u_across = rotate_along_across(df["Eas"].values, df["Nor"].values, bearing_flood)
    df["u_along"] = u_along
    df["u_across"] = u_across

    # Smooth
    n = max(1, int(args.smooth_n))
    df = df.sort_values("Time").reset_index(drop=True)
    df["u_along_smooth"] = df["u_along"].rolling(window=n, center=True, min_periods=1).mean()

    # Label
    df["Phase"] = df["u_along_smooth"].apply(lambda x: label_phase(x, Uthr=args.thr))

    # Save
    out_csv = args.out_csv or str(Path(args.csv).with_name(Path(args.csv).stem + "_labeled.csv"))
    keep_cols = [args.time_col, "Eas", "Nor", "u_along", "u_across", "Phase"]
    # filter to only existing
    keep_cols = [c for c in keep_cols if c in df.columns]
    df[keep_cols].to_csv(out_csv, index=False)
    print(f"[OK] Saved labeled file: {out_csv}")
    print(f"[Info] Bearing used (flood-ward): {bearing_flood:.2f}°T  ({bearing_source})")
    print(f"[Info] Slack threshold: ±{args.thr:.3f} m/s   Smooth window: {n} samples")

    # Plot (optional)
    # args.plot = True
    if args.plot:
        plt.figure(figsize=(9,4))
        plt.plot(df["Time"], df["u_along_smooth"], marker="o")
        plt.axhline(args.thr, linestyle="--")
        plt.axhline(-args.thr, linestyle="--")
        plt.title(f"Along-channel velocity (bearing={bearing_flood:.1f}°T flood-ward)")
        plt.xlabel("Time")
        plt.ylabel("u_along (m/s)")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    main()
