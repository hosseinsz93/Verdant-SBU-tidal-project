from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import dates as mdates
from matplotlib.ticker import FixedLocator
from matio import load_from_mat
import warnings
import pandas as pd


def flowVisualized(mainOutputDir, DMY, velSigned_depthAvg, velDir_depthAvg, floodSign_depthAvg, ebbSign_depthAvg,
                   Eas_depthAvg, Nor_depthAvg, floodVelMag_depthtimeAvg, ebbVelMag_depthtimeAvg, floodMeanDir,
                   ebbMeanDir, siteID):

    # Check whether output folder exists
    if not Path(mainOutputDir).is_dir():
        raise FileNotFoundError("The flow visualization function does not receive a valid main output directory.")

    # Sets default font to Times New Roman
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['savefig.dpi'] = 600

    generate_figure1(mainOutputDir, DMY, floodVelMag_depthtimeAvg, ebbVelMag_depthtimeAvg, floodMeanDir, ebbMeanDir,
                     siteID)
    generate_figure2(mainOutputDir, DMY, velSigned_depthAvg, floodSign_depthAvg, ebbSign_depthAvg, siteID)
    generate_figure3(mainOutputDir, DMY, Eas_depthAvg, Nor_depthAvg, siteID)
    generate_figure4(mainOutputDir, DMY, velSigned_depthAvg, floodSign_depthAvg, ebbSign_depthAvg, siteID)
    generate_figure5(mainOutputDir, DMY, velSigned_depthAvg, velDir_depthAvg, siteID)
    generate_figure6(mainOutputDir, DMY, velSigned_depthAvg, floodSign_depthAvg, ebbSign_depthAvg, siteID)


def calculate_midnight_ticks(DMY, freq='7D'):
    """
    Generates tick dates fixed at 00:00 (midnight) from a list/numpy array of date strings (DMY).
    Ensures the first tick starts beore or at the date start, and the last tick extends past the data end.
    """
    dmy_datetime = pd.to_datetime(DMY)

    start_date = dmy_datetime.min().floor('D')
    end_date = dmy_datetime.max().ceil('D')

    tick_dates = pd.date_range(start=start_date, end=end_date, freq=freq)

    if tick_dates[-1] < end_date:
        tick_dates = tick_dates.append(pd.DatetimeIndex(tick_dates[-1]) + pd.Timedelta(freq))

    return tick_dates


# Figure 1: Summary of Info about Tidal Flow (_FlowSummary.png)
def generate_figure1(mainOutputDir, DMY, floodVelMag_depthtimeAvg, ebbVelMag_depthtimeAvg, floodMeanDir, ebbMeanDir,
                     siteID):
    fig_title = f"{siteID} ({DMY[0]} to {DMY[-1]})"
    fig = plt.figure(num=fig_title, figsize=(1000.0, 450.0, 'px'))

    # Axis 1: Bar Chart
    ax_bar = fig.add_axes((0.1, 0.15, 0.35, 0.7))
    ax_bar.bar(['Flood', 'Ebb'], [floodVelMag_depthtimeAvg, ebbVelMag_depthtimeAvg], edgecolor='k', width=0.6)
    ax_bar.set_ylabel('Average Magnitude (m/s)')
    ax_bar.set_title(fig_title, fontweight='bold', fontsize=11)
    ax_bar.set_ylim(0, 0.6)
    ax_bar.set_axisbelow(True)
    ax_bar.grid(True, linestyle='-', color='#E0E0E0', alpha=0.7)


    # Axis 2: Polar Plot
    ax_polar = fig.add_axes((0.55, 0.1, 0.4, 0.8), polar=True)
    ax_polar.set_theta_zero_location('N')
    ax_polar.set_theta_direction(-1)

    floodRad = np.radians(floodMeanDir)
    ebbRad = np.radians(ebbMeanDir)

    line1, = ax_polar.plot([0, floodRad], [0, floodVelMag_depthtimeAvg], linewidth=2, color=[0.2, 0.4, 0.9],
                           label='Flood')
    line2, = ax_polar.plot([0, ebbRad], [0, ebbVelMag_depthtimeAvg], linewidth=2, color=[0.9, 0.2, 0.2],
                           label='Ebb')

    ax_polar.scatter(floodRad, floodVelMag_depthtimeAvg, s=70, color=[0.2, 0.4, 0.9], zorder=3)
    ax_polar.scatter(ebbRad, ebbVelMag_depthtimeAvg, s=70, color=[0.9, 0.2, 0.2])

    ax_polar.set_thetagrids(np.arange(0, 360, 45), ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
    ax_polar.set_rlim(0, 0.6)
    ax_polar.set_rticks([0.2, 0.4, 0.6])
    ax_polar.set_rlabel_position(0)
    ax_polar.grid(True, linestyle='-', color='#E0E0E0', alpha=0.7)
    ax_polar.text(floodRad-0.06, floodVelMag_depthtimeAvg-0.07, f'{floodMeanDir:.1f}\u00b0 T', color=[0.2, 0.4, 0.9],
                  weight='bold', ha='right')
    ax_polar.text(ebbRad+0.03, ebbVelMag_depthtimeAvg-0.07, f"{ebbMeanDir:.1f}\u00b0 T", color=[0.9, 0.2, 0.2],
                  weight='bold', ha='left')
    ax_polar.set_title(f'{siteID} ({DMY[0]} to {DMY[-1]}):\nDominant Current Directions (True North)', pad=20,
                       fontweight='bold')
    ax_polar.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07), ncol=1, frameon=True, edgecolor='k')

    plt.savefig(Path(mainOutputDir) / '_FlowSummary.png', bbox_inches='tight')
    plt.close(fig)
    print("fig1 (_FlowSummary.png) generated.")


# Figure 2: Time Series (_TimeSeries.png)
def generate_figure2(mainOutputDir, DMY, velSigned_depthAvg, floodSign_depthAvg, ebbSign_depthAvg, siteID):
    dmy_arr = pd.to_datetime(np.array(DMY)).to_numpy()
    vel_arr = np.array(velSigned_depthAvg)

    fig_title = f'{siteID} ({DMY[0]} to {DMY[-1]}):\nTidal Flow Analysis: Flood (+) and Ebb (-) Tides'
    fig = plt.figure(num=fig_title, figsize=(1200.0, 600.0, 'px'))

    ax = fig.add_subplot(111)

    ax.plot(dmy_arr, vel_arr, 'k', linewidth=0.5, label='Velocity')

    ax.plot(dmy_arr[floodSign_depthAvg], vel_arr[floodSign_depthAvg], 'b.', markersize=4, label='Flood Tide')
    ax.plot(dmy_arr[ebbSign_depthAvg], vel_arr[ebbSign_depthAvg], 'r.', markersize=4, label='Ebb Tide')

    ax.set_xlabel('Time')
    ax.set_ylabel('Velocity (m/s)')
    ax.set_title(fig_title, pad=15, fontweight='bold')

    tick_dates = calculate_midnight_ticks(DMY, freq='7D')

    ax.set_xlim(tick_dates[0], tick_dates[-1])

    ax.xaxis.set_major_locator(FixedLocator(mdates.date2num(tick_dates)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))

    ax.minorticks_off()
    plt.setp(ax.get_xticklabels(), rotation=0, ha='center')

    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5])

    ax.legend(loc='upper right', frameon=True, edgecolor='k')
    ax.grid(True, linestyle='-', color='#E0E0E0', alpha=0.7)

    plt.savefig(Path(mainOutputDir) / '_TimeSeries.png', bbox_inches='tight')
    plt.close(fig)
    print("fig2 (_TimeSeries.png) generated.")


# Figure 3: Velocity Components (_Componenets.png)
def generate_figure3(mainOutputDir, DMY, Eas_depthAvg, Nor_depthAvg, siteID):
    dmy_arr = pd.to_datetime(np.array(DMY)).to_numpy()
    east_arr = np.array(Eas_depthAvg)
    north_arr = np.array(Nor_depthAvg)

    fig_title = f'{siteID} ({DMY[0]} to {DMY[-1]}): Velocity Components'
    fig = plt.figure(num=fig_title, figsize=(1200.0, 400.0, 'px'))

    ax = fig.add_subplot(111)

    ax.plot(dmy_arr, east_arr, 'r-', linewidth=1.2, label='Eastward')
    ax.plot(dmy_arr, north_arr, 'g-', linewidth=1.2, label='Northward')

    ax.set_xlabel('Time')
    ax.set_ylabel('Velocity (m/s)')
    ax.set_title(fig_title, pad=15, fontweight='bold')

    tick_dates = calculate_midnight_ticks(DMY, freq='7D')

    ax.set_xlim(tick_dates[0], tick_dates[-1])

    ax.xaxis.set_major_locator(FixedLocator(mdates.date2num(tick_dates)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))

    ax.minorticks_off()

    ax.legend(loc='upper right', frameon=True, edgecolor='k')
    ax.grid(True, linestyle='-', color='#E0E0E0', alpha=0.7)

    plt.savefig(Path(mainOutputDir) / '_Components.png', bbox_inches='tight')
    plt.close(fig)
    print("fig3 (_Components) generated.")


# Figure 4: Statistical Analysis (_Statistics.png)
def generate_figure4(mainOutputDir, DMY, velSigned_depthAvg, floodSign_depthAvg, ebbSign_depthAvg, siteID):
    vel_arr = np.array(velSigned_depthAvg)
    flood_stats = vel_arr[floodSign_depthAvg]
    ebb_stats = vel_arr[ebbSign_depthAvg]

    fig_title = f'{siteID} ({DMY[0]} to {DMY[-1]}): Statistical Analysis'
    fig = plt.figure(num=fig_title, figsize=(800.0, 600.0, 'px'))

    # Axis 1: Velocity Distribution


    # Axis 2: Flood/Ebb Box Plots



# Figure 5: Current Rose
def generate_figure5(mainOutputDir, DMY, velSigned_depthAvg, velDir_depthAvg, SiteID):
    pass


# Figure 6: Tidal Cycles
def generate_figure6(mainOutputDir, DMY, velSigned_depthAvg, floodSign_depthAvg, ebbSign_depthAvg, SiteID):
    pass


if __name__ == '__main__':
    mat_file_path = "VertVel_Cornfield_2425_flowVisualized.mat"
    output_directory = "./test_output"
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    print(f"Loading data via mat-io from {mat_file_path}...")

    warnings.filterwarnings("ignore", message=".*mat_to_datetime: Ignoring 'fmt' property.*")
    mat_data = load_from_mat(mat_file_path, raw_data=False)

    siteID = mat_data['siteID'].item().strip()
    velSigned_depthAvg = np.array(mat_data['velSigned_depthAvg']).squeeze()
    velDir_depthAvg = np.array(mat_data['velDir_depthAvg']).squeeze()
    Eas_depthAvg = np.array(mat_data['Eas_depthAvg']).squeeze()
    Nor_depthAvg = np.array(mat_data['Nor_depthAvg']).squeeze()
    floodVelMag_depthtimeAvg = float(np.array(mat_data['floodVelMag_depthtimeAvg']).item())
    ebbVelMag_depthtimeAvg = float(np.array(mat_data['ebbVelMag_depthtimeAvg']).item())
    floodMeanDir = float(np.array(mat_data['floodMeanDir']).item())
    ebbMeanDir = float(np.array(mat_data['ebbMeanDir']).item())
    floodSign_depthAvg = np.array(mat_data['floodSign_depthAvg']).squeeze().astype(bool)
    ebbSign_depthAvg = np.array(mat_data['ebbSign_depthAvg']).squeeze().astype(bool)
    raw_dmy = mat_data['DMY'].squeeze()
    pandas_dates = pd.to_datetime(raw_dmy)
    DMY = pandas_dates.strftime('%m/%d/%Y %H:%M')

    print(f"Data from {mat_file_path} loaded successfully.\n")



    flowVisualized(
        mainOutputDir=output_directory,
        DMY=DMY,
        velSigned_depthAvg=velSigned_depthAvg,
        velDir_depthAvg=velDir_depthAvg,
        floodSign_depthAvg=floodSign_depthAvg,
        ebbSign_depthAvg=ebbSign_depthAvg,
        Eas_depthAvg=Eas_depthAvg,
        Nor_depthAvg=Nor_depthAvg,
        floodVelMag_depthtimeAvg=floodVelMag_depthtimeAvg,
        ebbVelMag_depthtimeAvg=ebbVelMag_depthtimeAvg,
        floodMeanDir=floodMeanDir,
        ebbMeanDir=ebbMeanDir,
        siteID=siteID
    )





