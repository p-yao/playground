# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 09:36:48 2020

Calculate freezing level at 3 locations in Yuba watershed

Save data for entire 2020-2021 water year
Include WWRF-NRT, CNRFC, and radar obs

Data:
    - WWRF-NRT GFS and ECMWF
    - CNRFC freeze obs and fcst
    - FMCW radar data
    - QPE, for reference

Details:
    - 6 hour time steps
    - currently for 2019-2020, 2020-2021 water years
    - radar data calculated on multiple "windows" centered on the 6 hour time steps
    
Window calculation:
    - FMCW data is on 10 min time steps, so there are up to 6*6=36 values at each step
    - Start with the 2 values centered on each time, then choose several window sizes
    - How about: 2, 4, 10, 20, 30, 36
    
Selected: include 4 window sizes:
    2 (10 minutes)
    6 (50 minutes)
    18 (170 minutes)
    36 (350 minutes)

@author: Peter Yao
"""

from pathlib import Path
from datetime import datetime, timedelta
import argparse
import numpy as np
import pandas as pd
import netCDF4
# import freeze_functions as ffun

#--------------------------------------------------------------------------
# User flags

# Parser controls in the input options from command line.
parser = argparse.ArgumentParser(description='Input options for WRF freezing level script')
parser.add_argument('-s', '--station', required=True, help="Station, ex. 'OVL'")
parser.add_argument('-y', '--year', required=True, help="Starting year")
parser.add_argument('-w', '--wrf', required=True, help="Include WRF data, y or n")
args = parser.parse_args()

start_year = int(args.year) # format: yyyy
end_year = start_year + 1

# Valid dates (entire water year)
start_date = datetime(start_year,11,1)
end_date = datetime(end_year,4,30)

# start_date = datetime(2020,11,1)
# end_date = datetime(2021,4,30)
num_days = (end_date - start_date).days + 1
valid_dates = [start_date + timedelta(days=x) for x in range(num_days)]

# Water year
year = str(start_year) + "-" + str(end_year)
# year = "2020-2021"

# Location to plot (for freezing level plots)
# station = "OVL"
station = args.station

# WRF flag
wrf_flag = args.wrf # 'y' = include, 'n' = skip

# Output file name
out_name = "Freezing_Level_Dataset_" + station + "_" + year + ".nc"

print(station)
print(year)
print(wrf_flag)


#--------------------------------------------------------------------------
# Station set up

station_names = {"NBB":"New Bullards Bar Dam", "OVL":"Oroville", "CFF":"Colfax"}
radar_types = {"NBB":"MRR", "OVL":"FMCW", "CFF":"FMCW"}
station_coords = {"NBB":[39.3964, -121.1438], "OVL":[39.5318, -121.4876], "CFF":[39.0798, -120.9379]}
elevations = {"NBB":632, "OVL":114, "CFF":644} # in meters


#--------------------------------------------------------------------------
# Path to script name
script = "/home/peyao/Python_Scripts/wrf_nrt_freeze_window_selected.py"

# WWRF-NRT data
wrfdir = Path("/data/downloaded/Forecasts/Model_Verification/WWRF-NRT/Z0C")

# QPE precipitation
qpedir = Path("/data/downloaded/Observations/CNRFC/QPE")

# CNRFC Freezing level data
# frz_obs_dir = Path("/home/peyao/CNRFC/FRZ_OBS")
# frz_fct_dir = Path("/home/peyao/CNRFC/FRZ_FCST")
frz_obs_dir = Path("/data/downloaded/Observations/CNRFC/FRZ_OBS")
frz_fct_dir = Path("/data/downloaded/Observations/CNRFC/FRZ_FCST")

# FMCW data
# fmcw_dir = Path("/home/peyao/HMT")
fmcw_dir = Path("/data/downloaded/SCRATCH/peyao_scratch/HMT/BrightBand") # Station is lowercase

# MRR data
mrr_dir = Path("/data/CW3E_data/CW3E_MRR_Archive")

# Output directory - nc file
ncdir = Path("/data/downloaded/SCRATCH/peyao_scratch/Freezing_Level_Dataset") / station
if not ncdir.is_dir():
    ncdir.mkdir(parents=True)

# Define number of lead/valid days - WRF in 3 hr, QPE in 6 hr steps
lead_days = 4
wrf_tsteps = 4 # Changed back to 4 - we only need every other time step
qpe_tsteps = 4

# Radar type
radar_type = radar_types.get(station)


# ---------------------------------------------------------------------------
# Import grid files
print("Importing grid files")

# Read in old (pre-2020) QPE grid
qpe_old = Path("/home/peyao/Grids/lat_lon_qpe_regrid.nc")

with netCDF4.Dataset(qpe_old, mode='r') as infile:
    
    # Lat/lon - subset: for faster calculations during 20 top events
    qpe_lat_old = infile.variables["lat"][:]
    qpe_lon_old = infile.variables["lon"][:]

# Read in QPE watershed file
qpe_shape = Path("/home/peyao/Grids/qpe_lat_lon_post2020.nc")

with netCDF4.Dataset(qpe_shape, mode='r') as infile:
    
    # Lat/lon - subset: for faster calculations during 20 top events
    qpe_lat = infile.variables["lat"][:]
    qpe_lon = infile.variables["lon"][:]

# ---------------------------------------------------------------------------
# Read in QPF watershed file
qpf_shape = Path("/data/downloaded/Observations/CNRFC/QPF/lat_lon_qpf.nc")

with netCDF4.Dataset(qpf_shape, mode='r') as infile:
    
    # Lat/lon
    qpf_lat = infile.variables["lat"][:]
    qpf_lon = infile.variables["lon"][:]

# ---------------------------------------------------------------------------
# Read in lat/lon (WRF)
# Domain 1 (9km)
grid_file1 = Path("/data/downloaded/WWRF-NRT") / "2020-2021" / "WRF_Output" / "2020111500" / "gfs" / "wrfout_d01_2020-11-15_00_00_00.nc"
# Domain 2 (3km)
grid_file2 = Path("/data/downloaded/WWRF-NRT") / "2020-2021" / "WRF_Output" / "2020111500" / "gfs" / "wrfout_d02_2020-11-15_00_00_00.nc"

# Read in 9km grid
with netCDF4.Dataset(grid_file1, mode='r') as infile:
    
    # Lat/lon
    xlat1 = infile.variables["XLAT"][0]
    xlon1 = infile.variables["XLONG"][0]

# Read in 3km grid
with netCDF4.Dataset(grid_file2, mode='r') as infile:
    
    # Lat/lon
    xlat2 = infile.variables["XLAT"][0]
    xlon2 = infile.variables["XLONG"][0]

# # Grid size
# lat_size1 = np.size(xlat1,0)
# lon_size1 = np.size(xlat1,1)
# lat_size2 = np.size(xlat2,0)
# lon_size2 = np.size(xlat2,1)


# ---------------------------------------------------------------------------
# Get index closest to each station
# Get lat/lon
lat = station_coords.get(station)[0]
lon = station_coords.get(station)[1]

# Find index of closest lat/lon in each grid
# WRF grid (domain 1)
dist = np.sqrt(np.square(xlat1 - lat) + np.square(xlon1 - lon))
wrf_idxs1 = np.unravel_index(dist.argmin(), dist.shape)

# WRF grid (domain 2)
dist = np.sqrt(np.square(xlat2 - lat) + np.square(xlon2 - lon))
wrf_idxs2 = np.unravel_index(dist.argmin(), dist.shape)

# QPE grid (2020 WY -)
dist = np.sqrt(np.square(qpe_lat - lat) + np.square(qpe_lon - lon))
qpe_idxs = np.unravel_index(dist.argmin(), dist.shape)

# Old QPE grid (pre 2020 WY)
dist = np.sqrt(np.square(qpe_lat_old - lat) + np.square(qpe_lon_old - lon))
qpe_idxs_old = np.unravel_index(dist.argmin(), dist.shape)

# QPF grid
dist = np.sqrt(np.square(qpf_lat - lat) + np.square(qpf_lon - lon))
qpf_idxs = np.unravel_index(dist.argmin(), dist.shape)


# ---------------------------------------------------------------------------
# Create list of valid dates

# Number of days
# num_days = len(valid_dates)

# Number of WRF time steps
# wrf_steps = wrf_tsteps*num_days # Changed back to 4 WRF steps per day (6 hr)
qpe_steps = qpe_tsteps*num_days
# fmcw_steps = 144*num_days
# mrr_steps = 24*num_days

# Time array start times (minutes after 10/1/2020)
# wrf_start = 0
# qpe_start = 0
# fmcw_start = 5
# mrr_start = 30

# wrf_start = start_date
qpe_start = start_date
# fmcw_start = start_date + timedelta(minutes=5)
# mrr_start = start_date + timedelta(minutes=30)

# Create time arrays (for output netCDF array)
# wrf_time = [wrf_start + 3*60*x for x in range(wrf_steps)]
# qpe_time = [qpe_start + 6*60*x for x in range(qpe_steps)]
# fmcw_time = [fmcw_start + 10*x for x in range(fmcw_steps)]
# mrr_time = [mrr_start + 60*x for x in range(mrr_steps)]

# wrf_time1 = [wrf_start + timedelta(minutes=3*60*x) for x in range(wrf_steps)]
qpe_time1 = [qpe_start + timedelta(minutes=6*60*x) for x in range(qpe_steps)]
# fmcw_time1 = [fmcw_start + timedelta(minutes=10*x) for x in range(fmcw_steps)]
# mrr_time1 = [mrr_start + timedelta(minutes=60*x) for x in range(mrr_steps)]

# wrf_time = [(t - datetime(1970, 1, 1)).total_seconds() for t in wrf_time1]
qpe_time = [(t - datetime(1970, 1, 1)).total_seconds() for t in qpe_time1]
# fmcw_time = [(t - datetime(1970, 1, 1)).total_seconds() for t in fmcw_time1]
# mrr_time = [(t - datetime(1970, 1, 1)).total_seconds() for t in mrr_time1]


# Create arrays to check for missing data
# gfs_missing1 = np.zeros((lead_days, qpe_steps))
# gfs_missing2 = np.zeros((lead_days, qpe_steps))

# ecm_missing1 = np.zeros((lead_days, qpe_steps))
# ecm_missing2 = np.zeros((lead_days, qpe_steps))

cnrfc_obs_missing = np.zeros(qpe_steps)
cnrfc_fct_missing = np.zeros((lead_days, qpe_steps))
cnrfc_qpe_missing = np.zeros(qpe_steps)

# fmcw_missing = np.zeros(qpe_steps)
mrr_missing = np.zeros(qpe_steps)


# ---------------------------------------------------------------------------
# Main loop
# Check WRF flag
if wrf_flag == "y":

    print("Main loop: calculating WRF output")
    
    # Create arrays to check for missing data
    gfs_missing1 = np.zeros((lead_days, qpe_steps))
    gfs_missing2 = np.zeros((lead_days, qpe_steps))
    
    ecm_missing1 = np.zeros((lead_days, qpe_steps))
    ecm_missing2 = np.zeros((lead_days, qpe_steps))
    
    # Create output arrays
    # Lead days: 4 day to 0 day
    # wrf_tsteps: number of total time steps in valid period
    
    # gfs
    freeze_gfs1 = np.full([lead_days, qpe_steps], np.nan) # 9km
    freeze_gfs2 = np.full([lead_days, qpe_steps], np.nan) # 3km
    
    # ecmwf
    freeze_ecm1 = np.full([lead_days, qpe_steps], np.nan) # 9km
    freeze_ecm2 = np.full([lead_days, qpe_steps], np.nan) # 3km
    
    # Counter
    count = 0
    
    # Loop through each valid date in newly created list
    for v in range(num_days):
        
        # Get valid date
        valid_date = valid_dates[v]
        
        # Convert date to string
        valid_str = datetime.strftime(valid_date, "%Y-%m-%d")
        
        # Import data
        print("WWRF-NRT: Valid " + valid_str)
        
        # List of init dates
        init_dates = [valid_date - timedelta(days=x) for x in range(lead_days,0,-1)]
        
        # Starting index depends on init date (4 lead days)
        start_idx = [96, 72, 48, 24]
        
        # Loop through each init date
        for i, init_date in enumerate(init_dates):
            
            # Convert date to string
            init_str = datetime.strftime(init_date, "%Y%m%d00")
            
            # Get list of forecast lead times
            f_times = [start_idx[i] + 6*x for x in range(wrf_tsteps)] # Changed back to 6 hrs
            
            # Loop through each 6 hr time step
            for f, f_time in enumerate(f_times):
            
                # Get path to WRF file
                # GFS 9km
                gfs_file1 = wrfdir / "gfs" / "9km" / year / init_str / ("WWRF-NRT_Z0C_d01_" + init_str + "_F" + str(f_time).zfill(3) + ".nc")
                # GFS 3km
                gfs_file2 = wrfdir / "gfs" / "3km" / year / init_str / ("WWRF-NRT_Z0C_d02_" + init_str + "_F" + str(f_time).zfill(3) + ".nc")
                
                # ECMWF 9km
                ecm_file1 = wrfdir / "ecmwf" / "9km" / year / init_str / ("WWRF-NRT_Z0C_d01_" + init_str + "_F" + str(f_time).zfill(3) + ".nc")
                # ECMWF 3km
                ecm_file2 = wrfdir / "ecmwf" / "3km" / year / init_str / ("WWRF-NRT_Z0C_d02_" + init_str + "_F" + str(f_time).zfill(3) + ".nc")
                
                # GFS
                # Check if files exist
                if gfs_file1.is_file():
                
                    # Read in data
                    with netCDF4.Dataset(gfs_file1, mode='r') as infile:
                        
                        # Total precipitation
                        gfs_temp1 = infile.variables["Z0C"][wrf_idxs1[0], wrf_idxs1[1]]
                        
                    # Save 6hr time step to output array
                    freeze_gfs1[i, count+f] = gfs_temp1
                    
                    # Mark data as present
                    gfs_missing1[i, count+f] = 1
                    
                else:
                    print(str(gfs_file1))
                
                if gfs_file2.is_file():
                    with netCDF4.Dataset(gfs_file2, mode='r') as infile:
                        
                        # Total precipitation
                        gfs_temp2 = infile.variables["Z0C"][wrf_idxs2[0], wrf_idxs2[1]]
                        
                    # Save 6hr time step to output array
                    freeze_gfs2[i, count+f] = gfs_temp2
                    
                    # Mark data as present
                    gfs_missing2[i, count+f] = 1
                    
                else:
                    print(str(gfs_file2))
                
                
                # ECMWF
                # Check if files exist
                if ecm_file1.is_file():
                
                    # Read in data
                    with netCDF4.Dataset(ecm_file1, mode='r') as infile:
                        
                        # Total precipitation
                        ecm_temp1 = infile.variables["Z0C"][wrf_idxs1[0], wrf_idxs1[1]]
                        
                    # Save 6hr time step to output array
                    freeze_ecm1[i, count+f] = ecm_temp1
                    
                    # Mark data as present
                    ecm_missing1[i, count+f] = 1
                    
                else:
                    print(str(ecm_file1))
                
                if ecm_file2.is_file():
                    with netCDF4.Dataset(ecm_file2, mode='r') as infile:
                        
                        # Total precipitation
                        ecm_temp2 = infile.variables["Z0C"][wrf_idxs2[0], wrf_idxs2[1]]
                        
                    # Save 6hr time step to output array
                    freeze_ecm2[i, count+f] = ecm_temp2
                    
                    # Mark data as present
                    ecm_missing2[i, count+f] = 1
                    
                else:
                    print(str(ecm_file2))
                
                
                # Remove any fill values (-30000)
                # wrf_temp[wrf_temp < 0] = np.nan
        
                # # Save 6hr time step to output array
                # freeze_nrt1[i, count+f] = wrf_temp1
                # freeze_nrt2[i, count+f] = wrf_temp2
                
        # Increment count
        count += wrf_tsteps # now 4

else:
    print("Skipping WRF output")


# ---------------------------------------------------------------------------
# Calculate freeze obs
print("Calculating freeze obs")

# # Create array to store output data (this is used to calculate freezing error)
# freeze_obs_fct = np.full([lead_days, qpe_steps], np.nan)

# Array for all 6 hr values (4 values per day)
freeze_obs_6hr = np.full([qpe_steps], np.nan)

# Counter
count = 0

# Loop through each valid date in newly created list
for v in range(num_days):
    
    # Get valid date
    valid_date = valid_dates[v]
    
    # Convert date to string
    valid_str = datetime.strftime(valid_date, "%Y-%m-%d")
    
    # Import data
    print("Freeze Obs: Valid " + valid_str)
    
    # Read in QPE watershed file
    # qrow1 = qpe_set[0]
    # qrow2 = qpe_set[1]
    # qcol1 = qpe_set[2]
    # qcol2 = qpe_set[3]

    # Create list of 6 hr time steps
    time_steps_obs = [datetime.strftime(valid_date + timedelta(hours=6*x), "%Y%m%d_%H") for x in range(qpe_tsteps)]
    hour_files_obs = ["oz." + t + "00.nc" for t in time_steps_obs]
    
    # Read in QPE
    for h, hour_file in enumerate(hour_files_obs):
        
        # Full path
        obs_file = frz_obs_dir / hour_file
        
        # Check if file is present
        if obs_file.is_file():
        
            # Read in data
            with netCDF4.Dataset(obs_file, mode='r') as infile:
                
                # Freezing level
                # try:
                #     obs_temp = infile.variables["FzLevel_SFC"][0, qpf_idxs[0], qpf_idxs[1]]
                # except:
                #     obs_temp = infile.variables["htfl"][0, qpf_idxs[0], qpf_idxs[1]]
                    
                # Old grid (pre-2020)
                if start_year < 2020:
                    obs_temp = infile.variables["htfl"][0, 0, qpe_idxs_old[0], qpe_idxs_old[1]]
                else:
                    obs_temp = infile.variables["FzLevel_SFC"][0, qpf_idxs[0], qpf_idxs[1]]
                
            # # Remove any fill values (-99999)
            # obs_temp[obs_temp < 0] = np.nan
            if obs_temp <= 0:
                obs_temp = np.nan
            
            # # Get station idx
            # idx = qpe_idxs
            
            # # Get array subset
            # stat_temp = obs_temp[idx[0], idx[1]]
        
            # Store each time step to temp array
            # qpe_temp[h] = np.nansum(obs_temp[qpe_grid == 1])/num_points2
            freeze_obs_6hr[count+h] = obs_temp/3.28084 # ft to m
            
            # # Obs arrays for calculating error
            # freeze_obs_fct[:, count+h] = obs_temp/3.28084
            
            # Mark file as present
            cnrfc_obs_missing[count+h] = 1
        
        else:
            print(str(obs_file))
        
    # Increment count
    count += qpe_tsteps # 4
    
    
# ---------------------------------------------------------------------------
# Calculate freeze forecast
print("Calculating freeze forecast")

# Add one more day to account for 12Z to 12Z offset
# Then at the end remove the first/last 2 time steps
valid_dates_fct = valid_dates + [valid_dates[-1] + timedelta(days=1)]

# # Create array to store output data - means
# freeze_fct = np.full([lead_days, qpe_steps], np.nan)

# Array for all 6 hr values (4 values per day)
# This is temp array that has extra 2 tsteps on each end
freeze_fct_temp = np.full([lead_days, qpe_steps+4], np.nan)

# Counter
count = 0

# Loop through each valid date in newly created list
# for v, valid_date in enumerate(valid_dates):
for v in range(num_days+1):
    
    # Get valid date
    valid_date = valid_dates_fct[v]
    
    # Convert date to string
    valid_str = datetime.strftime(valid_date, "%Y-%m-%d")
    
    # Import data
    print("Freeze Forecast: Valid " + valid_str)
    
    # List of init dates
    # qpf_lead = reversed([36 + 24*x for x in range(lead_days)])
    # init_dates = [valid_date - timedelta(hours=x) for x in qpf_lead]
    init_dates = [valid_date - timedelta(days=x) for x in range(lead_days,0,-1)]
    
    # Hour starting index for each init date
    # start_idx = [3,7,11,15,19]
    # start_idx = [15,11,7,3] # QPF is in 6 hr time steps
    start_idx = [12,8,4,0]
    
    # Loop through each init date
    for i, init_date in enumerate(init_dates):
        
        # Get start index for current init date
        start = start_idx[i]
        
        # Convert date to string
        init_str = datetime.strftime(init_date, "%Y%m%d")
        
        # QPF full path to file
        qpf_file = frz_fct_dir / ("fzlevel." + init_str + "_1200.nc")
        
        # Check if file is present
        if qpf_file.is_file():
        
            # Read in data
            with netCDF4.Dataset(qpf_file, mode='r') as infile:
                
                # Total precipitation
                fct_temp = infile.variables["FzLevel_SFC"][:, qpf_idxs[0], qpf_idxs[1]]
                
            # Remove any fill values (-30000)
            fct_temp[fct_temp <= 0] = np.nan
            
            # # Get station idx
            # idx = qpf_idxs
            
            # Get array slice
            stat_temp = fct_temp[start:start+4]
            
            # Save all 6hr time steps
            freeze_fct_temp[i, count:count+4] = stat_temp/3.28084 # Convert from ft to m
            
            # Mark file as present
            cnrfc_fct_missing[i, count:count+4] = 1
            
        else:
            print(str(qpf_file))
        
    # Increment count
    count += qpe_tsteps # 4

# Remove first and last 2 tsteps on each end
freeze_fct_6hr = freeze_fct_temp[:, 2:-2]


# ---------------------------------------------------------------------------
# Calculate radar obs
if radar_type == "FMCW":

    print("Calculating observed FMCW values") ###### Changed to go from 12Z to 12Z
    
    # Set window size (2, 4, 6, 8, 10, 20, 30, 36) in number of 10min time steps
    # windows = [2, 4, 8, 12, 16, 20, 24, 30, 36]
    # windows = range(2,38,2)
    windows = [2, 6, 18, 36]
    windows_minutes = [10 + (w-2)*10 for w in windows]
    num_windows = len(windows)
    
    # Missing data array
    fmcw_missing = np.zeros([num_windows, qpe_steps])
    
    # Get station elevation
    fmcw_elev = elevations.get(station)
    
    # # Number of time steps per day for FMCW data
    # fmcw_tsteps = 24*6
    
    # Create array to store output data
    # freeze_fmcw = np.full([fmcw_steps], np.nan)
    freeze_fmcw = np.full([num_windows, qpe_steps], np.nan)
    
    # Create array of times
    # minutes = (5, 15, 25, 35, 45, 55) # these are the 6 time steps within each hour
    
    # Loop through each hour
    count = 0
    
    # Loop through each valid date in newly created list
    # for v, valid_date in enumerate(valid_dates):
    for v in range(num_days):
        
        # Get valid date
        valid_date = valid_dates[v]
        
        # Convert date to string
        valid_str = datetime.strftime(valid_date, "%Y-%m-%d")
        
        # Import data
        print("FMCW Obs: Valid " + valid_str)
        
        # New method: windows
        # Need to loop through each 6 hr time step separately
        # Create list of the 4 6-hr time steps for the current day
        fmcw_time_steps = [valid_date + timedelta(hours=6*x) for x in range(4)]
        
        # Then, loop through each time step
        for s, step in enumerate(fmcw_time_steps):
            
            # print(s)
            
            # Temp counter within each time step
            temp_count = 0
            
            # Create temp array to store the full window of data
            window_temp = np.full(36, np.nan) # 4 six-hr time steps per day x 36 ten-min radar time steps per 6 hrs
            
            # Create list of the 3 hourly files on either side of the current time step
            fmcw_temp = [step + timedelta(hours=x) for x in range(-3,3)]
            # fmcw_files = [station.lower() + temp.strftime("%y%j%H") + ".snw" for temp in fmcw_temp]
            # fmcw_paths = [fmcw_dir / station / temp.strftime("%Y") / temp.strftime("%j") / temp for temp in fmcw_files]
            # fmcw_paths = [fmcw_dir / station / temp.strftime("%Y") / temp.strftime("%j") / (station.lower() + temp.strftime("%y%j%H") + ".snw") for temp in fmcw_temp]
            
            fmcw_paths = [fmcw_dir / station.lower() / temp.strftime("%Y") / temp.strftime("%j") / (station.lower() + temp.strftime("%y%j%H") + ".snw") for temp in fmcw_temp]
            
            # Missing flag
            # Check all 6 files
            fmcw_missing_flag = [0,0,0,0,0,0]
            
            # Loop through each of the 6 hourly files
            for f, fpath in enumerate(fmcw_paths):
                
                # Check if file exists
                if fpath.is_file():
                
                    # Read in current file
                    data = pd.read_fwf(fpath, header=None)
                    
                    # Extract mlh data
                    data2 = data.loc[:,1]
                    
                    # Find index range
                    idx = data2[data2 == "mlh"].index[0]
                    mlh_data = data.loc[idx+1:idx+6,1].values
                    
                    # Parse and append data to output array
                    for v_idx, value in enumerate(mlh_data):
                        # freeze_fmcw[count+v_idx] = float(value)
                        
                        # New line here
                        window_temp[temp_count+v_idx] = float(value)
                        
                    # # Mark file as present - come back to this later
                    # fmcw_missing[count:count+6] = 1
                    fmcw_missing_flag[f] = 1
                    
                else:
                    print(str(fpath))
                    
                # Increment index
                temp_count += 6
                
            
            # Replace -9.999 values with NaN
            window_temp[window_temp <= -9.9] = np.nan
            
            # Window size calculation is here!!!!!!!!
            #########################################
            # Loop through each window size
            for w, window in enumerate(windows):
                
                # Half width is half of window size
                half_width = int(window/2)
                
                # Starting index (halfway point of 36 elements)
                start_idx = 18 - half_width
                
                # Ending index
                end_idx = 18 + half_width
                
                # Slice array based on specified window size
                window_slice = window_temp[start_idx:end_idx]
                
                # Calculate median of window
                median_temp = np.nanmedian(window_slice)
                
                # Assign value to output array
                freeze_fmcw[w, count] = median_temp
                
                # Check missing values
                if w == 0 or w == 1:
                    if any(check == 1 for check in fmcw_missing_flag[2:4]):
                        fmcw_missing[w, count] = 1
                elif w == 2:
                    if any(check == 1 for check in fmcw_missing_flag[1:5]):
                        fmcw_missing[w, count] = 1
                elif w == 3:
                    if any(check == 1 for check in fmcw_missing_flag):
                        fmcw_missing[w, count] = 1
    
            # # Missing array
            # if fmcw_missing_flag > 0:
            #     fmcw_missing[count] = 1
    
            # Increment count (1 for each 6hr time step)
            count += 1
    
    # # Replace -9.999 values with NaN
    # freeze_fmcw[freeze_fmcw <= -9.9] = np.nan
    
    # Convert from km to m
    freeze_fmcw *= 1000
    
    # Convert from AGL to MSL
    freeze_fmcw += fmcw_elev


# Else: MRR
elif radar_type == "MRR":
    
    print("Calculating observed MRR values")
    
    # Windows
    windows = [2, 4, 6] # number of steps included
    windows_minutes = [w*30 for w in windows]
    num_windows = len(windows)
    
    # Get station elevation
    mrr_elev = elevations.get(station)
    
    # Raintype CSV header
    header = ["Date","BB Height","Raintype"]
    
    # Create array to store output data
    # freeze_mrr = np.full([mrr_steps], np.nan)
    freeze_mrr = np.full([num_windows, qpe_steps], np.nan)
    
    # Loop through each day
    count = 0 # count each 6hr time step
    
    # MRR indexs
    mrr_indexs = [24, 30, 36, 42]
    
    # Loop through each valid date in newly created list
    for v in range(num_days):
        
        # Get valid date
        valid_date = valid_dates[v]
        
        # Get preceding date as well - because for the time step at 0Z we want
        # the median of the 2 time steps at 23:30 and 00:30
        mrr_days = [valid_date - timedelta(days=1), valid_date]
        
        # Convert dates to string
        valid_str = [datetime.strftime(d, "%Y%m%d") for d in mrr_days]
        mrr_folder = [datetime.strftime(d, "%Y%m") for d in mrr_days]
        
        # Get 6 hr time steps
        # mrr_time_steps = [valid_date + timedelta(hours=6*x) for x in range(4)]
        # mrr_indexs = [23, 29, 35, 41]
        
        # Import data
        print("MRR Obs: Valid " + valid_str[1])
        
        # Get csv file containing brightband values
        # Need to get the file before this day as well
        mrr_files = [mrr_dir / station / "raintype_bbheight" / mrr_folder[m] / (station + "_" + valid_str[m] + "_raintype_hourly.csv") for m in range(2)]
        
        # Read in files
        if mrr_files[0].is_file():
            data_in = pd.read_csv(mrr_files[0], names=header)
            data1 = data_in.loc[:, "BB Height"]
        else:
            data1 = np.full([24], np.nan)
            print(str(mrr_files[0]))
        
        if mrr_files[1].is_file():
            data_in = pd.read_csv(mrr_files[1], names=header)
            data2 = data_in.loc[:, "BB Height"]
        else:
            data2 = np.full([24], np.nan)
            print(str(mrr_files[1]))
        
        # Concatenate
        data = np.concatenate((data1, data2))
        
        # Loop through each time step
        for s, step in enumerate(mrr_indexs):
            exist_flag = 0
            
            # Check whether file exists
            if s == 0: # check both files for first time step at 0Z
                if mrr_files[0].is_file() or mrr_files[1].is_file():
                    exist_flag = 1
            else: # just check second file
                if mrr_files[1].is_file():
                    exist_flag = 1
            
            # If exists, calculate median of 2 time steps
            if exist_flag == 1:
                
                # Loop through each window size
                for w, window in enumerate(windows):
                    
                    # Half width is half of window size
                    half_width = int(window/2)
                    
                    # Starting index
                    start_idx = step - half_width
                    
                    # Ending index
                    end_idx = step + half_width
                    
                    # Calculate median and assign to array
                    freeze_mrr[w, count] = np.nanmedian(data[start_idx:end_idx])
                
                # Mark file as present
                mrr_missing[count] = 1
            
            # Increment count
            count += 1
            
        #     # Read in csv
        #     # Check if file exists
        #     if mrr_file.is_file():
            
        #         # Read in current file
        #         data = pd.read_csv(mrr_file, names=header)
                
        #         # Extract brightband column
        #         bb_out = data.loc[:, "BB Height"]
                
        #         # Save to output array
        #         freeze_mrr[count:count+24] = bb_out
                
        #         # Mark file as present
        #         mrr_missing[count:count+24] = 1
                
        #     else:
        #         print(str(mrr_file))
            
        # # Increment count
        # count += 24
            
    # Convert from AGL to MSL
    freeze_mrr += mrr_elev

    
# ---------------------------------------------------------------------------
# Calculate QPE
print("Calculating QPE")

# Array for all 6 hr values (4 values per day)
qpe_obs_6hr = np.full([qpe_steps], np.nan)

# Counter
count = 0

# Loop through each valid date in newly created list
# for v, valid_date in enumerate(valid_dates):
for v in range(num_days):
    
    # Get valid date
    valid_date = valid_dates[v]
    
    # Convert date to string
    valid_str = datetime.strftime(valid_date, "%Y-%m-%d")
    
    # Import data
    print("QPE Obs: Valid " + valid_str)
    
    # Read in QPE watershed file
    # qrow1 = qpe_set[0]
    # qrow2 = qpe_set[1]
    # qcol1 = qpe_set[2]
    # qcol2 = qpe_set[3]

    # Create list of 6 hr time steps
    time_steps_qpe = [datetime.strftime(valid_date + timedelta(hours=6*x), "%Y%m%d_%H") for x in range(qpe_tsteps)]
    hour_files_qpe = ["qpe." + t + "00.nc" for t in time_steps_qpe]
    
    # Read in QPE
    for h, hour_file in enumerate(hour_files_qpe):
        
        # Full path
        qpe_file = qpedir / hour_file
        
        # Check if file exists
        if qpe_file.is_file():
        
            # Read in data
            with netCDF4.Dataset(qpe_file, mode='r') as infile:
                
                # Total precipitation
                # try:
                # obs_temp = infile.variables["qpe_grid"][0, qpe_idxs[0], qpe_idxs[1]]
                # except:
                #     obs_temp = infile.variables["tp"][0, 0, qrow1:qrow2, qcol1:qcol2]
                
                # Pre 2020 WY
                if start_year < 2020:
                    obs_temp = infile.variables["tp"][0, 0, qpe_idxs_old[0], qpe_idxs_old[1]]
                else:
                    obs_temp = infile.variables["qpe_grid"][0, qpe_idxs[0], qpe_idxs[1]]
                
            # Remove any fill values (-99999)
            # obs_temp[obs_temp < 0] = np.nan
            if obs_temp < 0:
                obs_temp = np.nan
            
            # # Get station idx
            # idx = qpe_idxs
            
            # # Get array subset
            # stat_temp = obs_temp[idx[0], idx[1]]
        
            # Store each time step to temp array
            qpe_obs_6hr[count+h] = obs_temp*25.4 # convert in to mm
            
            # Mark file as present
            cnrfc_qpe_missing[count+h] = 1
            
        else:
            print(str(qpe_file))
        
    # Increment count
    count += qpe_tsteps # 4
    
# Calculate cumulative sum
# qpe_obs_sum = np.cumsum(qpe_obs_6hr)


# ---------------------------------------------------------------------------
# Concatenate into all members
# freeze_all = np.concatenate((freeze_ecm, freeze_gefs), axis=0)

# # Calculate error (forecast - obs)
# freeze_nrt_error1 = freeze_nrt1 - freeze_obs_fct #9km
# freeze_nrt_error2 = freeze_nrt2 - freeze_obs_fct #3km
# cnrfc_error = freeze_fct_6hr - freeze_obs_fct


# ---------------------------------------------------------------------------
# Write to netCDF
print("Writing to netCDF")

# Get current time
today = datetime.utcnow()
date_str = datetime.strftime(today, "%Y-%m-%d %H:%M:%S")

# Variables to save:
    # freeze_thompson = np.full([num_stations, thompson_size, lead_days, num_days], np.nan)
    # freeze_morrison = np.full([num_stations, morrison_size, lead_days, num_days], np.nan)
    # freeze_ecm = np.full([num_stations, ecm_size, lead_days, num_days], np.nan)
    # freeze_gefs = np.full([num_stations, gefs_size, lead_days, num_days], np.nan)
    
    # freeze_obs_all = np.full([num_stations, ecm_size+gefs_size, lead_days, num_days], np.nan)
    # freeze_obs_thompson = np.full([num_stations, thompson_size, lead_days, num_days], np.nan)
    # freeze_obs_morrison = np.full([num_stations, morrison_size, lead_days, num_days], np.nan)
    # freeze_obs_ecm = np.full([num_stations, ecm_size, lead_days, num_days], np.nan)
    # freeze_obs_gefs = np.full([num_stations, gefs_size, lead_days, num_days], np.nan)
    # freeze_obs_fct = np.full([num_stations, lead_days, num_days], np.nan)
    # freeze_obs_6hr = np.full([num_stations, qpe_tsteps, num_days], np.nan)
    
    # freeze_fct = np.full([num_stations, lead_days, num_days], np.nan)
    # freeze_fmcw = np.full([num_stations, fmcw_tsteps, num_days], np.nan)

# Output file
outfile = ncdir / out_name

with netCDF4.Dataset(outfile, "w", format = "NETCDF4") as nc_out:
    
    # Define dimensions
    nc_out.createDimension("lead_days", lead_days)
    nc_out.createDimension("time", qpe_steps)
    nc_out.createDimension("num_windows", num_windows)
    
    # Add global attributes
    if wrf_flag == "y":
        nc_out.title = "Freezing level dataset: includes WWRF-NRT (both GFS and ECMWF forcings, for 2019-2021 only) vs. CNRFC freezing level forecasts, CNRFC freezing level obs, radar freezing level obs (FMCW or MRR), and CNRFC QPE precipitation"
    else:
        nc_out.title = "Freezing level dataset: includes CNRFC freezing level obs/forecast, radar freezing level obs (FMCW or MRR), and CNRFC QPE precipitation"
    nc_out.station = station + ": " + station_names.get(station)
    nc_out.valid_times = datetime.strftime(start_date, "%Y-%m-%d") + " to " + datetime.strftime(end_date, "%Y-%m-%d")
    nc_out.author = "Peter Yao"
    nc_out.history = "Created: " + date_str + " UTC"
    
    # Create variables and attributes
    # FMCW windows
    if radar_type == "FMCW":
        nc_window = nc_out.createVariable("fmcw_windows", "i4", ("num_windows"))
        nc_window.description = "Window size of FMCW observations for each 6 hr time step"
        nc_window.units = "minutes"
    elif radar_type == "MRR":
        nc_window = nc_out.createVariable("mrr_windows", "i4", ("num_windows"))
        nc_window.description = "Window size of MRR observations for each 6 hr time step"
        nc_window.units = "minutes"
    
    # Forecast lead time
    nc_time0 = nc_out.createVariable("lead_days", "i4", ("lead_days"))
    nc_time0.wrf_description = "For WWRF-NRT forecasts, ex. valid time June 4 0Z-24Z, lead day values of 4=init May 31 0Z, 3=init June 1 0Z, 2=init June 2 0Z, 1=init June 3 0Z"
    nc_time0.cnrfc_description = "For CNRFC forecasts, ex. valid time June 4 12Z-36Z, lead day values of 4=init May 31 12Z, 3=init June 1 12Z, 2=init June 2 12Z, 1=init June 3 12Z"
    nc_time0.units = "days"
    
    # Time arrays
    nc_time = nc_out.createVariable("time", "i4", ("time"))
    nc_time.description = "Valid time: time steps created to match CNRFC data (6 hour time steps)"
    nc_time.timezone = "UTC"
    nc_time.units = "seconds since 1970-01-01 00:00:00"
    
    # WRF outputs
    if wrf_flag == "y":
        nc1 = nc_out.createVariable("freeze_nrt_gfs_9km", "single", ("lead_days", "time"), fill_value = -9999)
        nc1.description = "WWRF-NRT GFS calculated freezing level: 9km domain"
        nc1.units = "m MSL"
        
        nc2 = nc_out.createVariable("freeze_nrt_gfs_3km", "single", ("lead_days", "time"), fill_value = -9999)
        nc2.description = "WWRF-NRT GFS calculated freezing level: 3km domain"
        nc2.units = "m MSL"
        
        nc1_2 = nc_out.createVariable("freeze_nrt_ecmwf_9km", "single", ("lead_days", "time"), fill_value = -9999)
        nc1_2.description = "WWRF-NRT ECMWF calculated freezing level: 9km domain"
        nc1_2.units = "m MSL"
        
        nc2_2 = nc_out.createVariable("freeze_nrt_ecmwf_3km", "single", ("lead_days", "time"), fill_value = -9999)
        nc2_2.description = "WWRF-NRT ECMWF calculated freezing level: 3km domain"
        nc2_2.units = "m MSL"
    
    # Obs outputs
    nc3 = nc_out.createVariable("freeze_obs", "single", ("time"), fill_value = -9999)
    nc3.description = "CNRFC observed freezing level"
    nc3.units = "m MSL"
    
    # CNRFC forecast output
    nc4 = nc_out.createVariable("freeze_fct", "single", ("lead_days", "time"), fill_value = -9999)
    nc4.description = "CNRFC forecast freezing level"
    nc4.units = "m MSL"
    
    # FMCW output
    if radar_type == "FMCW":
        nc13 = nc_out.createVariable("freeze_fmcw", "single", ("num_windows", "time"), fill_value = -9999)
        nc13.description = "FMCW freezing level"
        nc13.units = "m MSL"
    else:
        nc13 = nc_out.createVariable("freeze_mrr", "single", ("num_windows", "time"), fill_value = -9999)
        nc13.description = "MRR freezing level"
        nc13.units = "m MSL"
    
    # QPE
    nc14 = nc_out.createVariable("qpe_precip", "single", ("time"), fill_value = -9999)
    nc14.description = "CNRFC observed QPE precipitation"
    nc14.units = "mm"

    
    # Missing data outputs
    if wrf_flag == "y":
        nc_miss1 = nc_out.createVariable("freeze_nrt_gfs_9km_missing", "i4", ("lead_days", "time"))
        nc_miss1.description = "WWRF-NRT GFS calculated freezing level: 9km domain"
        nc_miss1.units = "1=present, 0=missing"
        
        nc_miss2 = nc_out.createVariable("freeze_nrt_gfs_3km_missing", "i4", ("lead_days", "time"))
        nc_miss2.description = "WWRF-NRT GFS calculated freezing level: 3km domain"
        nc_miss2.units = "1=present, 0=missing"
        
        nc_miss1_2 = nc_out.createVariable("freeze_nrt_ecmwf_9km_missing", "i4", ("lead_days", "time"))
        nc_miss1_2.description = "WWRF-NRT ECMWF calculated freezing level: 9km domain"
        nc_miss1_2.units = "1=present, 0=missing"
        
        nc_miss2_2 = nc_out.createVariable("freeze_nrt_ecmwf_3km_missing", "i4", ("lead_days", "time"))
        nc_miss2_2.description = "WWRF-NRT ECMWF calculated freezing level: 3km domain"
        nc_miss2_2.units = "1=present, 0=missing"
    
    # Obs outputs
    nc_miss3 = nc_out.createVariable("freeze_obs_missing", "i4", ("time"))
    nc_miss3.description = "CNRFC observed freezing level"
    nc_miss3.units = "1=present, 0=missing"
    
    # CNRFC forecast output
    nc_miss4 = nc_out.createVariable("freeze_fct_missing", "i4", ("lead_days", "time"))
    nc_miss4.description = "CNRFC forecast freezing level"
    nc_miss4.units = "1=present, 0=missing"
    
    # FMCW output
    if radar_type == "FMCW":
        nc_miss13 = nc_out.createVariable("freeze_fmcw_missing", "i4", ("num_windows", "time"))
        nc_miss13.description = "FMCW freezing level"
        nc_miss13.units = "1=present, 0=missing"
    else:
        nc_miss13 = nc_out.createVariable("freeze_mrr_missing", "i4", ("time"))
        nc_miss13.description = "MRR freezing level"
        nc_miss13.units = "1=present, 0=missing"
    
    # QPE
    nc_miss14 = nc_out.createVariable("qpe_precip_missing", "i4", ("time"))
    nc_miss14.description = "CNRFC observed QPE precipitation"
    nc_miss14.units = "1=present, 0=missing"
    
    
    # Latitude
    latitude = nc_out.createVariable('lat','f4')
    latitude.description = "Station latitude"
    latitude.units = "degrees_north"
    
    # Longitude
    longitude = nc_out.createVariable('lon','f4')
    longitude.description = "Station longitude"
    longitude.units = "degrees_east"
    
    # Station elevation
    elevation = nc_out.createVariable('site_elevation','i4')
    elevation.description = "Station elevation"
    elevation.units = "m MSL"
    
    
    # Assign data to variables
    nc_time0[:] = [4,3,2,1]
    nc_time[:] = qpe_time
    nc_window[:] = windows_minutes
    
    if wrf_flag == "y":
        nc1[:] = np.ma.masked_invalid(freeze_gfs1)
        nc2[:] = np.ma.masked_invalid(freeze_gfs2)
        nc1_2[:] = np.ma.masked_invalid(freeze_ecm1)
        nc2_2[:] = np.ma.masked_invalid(freeze_ecm2)
    nc3[:] = np.ma.masked_invalid(freeze_obs_6hr)
    nc4[:] = np.ma.masked_invalid(freeze_fct_6hr)
    
    if radar_type == "FMCW":
        nc13[:] = np.ma.masked_invalid(freeze_fmcw)
    else:
        nc13[:] = np.ma.masked_invalid(freeze_mrr)
    
    nc14[:] = np.ma.masked_invalid(qpe_obs_6hr)
    
    if wrf_flag == "y":
        nc_miss1[:] = gfs_missing1
        nc_miss2[:] = gfs_missing2
        nc_miss1_2[:] = ecm_missing1
        nc_miss2_2[:] = ecm_missing2
    nc_miss3[:] = cnrfc_obs_missing
    nc_miss4[:] = cnrfc_fct_missing
    if radar_type == "FMCW":
        nc_miss13[:] = fmcw_missing
    else:
        nc_miss13[:] = mrr_missing
    nc_miss14[:] = cnrfc_qpe_missing
    
    latitude[:] = lat
    longitude[:] = lon
    elevation[:] = elevations.get(station)


print("Processing finished")
    