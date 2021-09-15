# -*- coding: utf-8 -*-
"""
Created on Wed May 19 09:24:36 2021

Script to calculate the freezing level for the WWRF-NRT for a given water year
Saves outputs to netCDF files, matching the format from the GFS files created by Brian

V2: fixes error

@author: Peter Yao
"""

from pathlib import Path
from datetime import datetime, timedelta
import os
import argparse
import numpy as np
import netCDF4
import nrt_freeze_functions as ffun

# Time it
time1 = datetime.now()

#--------------------------------------------------------------------------
# User flags

# Water year (format: "YYYY-YYYY")
year = "2019-2020"

# Domain (1=9km, 2=3km)
# domain = 2

# Forcing (gfs or ecmwf)
# forcing = "gfs"

# Calculation flag (0=calculate temperature manually, 1=use wrf-python)
calc_flag = 1

# Input directory (containing WWRF-NRT files)
wrfdir = Path("/data/downloaded/WWRF-NRT") / year / "WRF_Output"

# Output directory (where to store output netCDF files)
outdir = Path("/data/downloaded/Forecasts/Model_Verification/WWRF-NRT/Z0C")


#--------------------------------------------------------------------------
# Initial setup

# Path to script name
script = "/home/peyao/Python_Scripts/wrf_nrt_calc_freeze.py"

# Parser controls in the input options from command line.
parser = argparse.ArgumentParser(description="Enter run number")
parser.add_argument('-d', '--domain', required=True, help="Domain: format 1 or 2")
parser.add_argument('-n', '--total', required=True, help="Total number of runs: format Y")
parser.add_argument('-r', '--run', required=True, help="Run number: format X of Y, starting at 1")
parser.add_argument('-f', '--forcing', required=True, help="Forcing: gfs or ecmwf")
args = parser.parse_args()

domain = int(args.domain) # format: 1 or 2
total = int(args.total) # format: number
run = int(args.run) # format: number
forcing = args.forcing # format: string

print("Year: " + year)
print("Domain: " + str(domain))
print("Processing: run " + str(run) + " of " + str(total))

# Specify domain name, number of forecast days
# Also specify grid file - depends on domain
if domain == 1:
    dlen = 7
    dname = "9km"
    if forcing == "gfs":
        grid_file = wrfdir / "2019112400" / "gfs" / "wrfout_d01_2019-11-24_00_00_00.nc"
    else:
        grid_file = wrfdir / "2019112400" / "ecmwf" / "wrfout_d01_2019-11-24_00_00_00.nc"
else:
    dlen = 5
    dname = "3km"
    if forcing == "gfs":
        grid_file = wrfdir / "2019112400" / "gfs" / "wrfout_d02_2019-11-24_00_00_00.nc"
    else:
        grid_file = wrfdir / "2019112400" / "ecmwf" / "wrfout_d02_2019-11-24_00_00_00.nc"
    
# Number of 3 hour time steps
tsteps = dlen*8 + 1

# Specify output path
out_path = outdir / forcing / dname / year

# Create output directory if it does not exist
if not out_path.is_dir():
    out_path.mkdir(parents=True, exist_ok=True)
    
# Import lat/lon array
with netCDF4.Dataset(grid_file, mode='r') as infile:
    
    # Lat/lon
    xlat = infile.variables["XLAT"][0]
    xlon = infile.variables["XLONG"][0]
    
# lat_size = np.size(xlat,0)
# lon_size = np.size(xlat,1)

# Output dimensions
# south_north = np.size(xlat, 0)
# west_east = np.size(xlat, 1)
#bottom_top = np.size(elev_theta, 1) # need to update these to make robust


# --------------------------------------------------------------------------
# Main loop
print("Main loop: calculating WRF freezing level")

# Get subfolder list (list of init dates)
folder_list = sorted(os.listdir(wrfdir))
# folder_list = ["2020111900"]

# Divide list of dates (parallelize)
split_dates = np.array_split(folder_list, total)

# # Convert to list of lists
# split_dates = []
# for i in split:
#    split_dates.append(list(i))

# # Time it
# time1 = datetime.now()

# Loop through each init date
for date in split_dates[run-1]:
    
    # # Time it
    # time3 = datetime.now()

    # Check correct folder name format (e.g. "2020111500")
    if len(date) == 10 and date.isdigit():
        
        # Only proceed if specified forcing exists
        forcing_folder = wrfdir / date / forcing
        if forcing_folder.is_dir():
        
            # Get current init date
            year = int(date[0:4])
            month = int(date[4:6])
            day = int(date[6:8])
            init_date = datetime(year, month, day)
            date_str = datetime.strftime(init_date, "%Y-%m-%d")
            folder_str = datetime.strftime(init_date, "%Y%m%d00")
            
            # Print to console
            print("Processing init date: " + date_str)
            
            # Create list of file names for each valid date
            file_names = ["wrfout_d0" + str(domain) + "_" + datetime.strftime(init_date + timedelta(hours=3*x), "%Y-%m-%d_%H") + "_00_00.nc" for x in range(tsteps)]
            valid_dates = [init_date + timedelta(hours=3*x) for x in range(tsteps)]
            
            # Specify output folder
            out_folder = out_path / folder_str
            
            # Create output directory if it does not exist
            if not out_folder.is_dir():
                out_folder.mkdir(parents=True, exist_ok=True)
            
            # Loop through each valid date
            for f, fname in enumerate(file_names):
                
                # Get input file path
                wrf_file = forcing_folder / fname
                
                # Output file path
                outfile = out_folder / ("WWRF-NRT_Z0C_d0" + str(domain) + "_" + folder_str + "_F" + str(f*3).zfill(3) + ".nc")
                
                # Get valid time
                valid_time = valid_dates[f]
                
                # Check if file exists - and that output does not exist
                if wrf_file.is_file() and not outfile.is_file():
                    
                    # Time it
                    time3 = datetime.now()
                    
                    # Use try/except block to catch errors
                    try:
                        # Calculate freezing level
                        freeze_out = ffun.calc_freezing_level(wrf_file, calc_flag)
                        
                        # Write output to netCDF
                        # outfile = out_folder / ("WWRF-NRT_Z0C_d0" + str(domain) + "_" + folder_str + "_F" + str(f*3).zfill(3) + ".nc")
                        ffun.write_netcdf(freeze_out, xlat, xlon, outfile, init_date, valid_time, dname, f, forcing)
                    
                    except:
                        print("Error: " + str(wrf_file))
                    
                    # End time
                    time4 = datetime.now()
                    runtime = (time4-time3).total_seconds()
                    
                    print("Time: " + str(runtime) + " sec")    
                
                else:
                    print("File not found")

# End time
time2 = datetime.now()
runtime = (time2-time1).total_seconds()
    
print("Total time: " + str(runtime) + " sec")            
print("Processing finished")

