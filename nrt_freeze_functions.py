# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:39:37 2021

Functions for calculating WWRF-NRT freezing level

@author: Peter Yao
"""

from datetime import datetime
# from scipy.interpolate import interp1d
from scipy.ndimage import zoom
import numpy as np
import netCDF4
import wrf

# ---------------------------------------------------------------------------
# Calculate freezing level at specified time step
def calc_freezing_level(wrf_file, calc_flag):
    
    # Calc_flag (0 = manual, 1 = wrf-python)
    
    # Read in data
    with netCDF4.Dataset(wrf_file, mode='r') as infile:
        
        # Check calculation type
        if calc_flag == 0: # Manual
            # Potential temperature
            th_temp1 = infile.variables["T"][0]
            
            # Pertubation pressure
            p_temp1 = infile.variables["P"][0]
            
            # Base state pressure
            pb_temp1 = infile.variables["PB"][0]
        
        else: # Auto
            # Temperature (degC)
            tc = wrf.getvar(infile, "tc").values
            
        # Variables used to calculate height
        # Pertubation geopotential
        ph_temp1 = infile.variables["PH"][0]
        
        # Base state geopotential
        phb_temp1 = infile.variables["PHB"][0]
        
        
    # Get dimensions
    bottom_top = np.size(ph_temp1, 0) - 1
    south_north = np.size(ph_temp1, 1)
    west_east = np.size(ph_temp1, 2)
    
    # # Interpolate to 1000 points
    # len1 = np.linspace(1,1000,bottom_top)
    # len2 = np.linspace(1,1000,1000)
    
    # Output array
    var_out = np.full([south_north, west_east], np.nan)
    
    # Convert geopotential height to height in meters
    # Calculate height (these are boundary edges, so n+1 from dbz)
    elev = (ph_temp1 + phb_temp1)/9.81
    
    # Interpolate to get elevation of the cell theta-points (midplane)
    # Height is from low to high
    elev_theta = (elev[1:] + elev[:-1])/2
            
    # Interpolate to 1000 points (cubic spline)
    tc_interp = zoom(tc, (1000/bottom_top,1,1))
    elev_interp = zoom(elev_theta, (1000/bottom_top,1,1))
    
    # Loop through each lat/lon coordinate
    for idx1 in range(south_north):
        for idx2 in range(west_east):
            
            # Get values at current coordinates
            tc_temp = tc_interp[:, idx1, idx2]
            elev_temp = elev_interp[:, idx1, idx2]
            
            # # Get values at current coordinates
            # ph_temp = ph_temp1[:, idx1, idx2]
            # phb_temp = phb_temp1[:, idx1, idx2]
    
            # # Manually calculate temperature in degC
            # if calc_flag == 0:
            #     th_temp = th_temp1[:, idx1, idx2]
            #     p_temp = p_temp1[:, idx1, idx2]
            #     pb_temp = pb_temp1[:, idx1, idx2]
                
            #     # Manually calculate temperature in degC
            #     tc_temp = (th_temp+300)/((1000/((p_temp+pb_temp)/100))**(2/7)) - 273.15
            
            # # Or retrieve temperature at current coordinates
            # else:
            #     tc_temp = tc[:, idx1, idx2]
            
            # # Convert geopotential height to height in meters
            # # Calculate height (these are boundary edges, so n+1 from dbz)
            # elev = (ph_temp + phb_temp)/9.81
            
            # # Interpolate to get elevation of the cell theta-points (midplane)
            # # Height is from low to high
            # elev_theta = (elev[1:] + elev[:-1])/2

            # # Interpolate to 1000 points
            # # len1 = np.linspace(1,1000,bottom_top)
            # # len2 = np.linspace(1,1000,1000)
            # tc_intpfunc = interp1d(len1, tc_temp, kind="cubic")
            # tc_interp = tc_intpfunc(len2)
            # elev_intpfunc = interp1d(len1, elev_theta, kind="cubic")
            # elev_interp = elev_intpfunc(len2)
    
            # Calculate freezing level
            # Find lowest height when temperature = 0 or crosses zero
            # First, check if any values = 0
            zero_idx = np.where(tc_temp == 0)[0]
            if len(zero_idx) > 0:
                zero_idx = zero_idx[0]
                
                var_out[idx1, idx2] = elev_temp[zero_idx]
            
            # More likely case that we need to find the sign change at the lowest height
            else:
                tc_cross = np.where(np.sign(tc_temp[:-1]) != np.sign(tc_temp[1:]))[0]
                if len(tc_cross) > 0:
                    zero_idx = tc_cross[0]
                    
                    # Assign 3-hr value to temp array
                    # freeze_temp[h, s_idx] = elev_interp[zero_idx]
                    var_out[idx1, idx2] = elev_temp[zero_idx]
                
                else:
                    # # Get file name
                    # temp_file = str(wrf_file)
                    
                    # # Find min temp/elevation
                    # temp_min = np.round(np.nanmin(tc_interp), 2)
                    # temp_idx = np.argmin(tc_interp)
                    # temp_elev = np.round(elev_interp[temp_idx], 2)
                    
                    # print("Min: " + str(temp_min) + "C at " + str(temp_elev) + " m: " + temp_file)
                    var_out[idx1, idx2] = np.nan

    return var_out


# ---------------------------------------------------------------------------
# Write freezing level output to netCDF
def write_netcdf(freeze_out, xlat, xlon, outfile, init_date, valid_time, domain, fhours, forcing):
    
    # Get current time
    time = datetime.utcnow()
    time_str = datetime.strftime(time, "%Y-%m-%d %H:%M:%S UTC")
    
    # Get init time
    init_str = datetime.strftime(init_date, "%Y%m%d_%H0000")
    
    # Get valid time
    valid_str = datetime.strftime(valid_time, "%Y%m%d_%H0000")
    
    # Define fill value
    fillV = 9.999000260554009e20
    
    with netCDF4.Dataset(outfile, "w", format = "NETCDF4") as nc_out:
    
        # Define dimensions
        nc_out.createDimension("south_north", np.size(xlat,0))
        nc_out.createDimension("west_east", np.size(xlat,1))
        nc_out.createDimension("time", 1)
        
        # Add global attributes
        nc_out.title = "Calculated WWRF-NRT freezing level"
        nc_out.institution = "Center for Western Weather and Water Extremes (CW3E), Scripps Institution of Oceanography"
        nc_out.history = "Created on " + time_str
        nc_out.author = "Peter Yao"
        nc_out.domain = domain
        nc_out.init_time = init_str
        nc_out.valid_time = valid_str
        nc_out.forcing = forcing
        
        # Create variables and attributes
        nc_freeze = nc_out.createVariable("Z0C", "single", ("south_north", "west_east"), fill_value = fillV)
        nc_freeze.standard_name = "freezing_level_altitude"
        nc_freeze.description = "Calculated freezing level, or zero-degree isotherm (in m MSL)"
        nc_freeze.units = "m"
        nc_freeze.coordinates = "XLONG XLAT"
        
        nc_lat = nc_out.createVariable("XLAT", "single", ("south_north", "west_east"))
        nc_lat.standard_name = "latitude"
        nc_lat.description = "Latitude, south is negative"
        nc_lat.units = "degrees_north"
        nc_lat.coordinates = "XLONG XLAT"
        
        nc_lon = nc_out.createVariable("XLONG", "single", ("south_north", "west_east"))
        nc_lon.standard_name = "longitude"
        nc_lon.description = "Longitude, west is negative"
        nc_lon.units = "degrees_east"
        nc_lon.coordinates = "XLONG XLAT"

        nc_time = nc_out.createVariable("time", "i4", ("time"))
        nc_time.standard_name = "time"
        nc_time.long_name = "Number of hours since beginning of valid date"
        nc_time.description = "hours since " + datetime.strftime(valid_time,'%Y-%m-%d') + " 00:00:00"
        nc_time.units = "hours since " + datetime.strftime(valid_time,'%Y-%m-%d') + " 00:00:00"
        
        nc_ftime = nc_out.createVariable("forecast_lead_time", "i4", ("time"))
        nc_ftime.long_name = "Number of hours since beginning of initialization date"
        nc_ftime.description = "hours since " + datetime.strftime(init_date,'%Y-%m-%d') + " 00:00:00"
        nc_ftime.units = "hours since " + datetime.strftime(init_date,'%Y-%m-%d') + " 00:00:00"
        
        
        # Assign data to variables
        nc_freeze[:] = np.ma.masked_invalid(freeze_out)
        nc_lat[:] = xlat
        nc_lon[:] = xlon
        nc_time[:] = valid_time.hour
        nc_ftime[:] = fhours*3
        
        
        