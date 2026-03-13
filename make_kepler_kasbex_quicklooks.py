#!/usr/bin/env python

import getopt, sys, os
import re
import copy

import datetime

import netCDF4 as nc4

import pyart
import numpy as np
import numpy.ma as ma
import shutil
import glob
import gzip

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib import colors  # Add this import
import cmocean
import getpass, socket

import pandas as pd

import cftime

# Import kepler utilities
from kepler_utils import get_valid_sweep_indices

version = 0.1

from pathlib import Path
homepath = Path.home()

matplotlib.use('Agg')

import cartopy

import geopandas as gpd
import fiona
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import contextily as ctx
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import pyproj

from pyproj import Geod

from contextily.tile import warp_img_transform, warp_tiles, _warper

os_key = "ngaMRdeOtnGTnXCb378JLN29j3H8AWAo"

from cartopy.io.img_tiles import GoogleTiles as GT

# Define KML file paths
kml_paths = {
    "verticals": "/home/users/cjwalden/WesConGrid_Verticals.kml",
    "horizontals": "/home/users/cjwalden/WesConGrid_Horizontals.kml",
    "outline": "/home/users/cjwalden/WesConGrid_Outline.kml",
    "surface_sites": "/home/users/cjwalden/surface_sites.kml"
}

# Define colormaps globally so all functions can access them
COLORMAPS = {
    'dbz': 'HomeyerRainbow',
    'vel': 'balance',
    'width': 'SpectralExtended',
    'ldr': 'SpectralExtended'
}

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:o:bkzm:a:t:", 
                              ["date=","inpath=","outpath=","dbz-only","basemap=","azimuth-offset=",
                               "min-lat=","max-lat=","min-lon=","max-lon=","max-age=","skip-all-transition"])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.datetime.now()
datestr = data_date.strftime('%Y%m%d')

tracking_tag = 'AMOF_20250508133639';

campaign = 'kasbex';

yaml_project_file = os.path.join(homepath,'amof_campaigns',f'{campaign}_project.yml')
yaml_instrument_file = os.path.join(homepath,'amof_instruments','amof_instruments.yml')

# Set default paths (will be updated after parsing date)
inpath = None  # Will be set after parsing datestr
figpath = f'/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/kasbex/L1c/quicklooks'

blflag = False;
darkmode = False;
dbz_only = False
basemap_type = "opentopo"  # Default
basemap_type_options = ["satellite", "opentopo", "cartodb"]
azimuth_offset = 0.0  # Default azimuth offset
max_file_age_hours = None  # Default: no age limit
skip_all_transition = False  # Default: don't skip all-transition sweeps

# Default lat/lon bounds for KASBEX at CAO
lat_min = 50.75
lat_max = 51.5
lon_min = -2.0
lon_max = -0.85

for o, a in opts:
    if o == "-d":
        datestr = a
    elif o == "-i":
        inpath = a;
    elif o == "-o":
        figpath = a;
    elif o == "-b":
        blflag = True;
    elif o == "-k":
        darkmode = True;
    elif o == "-z" or o == "--dbz-only":
        dbz_only = True
    elif o == "-m" or o == "--basemap":
        basemap_type = a  # Options: "satellite", "opentopo", "cartodb"
    elif o == "-a" or o == "--azimuth-offset":
        azimuth_offset = float(a)  # Convert to float
    elif o == "-t" or o == "--max-age":
        max_file_age_hours = float(a)  # Convert to float (hours)
    elif o == "--min-lat":
        lat_min = float(a)
    elif o == "--max-lat":
        lat_max = float(a)
    elif o == "--min-lon":
        lon_min = float(a)
    elif o == "--max-lon":
        lon_max = float(a)
    elif o == "--skip-all-transition":
        skip_all_transition = True
    else:
        assert False, "unhandled option"

# Set default inpath after parsing datestr (if not overridden by -i)
if inpath is None:
    inpath = f'/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/kasbex/L1c/{datestr}'

print(f"Processing date: {datestr}")
print(f"Input path: {inpath}")
print(f"Output path: {figpath}")
print(f"PPI map bounds: lat [{lat_min:.3f}, {lat_max:.3f}], lon [{lon_min:.3f}, {lon_max:.3f}]")
if max_file_age_hours is not None:
    print(f"Maximum file age: {max_file_age_hours} hours")
else:
    print("No maximum file age limit set")
print(f"Skip all-transition sweeps: {'enabled' if skip_all_transition else 'disabled'}")

def filter_files_by_age(files, max_age_hours=None):
    """
    Filter files by their modification time.
    
    Args:
        files (list): List of file paths
        max_age_hours (float): Maximum age in hours, None for no limit
        
    Returns:
        list: Filtered list of files within the age limit
    """
    if max_age_hours is None:
        return files
    
    current_time = datetime.datetime.now()
    max_age_seconds = max_age_hours * 3600
    filtered_files = []
    
    for file_path in files:
        try:
            # Get file modification time
            file_mtime = datetime.datetime.fromtimestamp(os.path.getmtime(file_path))
            file_age_seconds = (current_time - file_mtime).total_seconds()
            
            if file_age_seconds <= max_age_seconds:
                filtered_files.append(file_path)
            else:
                print(f"Skipping old file: {os.path.basename(file_path)} (age: {file_age_seconds/3600:.1f} hours)")
                
        except Exception as e:
            print(f"Error checking age of {file_path}: {e}")
            # Include file if we can't check its age
            filtered_files.append(file_path)
    
    if max_age_hours is not None:
        print(f"Age filter: {len(filtered_files)}/{len(files)} files within {max_age_hours} hours")
    
    return filtered_files

def filter_files_by_scan_age(files, max_age_hours=None):
    """
    Filter files by their radar scan time (extracted from filename or file data).
    
    Args:
        files (list): List of file paths
        max_age_hours (float): Maximum age in hours, None for no limit
        
    Returns:
        list: Filtered list of files within the age limit
    """
    if max_age_hours is None:
        return files
    
    current_time = datetime.datetime.utcnow()
    max_age_seconds = max_age_hours * 3600
    filtered_files = []
    
    for file_path in files:
        try:
            # First try to extract scan time from filename
            # Expected format: *YYYYMMDD-HHMMSS*
            filename = os.path.basename(file_path)
            
            # Look for YYYYMMDD-HHMMSS pattern in filename
            import re
            time_pattern = r'(\d{8})-(\d{6})'
            match = re.search(time_pattern, filename)
            
            scan_time = None
            
            if match:
                date_str = match.group(1)  # YYYYMMDD
                time_str = match.group(2)  # HHMMSS
                
                try:
                    scan_time = datetime.datetime.strptime(f"{date_str}-{time_str}", "%Y%m%d-%H%M%S")
                except ValueError:
                    pass
            
            # If filename parsing failed, try to read from netCDF file
            if scan_time is None:
                try:
                    with nc4.Dataset(file_path, 'r') as ds:
                        # Try to get time from netCDF attributes or variables
                        if hasattr(ds, 'time_coverage_start'):
                            # ISO format: YYYY-MM-DDTHH:MM:SSZ
                            time_str = ds.time_coverage_start
                            if time_str.endswith('Z'):
                                time_str = time_str[:-1]
                            scan_time = datetime.datetime.fromisoformat(time_str.replace('T', ' '))
                        elif 'time' in ds.variables:
                            # Use first time value
                            import cftime
                            time_data = ds.variables['time']
                            if len(time_data) > 0:
                                scan_time = cftime.num2pydate(time_data[0], time_data.units)
                                # Convert to naive datetime if it's timezone-aware
                                if hasattr(scan_time, 'tzinfo') and scan_time.tzinfo is not None:
                                    scan_time = scan_time.replace(tzinfo=None)
                except Exception as e:
                    print(f"Could not read scan time from {filename}: {e}")
            
            # If we still don't have scan time, fall back to file modification time
            if scan_time is None:
                print(f"Using file modification time for {filename} (could not extract scan time)")
                scan_time = datetime.datetime.fromtimestamp(os.path.getmtime(file_path))
            
            # Calculate age based on scan time
            scan_age_seconds = (current_time - scan_time).total_seconds()
            
            if scan_age_seconds <= max_age_seconds:
                filtered_files.append(file_path)
            else:
                print(f"Skipping old scan: {filename} (scan age: {scan_age_seconds/3600:.1f} hours, scan time: {scan_time.strftime('%Y-%m-%d %H:%M:%S')})")
                
        except Exception as e:
            print(f"Error checking scan age of {file_path}: {e}")
            # Include file if we can't check its age
            filtered_files.append(file_path)
    
    if max_age_hours is not None:
        print(f"Scan age filter: {len(filtered_files)}/{len(files)} files within {max_age_hours} hours of current time")
    
    return filtered_files

def ensure_increasing_azimuths(radar):
    """
    Checks if azimuths are increasing and reverses them if necessary.

    Parameters:
    radar (pyart.core.Radar): The original Py-ART radar object.

    Returns:
    pyart.core.Radar: A new radar object with azimuths in decreasing order.
    """
    # Deep copy to avoid modifying the original radar object
    radar_adjusted = copy.deepcopy(radar)

    # Check if azimuths are decreasing
    #if np.all(np.diff(radar.azimuth["data"]) > 0):  
    if radar.azimuth["data"][1] < radar.azimuth["data"][0]:
        print("Azimuths are decreasing, reversing them...")

        # Reverse azimuth array
        radar_adjusted.azimuth["data"] = radar.azimuth["data"][::-1]

        # Reverse all data fields to maintain alignment
        for field_name in radar.fields:
            radar_adjusted.fields[field_name]["data"] = radar.fields[field_name]["data"][::-1]

    else:
        print("Azimuths are already in increasing order. No changes made.")

    return radar_adjusted


def get_average_ray_duration(radar):
    """
    Calculate the average ray duration in milliseconds from radar time data.
    
    Args:
        radar: PyART radar object
        
    Returns:
        float: Average ray duration in milliseconds, or None if cannot calculate
    """
    try:
        # Get time data in seconds
        time_data = radar.time['data']
        
        if len(time_data) < 2:
            return None
            
        # Calculate time differences between consecutive rays
        time_diffs = np.diff(time_data)
        
        # Remove any negative differences (shouldn't happen but just in case)
        time_diffs = time_diffs[time_diffs > 0]
        
        if len(time_diffs) == 0:
            return None
            
        # Calculate average ray duration in milliseconds
        avg_ray_duration_ms = np.mean(time_diffs) * 1000
        
        return avg_ray_duration_ms
        
    except Exception as e:
        print(f"Could not calculate ray duration: {e}")
        return None

def make_dymecs_rhi_plot(ncfile, figpath, region, blflag=False, darkmode=False, azimuth_offset=0.0, skip_all_transition=False):
    if darkmode:
        plt.style.use('dark_background')

    if blflag:
        hmax = 4
        xmin = 0
        xmax = 25
    else:
        hmax = 12
        xmin = 0
        xmax = 40

    dbz_cmap = 'HomeyerRainbow'
    vel_cmap = 'balance'
    ldr_cmap = 'SpectralExtended'
    spw_cmap = 'SpectralExtended'
    snr_cmap = 'Carbone42'  # Good colormap for SNR

    # Validate the time dimension before reading the file
    try:
        with nc4.Dataset(ncfile, 'r') as ds:
            if 'time' not in ds.dimensions or ds.dimensions['time'].size == 0:
                print(f"Skipping file {ncfile}: 'time' dimension is missing or has length 0.")
                return
    except Exception as e:
        print(f"Error validating 'time' dimension in {ncfile}: {e}")
        return

    # Read the file with Py-ART
    try:
        RadarDS = pyart.io.read_cfradial(ncfile, delay_field_loading=True)
    except Exception as e:
        print(f"Error reading file {ncfile} with Py-ART: {e}")
        return

    # Validate elevation data
    if RadarDS.elevation['data'].size == 0:
        print(f"Skipping file {ncfile}: Elevation data is missing or has size 0.")
        return

    # Skip files where all values of antenna_transition are 1
    if RadarDS.antenna_transition is not None:
        if np.all(RadarDS.antenna_transition == 1):
            print(f"Skipping file {ncfile}: All values of 'antenna_transition' are 1.")
            return

    # Calculate average ray duration
    avg_ray_duration = get_average_ray_duration(RadarDS)
    if avg_ray_duration is not None:
        ray_duration_str = f"{avg_ray_duration:.1f}ms"
    else:
        ray_duration_str = "unknown"

    # Apply azimuth offset
    if azimuth_offset != 0.0:
        RadarDS.azimuth['data'] += azimuth_offset
        RadarDS.azimuth['data'] = RadarDS.azimuth['data'] % 360.0

    # Assign dtime_sweep before using it
    #dtime_sweep = cftime.num2pydate(RadarDS.time['data'][0], RadarDS.time['units'])
    #dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S")
    nsweeps = RadarDS.nsweeps
    
    # Get valid sweep indices (optionally filtering out all-transition sweeps)
    valid_sweep_indices = get_valid_sweep_indices(RadarDS, skip_all_transition=skip_all_transition)
    
    if skip_all_transition and len(valid_sweep_indices) < nsweeps:
        print(f"Skipping {nsweeps - len(valid_sweep_indices)} sweep(s) that are 100% antenna transitions")

    vel_field = RadarDS.fields['VEL']
    vel_limit_lower = vel_field['field_limit_lower']
    vel_limit_upper = vel_field['field_limit_upper']

    figpath = os.path.join(figpath, 'rhi', datestr)
    if not os.path.isdir(figpath):
        os.makedirs(figpath)

    if nsweeps > 0:
        for s in valid_sweep_indices:
            Radar = RadarDS.extract_sweeps([s])
            display = pyart.graph.RadarDisplay(Radar)
            gatefilter = pyart.correct.GateFilter(Radar)

            dtime_sweep = cftime.num2pydate(Radar.time['data'][0], Radar.time['units'])
            dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S")

            print(f"sweep {s+1}/{nsweeps}")
            rhi_az = Radar.get_azimuth(0)[0]

            fig, ax = plt.subplots(5, 1, figsize=(15, 22), constrained_layout=True)
            fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2, wspace=0.2)

            display.plot_rhi("DBZ", ax=ax[0], sweep=0, vmin=-60, vmax=40, gatefilter=gatefilter,
                             cmap=dbz_cmap, colorbar_orient='horizontal',
                             title=f"Radar equivalent reflectivity factor (dBZ)\nAzimuth: {rhi_az:.1f}° | Date: {dtime_sweep.strftime('%Y-%m-%d')} | Time: {dtime_sweep.strftime('%H:%M:%S')} | Ray duration: {ray_duration_str}")
            ax[0].set_xlabel("Distance from radar (km)")
            ax[0].set_ylabel("Distance above radar (km)")
            ax[0].set_ylim(0, hmax)
            ax[0].set_xlim(xmin, xmax)
            ax[0].grid(True)
            ax[0].set_aspect('equal', 'box')

            display.plot_rhi("VEL", ax=ax[1], sweep=0, vmin=vel_limit_lower, vmax=vel_limit_upper, gatefilter=gatefilter,
                             cmap=vel_cmap, colorbar_orient='horizontal',
                             title=f"Radial velocity of scatterers (m/s)\nAzimuth: {rhi_az:.1f}° | Date: {dtime_sweep.strftime('%Y-%m-%d')} | Time: {dtime_sweep.strftime('%H:%M:%S')} | Ray duration: {ray_duration_str}")
            ax[1].set_xlabel("Distance from radar (km)")
            ax[1].set_ylabel("Distance above radar (km)")
            ax[1].set_ylim(0, hmax)
            ax[1].set_xlim(xmin, xmax)
            ax[1].grid(True)
            ax[1].set_aspect('equal', 'box')

            display.plot_rhi("WIDTH", ax=ax[2], sweep=0,
                             norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1), vmax=np.sqrt(1e1)), gatefilter=gatefilter,
                             cmap=spw_cmap, colorbar_orient='horizontal',
                             title=f"Radar Doppler spectrum width (m/s)\nAzimuth: {rhi_az:.1f}° | Date: {dtime_sweep.strftime('%Y-%m-%d')} | Time: {dtime_sweep.strftime('%H:%M:%S')} | Ray duration: {ray_duration_str}")
            ax[2].set_ylim(0, hmax)
            ax[2].set_xlim(xmin, xmax)
            ax[2].grid(True)
            ax[2].set_aspect('equal', 'box')

            display.plot_rhi("LDR", ax=ax[3], sweep=0, vmin=-35, vmax=5, gatefilter=gatefilter,
                             cmap=ldr_cmap, colorbar_orient='horizontal',
                             title=f"Radar linear depolarization ratio (dB)\nAzimuth: {rhi_az:.1f}° | Date: {dtime_sweep.strftime('%Y-%m-%d')} | Time: {dtime_sweep.strftime('%H:%M:%S')} | Ray duration: {ray_duration_str}")
            ax[3].set_ylim(0, hmax)
            ax[3].set_xlim(xmin, xmax)
            ax[3].grid(True)
            ax[3].set_aspect('equal', 'box')

            display.plot_rhi("SNR", ax=ax[4], sweep=0, vmin=-20, vmax=60,
                             cmap=snr_cmap, colorbar_orient='horizontal',
                             title=f"Radar signal-to-noise ratio (dB)\nAzimuth: {rhi_az:.1f}° | Date: {dtime_sweep.strftime('%Y-%m-%d')} | Time: {dtime_sweep.strftime('%H:%M:%S')} | Ray duration: {ray_duration_str}")
            ax[4].set_ylim(0, hmax)
            ax[4].set_xlim(xmin, xmax)
            ax[4].grid(True)
            ax[4].set_aspect('equal', 'box')

            figname = f'ncas-mobile-ka-band-radar-1_cao_{dtime_sweep_str}_rhi_az{rhi_az:0.2f}_l1_v1.0.0.png'
            plt.savefig(os.path.join(figpath, figname), dpi=300)
            plt.close();



def make_dymecs_rhi_plots_day(datestr, inpath, figpath, blflag=False, darkmode=False, azimuth_offset=0.0, max_file_age_hours=None, skip_all_transition=False):
    inpath_date = inpath

    os.chdir(inpath_date)
    print(inpath_date)
    
    # Display max age setting
    if max_file_age_hours is not None:
        print(f"Maximum file age: {max_file_age_hours} hours")
    else:
        print("No maximum file age limit set")
    
    # Filter for L1 data version files only
    import re
    
    # Find all potential RHI files
    all_rhi_files = [os.path.join(inpath_date, f) for f in glob.glob('*{}*_rhi*.nc'.format(datestr))]
    
    # Filter for L1 data version pattern: _l1_v\d+\.\d+\.\d+\.nc
    l1_pattern = re.compile(r'_l1_v\d+\.\d+\.\d+\.nc$')
    woest_rhi_files = [f for f in all_rhi_files if l1_pattern.search(f)]
    
    # Filter files by filename timestamp if max_file_age_hours is specified
    if max_file_age_hours is not None:
        current_time_utc = datetime.datetime.utcnow()  # Explicitly use UTC
        max_age_seconds = max_file_age_hours * 3600
        filtered_files = []
        
        print(f"Current UTC time: {current_time_utc.strftime('%Y-%m-%d %H:%M:%S')} UTC")
        print(f"Filtering files older than {max_file_age_hours} hours...")
        
        # Pattern to match your filename format: YYYYMMDD-HHMMSS
        time_pattern = r'(\d{8})-(\d{6})'  # Match YYYYMMDD-HHMMSS in filenames
        
        for f in woest_rhi_files:
            filename = os.path.basename(f)
            match = re.search(time_pattern, filename)
            
            if match:
                date_str = match.group(1)  # YYYYMMDD
                time_str = match.group(2)  # HHMMSS
                try:
                    # Parse the timestamp as UTC (since filenames contain UTC timestamps)
                    file_time_utc = datetime.datetime.strptime(f"{date_str}-{time_str}", "%Y%m%d-%H%M%S")
                    file_age_seconds = (current_time_utc - file_time_utc).total_seconds()
                    file_age_hours = file_age_seconds / 3600
                    
                    if file_age_seconds <= max_age_seconds:
                        filtered_files.append(f)
                        print(f"✓ Including file {filename}: age {file_age_hours:.1f} hours")
                    else:
                        print(f"✗ Skipping file {filename}: age {file_age_hours:.1f} hours (older than {max_file_age_hours} hours)")
                except ValueError as e:
                    print(f"⚠ Could not parse timestamp from {filename}: {e}")
                    # Include file if timestamp parsing fails
                    filtered_files.append(f)
            else:
                print(f"⚠ Skipping file {filename}: Could not extract timestamp from filename")
                # Optionally include files without timestamps
                # filtered_files.append(f)
        
        woest_rhi_files = filtered_files
        print(f"After age filtering: {len(woest_rhi_files)} files remain out of {len(all_rhi_files)} total")
    
    print(f'All RHI files found: {len(all_rhi_files)}')
    print(f'L1 RHI files (matching _l1_v*.*.*.nc): {len(woest_rhi_files)}')
    
    if woest_rhi_files:
        print(f'Sample L1 RHI files to process: {[os.path.basename(f) for f in woest_rhi_files[:3]]}')
    else:
        print('No L1 RHI files found matching the pattern and age criteria')
        return

    for f in woest_rhi_files:
        DS = nc4.Dataset(f, 'r')
        comment_attribute = DS.getncattr('comment') if hasattr(DS, 'comment') else 'No comment'
        
        # Check if antenna_transition exists and validate its values
        if 'antenna_transition' in DS.variables:
            antenna_transition_data = DS.variables['antenna_transition'][:]
            # Skip if all values are 1 or only one value is not 1
            if np.all(antenna_transition_data == 1) or np.sum(antenna_transition_data != 1) == 1:
                print(f"Skipping file {os.path.basename(f)}: 'antenna_transition' is invalid (all 1s or only one value not 1).")
                DS.close()  # Close the dataset
                continue
    
        DS.close()  # Close the dataset
    
        print(f"Processing: {os.path.basename(f)} - {comment_attribute}")
        make_dymecs_rhi_plot(f, figpath, '', blflag=blflag, darkmode=darkmode, azimuth_offset=azimuth_offset, skip_all_transition=skip_all_transition)

    return

def calculate_zoom_level(lat_min, lat_max, lon_min, lon_max, map_width_pixels=1000):
    """
    Calculate appropriate zoom level based on bounding box.
    
    Args:
        lat_min, lat_max, lon_min, lon_max: Bounding box coordinates
        map_width_pixels: Expected map width in pixels
        
    Returns:
        int: Appropriate zoom level (1-20)
    """
    import math
    
    # Calculate bounding box dimensions
    lat_span = lat_max - lat_min
    lon_span = lon_max - lon_min
    
    # Use the larger of the two ranges to determine zoom
    max_span = max(lat_span, lon_span)

    # Rough mapping from degree span to zoom level
    # These values are empirically determined for good tile resolution
    if max_span > 5.0:
        zoom = 8
    elif max_span > 2.0:
        zoom = 9
    elif max_span > 1.0:
        zoom = 10
    elif max_span > 0.5:
        zoom = 11
    elif max_span > 0.25:
        zoom = 12
    elif max_span > 0.125:
        zoom = 13
    elif max_span > 0.0625:
        zoom = 14
    elif max_span > 0.03125:
        zoom = 15
    else:
        zoom = 16
    
    # Cap zoom level to avoid overly detailed tiles
    zoom = min(zoom, 16)
    zoom = max(zoom, 8)

    print(f"Bounding box: {lat_span:.4f}° lat x {lon_span:.4f}° lon")
    print(f"Max range: {max_span:.4f}°, calculated zoom level: {zoom}")

    return zoom

def make_dymecs_ppi_map_plot(radar, os_key, kml_paths, save_path, sweep=0, vmin_vmax=None, 
                           dbz_only=False, basemap_type="opentopo", azimuth_offset=0.0,
                           lat_bounds=None, lon_bounds=None):
    """
    Plots a 2x2 PPI map of radar variables using different basemap backgrounds.
    Now automatically calculates zoom level based on bounding box and includes ray duration in titles.
    """

    ppi_el = radar.get_elevation(0)[0];

    # Calculate average ray duration
    avg_ray_duration = get_average_ray_duration(radar)
    if avg_ray_duration is not None:
        ray_duration_str = f"{avg_ray_duration:.1f}ms"
    else:
        ray_duration_str = "unknown"

    # Setup OpenStreetMap tiles
    tiles = cimgt.GoogleTiles(style="satellite")

    lat_0=radar.latitude["data"][0];
    lon_0=radar.longitude["data"][0];

    projection = ccrs.AzimuthalEquidistant(central_longitude=lon_0, central_latitude=lat_0)

    # Setup PyART gate filtering
    gatefilter = pyart.filters.GateFilter(radar)

    # Radar display setup
    display = pyart.graph.RadarMapDisplay(radar)

    # Define color maps for different radar variables
    colormaps = {
        "DBZ": "HomeyerRainbow",
        "VEL": "balance",
        "WIDTH": "viridis",
        "LDR": "viridis"
    }

    VEL_min = radar.fields['VEL']['field_limit_lower'];
    VEL_max = radar.fields['VEL']['field_limit_upper'];
    
    # Default vmin and vmax values if not provided
    default_vmin_vmax = {
        "DBZ": (-20, 50),
        "VEL": (VEL_min, VEL_max),
        "WIDTH": (1e-1 * np.sqrt(1e-1), np.sqrt(1e1)),  # Updated WIDTH limits
        "LDR": (-35, 5)
    }

    # Use provided bounds or defaults
    if lat_bounds is not None:
        lat_min, lat_max = lat_bounds
    else:
        # Default bounds for KASBEX at CAO
        lat_min, lat_max = 51.25, 51.35
        
    if lon_bounds is not None:
        lon_min, lon_max = lon_bounds
    else:
        # Default bounds for KASBEX at CAO
        lon_min, lon_max = -1.3, -1.15

    print(f"Using map bounds: lat [{lat_min:.3f}, {lat_max:.3f}], lon [{lon_min:.3f}, {lon_max:.3f}]")

    # Calculate appropriate zoom level based on bounding box
    zoom_level = calculate_zoom_level(lat_min, lat_max, lon_min, lon_max)

    if dbz_only:
        # Plot only reflectivity in a single panel
        variables = ["DBZ"]
        positions = [(0, 0)]
        fig, ax = plt.subplots(1, 1, figsize=(15, 15), subplot_kw=dict(projection=projection))
        # Convert single axis to array format for consistent indexing
        ax = np.array([[ax]])
    else:
        # Plot all variables in 2x2 grid
        variables = ["DBZ", "VEL", "WIDTH", "LDR"]
        positions = [(0, 0), (1, 0), (0, 1), (1, 1)]
        fig, ax = plt.subplots(2, 2, figsize=(30, 30), subplot_kw=dict(projection=projection))

        # Adjust spacing between plots
        fig.subplots_adjust(wspace=0.1, hspace=0.1)  # Increase horizontal and vertical spacing

    # Define RHI azimuths for radial lines
    rhi_azimuths = [0.0, 25.0, 44.0, 60.0, 80.0, 99.0, 120.0, 143.0, 168.0, 
                    189.0, 215.0, 232.0, 244.0, 255.0, 270.0, 283.0, 293.0, 
                    310.0, 331.0, 348.0]

    # Maximum radar range in km
    max_range_km = 30.0;
    
    # Get radar location
    radar_lon, radar_lat = radar.longitude['data'][0], radar.latitude['data'][0]

    # Define the geodetic projection (WGS84 Ellipsoid)
    geod = Geod(ellps="WGS84")
    
    # Get sweep start time for titles and filename
    dtime_sweep = cftime.num2pydate(radar.time['data'][0], radar.time['units'])
    dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S")
    
    for var, pos in zip(variables, positions):
        # Use the original colormap
        cmap_orig = plt.get_cmap(colormaps[var])
        
        # Use LogNorm for WIDTH variable
        if var == "WIDTH":
            norm = colors.LogNorm(vmin=default_vmin_vmax[var][0], vmax=default_vmin_vmax[var][1])
        else:
            norm = plt.Normalize(vmin=default_vmin_vmax[var][0], vmax=default_vmin_vmax[var][1])
        
        # Get the standard name (long_name) for the variable
        if var in radar.fields and 'long_name' in radar.fields[var]:
            var_name = radar.fields[var]['long_name']
        else:
            var_name = var  # Fallback to the variable name if long_name is not available
    
        # Plot the PPI map with the specified alpha value
        display.plot_ppi_map(
            var, 0, 
            vmin=default_vmin_vmax[var][0], vmax=default_vmin_vmax[var][1],
            min_lat=40, max_lat=53,  # Keep PyART's internal limits wide
            min_lon=-4, max_lon=2,   # Keep PyART's internal limits wide
            lon_lines=None,
            lat_lines=None,
            projection=projection,
            fig=fig,
            ax=ax[pos],
            lat_0=lat_0,
            lon_0=lon_0,
            cmap=cmap_orig,  # Use original colormap
            norm=norm,  # Apply the appropriate normalization
            alpha=0.5,  # Keep the plot transparent
            colorbar_flag=False,  # Disable automatic colorbar
            gatefilter=gatefilter,
            edges=False,
            embellish=False,
            resolution="10m",
            title=f"{var_name}\nElevation: {ppi_el:.2f}° | Date: {dtime_sweep.strftime('%Y-%m-%d')} | Time: {dtime_sweep.strftime('%H:%M:%S')} | Ray duration: {ray_duration_str}",
        )

        # Create a separate ScalarMappable for the colorbar with alpha=1.0
        colorbar_mappable = ScalarMappable(norm=norm, cmap=cmap_orig)
        colorbar_mappable.set_array([])  # Required for the colorbar to work
        cbar = fig.colorbar(
            colorbar_mappable,
            ax=ax[pos],
            orientation='horizontal',
            shrink=0.6,  # Decrease shrink value to make the colorbar less wide
            pad=0.12  # Adjust padding as needed
        )
        cbar.set_alpha(1.0)  # Ensure the colorbar is fully opaque

        # Add a label to the colorbar with long_name and units
        if var in radar.fields:
            long_name = radar.fields[var].get('long_name', var)  # Fallback to variable name if long_name is missing
            units = radar.fields[var].get('units', '')  # Fallback to empty string if units are missing
            cbar.set_label(f"{long_name} ({units})", fontsize=12)  # Combine long_name and units
        else:
            cbar.set_label(f"{var}", fontsize=12)  # Fallback to variable name if metadata is missing

        display.plot_range_rings([10, 20, 30], '0.5', lw=0.5)


    # Set map extents and add basemap for all panels
    if dbz_only:
        panel_list = [(0, 0)]
    else:
        panel_list = positions
        
    for pos in panel_list:
        # Use the configurable bounds here
        ax[pos].set_extent([lon_min, lon_max, lat_min, lat_max], ccrs.PlateCarree())
        
        # Choose basemap based on type with calculated zoom level
        if basemap_type == "satellite":
            try:
                # Use Google satellite tiles with calculated zoom
                ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.Esri.WorldImagery, crs=projection)
                print(f"Added satellite basemap with zoom level {zoom_level}")
            except Exception as e:
                print(f"Failed to load satellite imagery, falling back to OpenTopoMap: {e}")
                try:
                    ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.OpenTopoMap, crs=projection)
                except Exception as e2:
                    print(f"OpenTopoMap also failed, using lower zoom: {e2}")
                    ctx.add_basemap(ax[pos], zoom=max(1, zoom_level-2), source=ctx.providers.OpenTopoMap, crs=projection)
        elif basemap_type == "cartodb":
            try:
                ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.CartoDB.Positron, crs=projection)
                print(f"Added CartoDB basemap with zoom level {zoom_level}")
            except Exception as e:
                print(f"Failed to load CartoDB, trying lower zoom: {e}")
                ctx.add_basemap(ax[pos], zoom=max(1, zoom_level-2), source=ctx.providers.CartoDB.Positron, crs=projection)
        else:  # Default to opentopo
            try:
                ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.OpenTopoMap, crs=projection)
                print(f"Added OpenTopoMap basemap with zoom level {zoom_level}")
            except Exception as e:
                print(f"Failed to load OpenTopoMap, falling back to CartoDB: {e}")
                try:
                    ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.CartoDB.Positron, crs=projection)
                except Exception as e2:
                    print(f"CartoDB also failed, using lower zoom: {e2}")
                    ctx.add_basemap(ax[pos], zoom=max(1, zoom_level-2), source=ctx.providers.CartoDB.Positron, crs=projection)

        # Add gridlines
        gl = ax[pos].gridlines(draw_labels=False, linestyle="None")
        gridliner = ax[pos].gridlines(draw_labels=True, linestyle="--", linewidth=0.5, color='gray')
        gridliner.xformatter = LONGITUDE_FORMATTER
        gridliner.yformatter = LATITUDE_FORMATTER
        gridliner.xlabel_style = {'size': 12, 'color': 'black'}
        gridliner.ylabel_style = {'size': 12, 'color': 'black'}

    # Save with appropriate filename
    if dbz_only:
        figname = f'ncas-mobile-ka-band-radar-1_cao_{dtime_sweep_str}_ppi_el{ppi_el:0.2f}_dbz_map_l1_v1.0.0.png';
    else:
        figname = f'ncas-mobile-ka-band-radar-1_cao_{dtime_sweep_str}_ppi_el{ppi_el:0.2f}_map_l1_v1.0.0.png';

    plt.savefig(os.path.join(save_path,figname),dpi=300, bbox_inches='tight');
    plt.close();

def make_dymecs_ppi_map_plots_day(datestr, inpath, figpath, dbz_only=False, basemap_type="opentopo", 
                                azimuth_offset=0.0, lat_bounds=None, lon_bounds=None, max_file_age_hours=None):
    inpath_date = inpath
    
    os.chdir(inpath_date)
    print(inpath_date)
    
    import re
    
    # Find all potential PPI files
    all_ppi_files = glob.glob('*{}*_ppi*.nc'.format(datestr))
    
    # Filter for L1 data version pattern: _l1_v*.*.*.nc (excludes _old and numbered files)
    l1_pattern = re.compile(r'_l1_v\d+\.\d+\.\d+\.nc$')
    dymecs_ppi_files = [f for f in all_ppi_files if l1_pattern.search(f)]
    
    # Apply scan age filter
    dymecs_ppi_files = filter_files_by_scan_age(dymecs_ppi_files, max_file_age_hours)
    
    print(f'All PPI files found: {len(all_ppi_files)}')
    print(f'L1 PPI files (matching _l1_v*.*.*.nc): {len(dymecs_ppi_files)}')
    print(f'ppi files = {all_ppi_files}')
    print(f'filtered L1 ppi files = {dymecs_ppi_files}')
    
    if dymecs_ppi_files:
        print(f'Sample L1 PPI files: {[f for f in dymecs_ppi_files[:3]]}')
    else:
        print('No L1 PPI files found matching the pattern and age criteria')
        return
    
    # Create output directory for PPI maps
    save_path = os.path.join(figpath, 'ppi', datestr)
    if not os.path.isdir(save_path): 
        os.makedirs(save_path)
    
    for f in dymecs_ppi_files:
        print(f"Processing PPI file: {f}")
        radar = pyart.io.read(os.path.join(inpath_date, f))
        
        # Apply azimuth offset if specified
        if azimuth_offset != 0.0:
            radar.azimuth['data'] += azimuth_offset
            # Ensure azimuths stay within 0-360 range
            radar.azimuth['data'] = radar.azimuth['data'] % 360.0
        
        for sweep in range(radar.nsweeps):
            radar_sweep = radar.extract_sweeps([sweep])
            make_dymecs_ppi_map_plot(radar_sweep, os_key, kml_paths, save_path, 
                                   sweep=sweep, dbz_only=dbz_only, basemap_type=basemap_type, 
                                   azimuth_offset=azimuth_offset, lat_bounds=lat_bounds, lon_bounds=lon_bounds)

def make_cobalt_vpt_plot_day(datestr, inpath, figpath, blflag=False, max_file_age_hours=None):
    """
    Create VPT (Vertical Pointing) quicklook plot for a full day.
    Now filters for L1 data version files only and includes scan age filtering.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path (already includes date directory)
        figpath: Output directory path
        blflag: Boundary layer flag for plot limits
        max_file_age_hours: Maximum scan age in hours, None for no limit
    """
    print(f"Creating VPT plot for date: {datestr}")
    
    # Set height limits
    hmax = 4 if blflag else 12
    
    # Find VPT file - inpath already includes the date directory
    inpath_date = inpath  # Don't add datestr again!
    
    # Check if directory exists
    if not os.path.exists(inpath_date):
        print(f"Input directory does not exist: {inpath_date}")
        return
        
    os.chdir(inpath_date)
    print(f"Looking for VPT files in: {inpath_date}")
    
    import re
    
    try:
        # Find all potential VPT files
        all_vpt_files = glob.glob(f'*{datestr}*vpt*.nc')
        
        # Filter for L1 data version pattern: _l1_v*.*.*.nc
        l1_pattern = re.compile(r'_l1_v\d+\.\d+\.\d+\.nc$')
        vpt_files = [f for f in all_vpt_files if l1_pattern.search(f)]
        
        # Apply scan age filter
        vpt_files_with_paths = [os.path.join(inpath_date, f) for f in vpt_files]
        vpt_files_filtered = filter_files_by_scan_age(vpt_files_with_paths, max_file_age_hours)
        vpt_files = [os.path.basename(f) for f in vpt_files_filtered]
        
        print(f'All VPT files found: {len(all_vpt_files)}')
        print(f'L1 VPT files (matching _l1_v*.*.*.nc): {len(vpt_files)}')
        
        if not vpt_files:
            print(f"No L1 VPT files found for {datestr} matching scan age criteria")
            if all_vpt_files:
                print(f"Available VPT files (not L1 or scans too old): {all_vpt_files}")
            return
            
        vpt_file = os.path.join(inpath_date, vpt_files[0])
        print(f"Processing VPT file: {vpt_file}")
        
    except Exception as e:
        print(f"Error finding VPT file: {e}")
        return
    
    # Read data
    try:
        with nc4.Dataset(vpt_file) as ds:
            product_version = ds.product_version if hasattr(ds, 'product_version') else 'v1.0.0'
            
        radar_ds = pyart.io.read_cfradial(vpt_file)
        
    except Exception as e:
        print(f"Error reading VPT file: {e}")
        return
    
    # Get velocity limits
    vel_field = radar_ds.fields['VEL']
    vel_limit_lower = vel_field.get('field_limit_lower', -20)
    vel_limit_upper = vel_field.get('field_limit_upper', 20)
    
    # Try to load previous day's data for continuity
    nsweeps_prev = 0
    radar_ds_prev = None
    
    try:
        prev_date = datetime.datetime.strptime(datestr, '%Y%m%d') - datetime.timedelta(days=1)
        prevstr = prev_date.strftime('%Y%m%d')
        # Construct previous day path properly
        base_path = os.path.dirname(inpath_date)  # Remove current date from path
        inpath_prev = os.path.join(base_path, prevstr)
        
        if os.path.exists(inpath_prev):
            os.chdir(inpath_prev)
            
            # Find L1 VPT files from previous day
            all_vpt_files_prev = glob.glob(f'*{prevstr}*vpt*.nc')
            vpt_files_prev = [f for f in all_vpt_files_prev if l1_pattern.search(f)]
            
            # Apply scan age filter to previous day files too
            vpt_files_prev_with_paths = [os.path.join(inpath_prev, f) for f in vpt_files_prev]
            vpt_files_prev_filtered = filter_files_by_scan_age(vpt_files_prev_with_paths, max_file_age_hours)
            vpt_files_prev = [os.path.basename(f) for f in vpt_files_prev_filtered]
            
            if vpt_files_prev:
                vpt_file_prev = os.path.join(inpath_prev, vpt_files_prev[0])
                radar_ds_prev = pyart.io.read_cfradial(vpt_file_prev)
                nsweeps_prev = radar_ds_prev.nsweeps
                print(f"Loaded previous day L1 data: {nsweeps_prev} sweeps")
        else:
            print(f"Previous day directory does not exist: {inpath_prev}")
            
    except Exception as e:
        print(f"Could not load previous day data: {e}")
    
    # Create figure
    fig, axes = plt.subplots(4, 1, figsize=(12, 18), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=1/72, h_pad=1/72, hspace=0.2, wspace=0.2)
    
    # Set up time limits for full day
    dtime = cftime.num2pydate(radar_ds.time['data'], radar_ds.time['units'])
    dt_min = dtime[0].replace(hour=0, minute=0, second=0)
    dt_max = dt_min + datetime.timedelta(days=1)
    time_str = dtime[0].strftime("%Y-%m-%d")
    
    # Plot first sweep to establish colorbars and titles
    radar_sweep_ds = radar_ds.extract_sweeps([0])
    gatefilter = pyart.correct.GateFilter(radar_sweep_ds)
    #gatefilter.exclude_below('SNR', -20)
    display = pyart.graph.RadarDisplay(radar_sweep_ds)
    
    _plot_vpt_fields(display, axes, radar_sweep_ds, gatefilter, time_str, hmax,
                    vel_limit_lower, vel_limit_upper)
    
    # Plot remaining sweeps from current day
    nsweeps = radar_ds.nsweeps
    for sweep_idx in range(1, nsweeps):
        print(f"Processing sweep {sweep_idx}/{nsweeps}")
        
        radar_sweep_ds = radar_ds.extract_sweeps([sweep_idx])
        gatefilter = pyart.correct.GateFilter(radar_sweep_ds)
        gatefilter.exclude_below('SNR', -20)
        display = pyart.graph.RadarDisplay(radar_sweep_ds)
        
        _plot_vpt_fields_overlay(display, axes, gatefilter, vel_limit_lower, vel_limit_upper)
    
    # Plot last sweep from previous day if available
    if radar_ds_prev and nsweeps_prev > 0:
        print("Adding previous day data")
        radar_sweep_ds = radar_ds_prev.extract_sweeps([nsweeps_prev - 1])
        gatefilter = pyart.correct.GateFilter(radar_sweep_ds)
        #gatefilter.exclude_below('SNR', -20)
        display = pyart.graph.RadarDisplay(radar_sweep_ds)
        
        _plot_vpt_fields_overlay(display, axes, gatefilter, vel_limit_lower, vel_limit_upper)
    
    # Final plot setup
    for ax in axes:
        ax.set_xlim(dt_min, dt_max)
        ax.grid(True)
        ax.set_xlabel('Time (UTC)')
    
    # Save figure
    figname = f'ncas-mobile-ka-band-radar-1_cao_{datestr}_vpt_l1_{product_version}.png'
    if blflag:
        figname = figname.replace('.png', '_bl.png')
    
    vpt_figpath = os.path.join(figpath, 'vpt')
    os.makedirs(vpt_figpath, exist_ok=True)
    
    plt.savefig(os.path.join(vpt_figpath, figname), dpi=300)
    plt.close()
    
    print(f"VPT plot saved to: {vpt_figpath}")

def _plot_vpt_fields(display, axes, radar_ds, gatefilter, time_str, hmax, 
                    vel_min, vel_max):
    """Plot VPT fields with titles and colorbars including ray duration."""
    
    # Calculate average ray duration for VPT plots
    avg_ray_duration = get_average_ray_duration(radar_ds)
    if avg_ray_duration is not None:
        ray_duration_str = f"{avg_ray_duration:.1f}ms"
    else:
        ray_duration_str = "unknown"
    
    # DBZ
    field_name = pyart.graph.common.generate_field_name(radar_ds, "DBZ")
    title = f"{pyart.graph.common.generate_radar_name(radar_ds)} {time_str}\n{field_name} | Date: {dtime[0].strftime('%Y-%m-%d')} | Time: {dtime[0].strftime('%H:%M:%S')} | Ray duration: {ray_duration_str}"
    display.plot_vpt("DBZ", ax=axes[0], time_axis_flag=True, title=title, edges=False,
                     gatefilter=gatefilter, vmin=-60, vmax=40, cmap=COLORMAPS['dbz'],
                     colorbar_orient='horizontal', filter_transitions=True)
    axes[0].set_ylim(0, hmax)
    
    # VEL  
    field_name = pyart.graph.common.generate_field_name(radar_ds, "VEL")
    title = f"{pyart.graph.common.generate_radar_name(radar_ds)} {time_str}\n{field_name} ray_dur = {ray_duration_str}"
    display.plot_vpt("VEL", ax=axes[1], time_axis_flag=True, title=title, edges=False,
                     gatefilter=gatefilter, vmin=vel_min, vmax=vel_max, 
                     cmap=COLORMAPS['vel'], colorbar_orient='horizontal')
    axes[1].set_ylim(0, hmax)
    
    # WIDTH
    field_name = pyart.graph.common.generate_field_name(radar_ds, "WIDTH")
    title = f"{pyart.graph.common.generate_radar_name(radar_ds)} {time_str}\n{field_name} ray_dur = {ray_duration_str}"
    display.plot_vpt("WIDTH", ax=axes[2], time_axis_flag=True, title=title, edges=False,
                     gatefilter=gatefilter, norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                     cmap=COLORMAPS['width'], colorbar_orient='horizontal')
    axes[2].set_ylim(0, hmax)
    
    # LDR
    field_name = pyart.graph.common.generate_field_name(radar_ds, "LDR")
    title = f"{pyart.graph.common.generate_radar_name(radar_ds)} {time_str}\n{field_name} ray_dur = {ray_duration_str}"
    display.plot_vpt("LDR", ax=axes[3], time_axis_flag=True, title=title, edges=False,
                     gatefilter=gatefilter, vmin=-35, vmax=5, cmap=COLORMAPS['ldr'],
                     colorbar_orient='horizontal')
    axes[3].set_ylim(0, hmax)

def _plot_vpt_fields_overlay(display, axes, gatefilter, vel_min, vel_max):
    """Plot VPT fields as overlay (no titles or colorbars)."""
    
    display.plot_vpt("DBZ", ax=axes[0], time_axis_flag=True, edges=False,
                     gatefilter=gatefilter, vmin=-60, vmax=40, cmap=COLORMAPS['dbz'],
                     colorbar_flag=False, title_flag=False, filter_transitions=True)
    
    display.plot_vpt("VEL", ax=axes[1], time_axis_flag=True, edges=False,
                     gatefilter=gatefilter, vmin=vel_min, vmax=vel_max, 
                     cmap=COLORMAPS['vel'], colorbar_flag=False, title_flag=False)
    
    display.plot_vpt("WIDTH", ax=axes[2], time_axis_flag=True, edges=False,
                     gatefilter=gatefilter, norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                     cmap=COLORMAPS['width'], colorbar_flag=False, title_flag=False)
    
    display.plot_vpt("LDR", ax=axes[3], time_axis_flag=True, edges=False,
                     gatefilter=gatefilter, vmin=-35, vmax=5, cmap=COLORMAPS['ldr'],
                     colorbar_flag=False, title_flag=False)
    
# Enable the function calls at the bottom
basemap_type="cartodb"  # Options: "opentopo", "satellite", "cartodb"
make_dymecs_ppi_map_plots_day(datestr, inpath, figpath, dbz_only=dbz_only, basemap_type=basemap_type,
                             azimuth_offset=azimuth_offset, lat_bounds=(lat_min, lat_max), lon_bounds=(lon_min, lon_max),
                             max_file_age_hours=max_file_age_hours);
make_dymecs_rhi_plots_day(datestr,inpath,figpath,blflag=blflag,darkmode=darkmode,azimuth_offset=azimuth_offset,
                         max_file_age_hours=max_file_age_hours,skip_all_transition=skip_all_transition);
#make_cobalt_vpt_plot_day(datestr, inpath, figpath, blflag=blflag, max_file_age_hours=max_file_age_hours)

