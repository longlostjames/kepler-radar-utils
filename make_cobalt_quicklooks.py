#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
make_cobalt_quicklooks.py

Generate quicklook plots for COBALT campaign radar data.

This script creates visualization plots for NCAS Mobile Ka-band Radar data
from the COBALT (Contrail Observations and Lifecycle Tracking) campaign.
It generates RHI (Range Height Indicator) and VPT (Vertical Pointing) plots
with overlaid aircraft position information.

Usage:
    python make_cobalt_quicklooks.py -d YYYYMMDD -i input_path -o output_path [-b] [-l] [--no-aircraft] [--skip-all-transition]

Arguments:
    -d, --date:              Date string in YYYYMMDD format
    -i, --inpath:            Input directory containing CF-Radial files
    -o, --outpath:           Output directory for quicklook images
    -b:                      Boundary layer flag (limits plots to 4km height)
    -l:                      Latest flag (process most recent data)
    --no-aircraft:           Skip aircraft markers and related processing
    --skip-all-transition:   Skip sweeps where 100% of rays are antenna_transition=1

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
Last modified: 10-08-2025
Version: 0.1
"""

import getopt
import sys
import os
import datetime
import glob
import re
from pathlib import Path

import netCDF4 as nc4
import pyart
import numpy as np
import numpy.ma as ma
import cftime
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
import cmocean

# Import kepler utilities
from kepler_utils import get_valid_sweep_indices

# Configuration
VERSION = 0.1
TRACKING_TAG = 'AMOF_20231120125118'
CAMPAIGN = 'cobalt'

# Default paths
NCAS_RADAR_PATH = '/gws/nopw/j04/ncas_radar_vol2'
NCAS_OBS_PROC_PATH = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1'
NCAS_OBS_RAW_PATH = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-mobile-ka-band-radar-1'

# Plot configuration
COLORMAPS = {
    'dbz': 'HomeyerRainbow',
    'vel': 'balance', 
    'ldr': 'viridis',
    'width': 'SpectralExtended'
}

# Earth radius in km for distance calculations
EARTH_RADIUS_KM = 6371.0

def parse_command_line():
    """
    Parse command line arguments.
    
    Returns:
        tuple: (datestr, inpath, outpath, blflag, latest, no_aircraft)
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:i:o:bl", 
                                   ["date=", "inpath=", "outpath=", "no-aircraft", "skip-all-transition"])
    except getopt.GetoptError as err:
        print(f"Error: {err}")
        print("Usage: python make_cobalt_quicklooks.py -d YYYYMMDD -i input_path -o output_path [-b] [-l] [--no-aircraft] [--skip-all-transition]")
        sys.exit(2)

    # Default values
    data_date = datetime.datetime.now()
    datestr = data_date.strftime('%Y%m%d')
    inpath = None
    outpath = None
    blflag = False
    latest = False
    no_aircraft = False
    skip_all_transition = False

    for option, argument in opts:
        if option in ("-d", "--date"):
            datestr = argument
        elif option in ("-i", "--inpath"):
            inpath = argument
        elif option in ("-o", "--outpath"):
            outpath = argument
        elif option == "-b":
            blflag = True
        elif option == "-l":
            latest = True
        elif option == "--no-aircraft":
            no_aircraft = True
        elif option == "--skip-all-transition":
            skip_all_transition = True
        else:
            print(f"Unhandled option: {option}")
            sys.exit(2)

    return datestr, inpath, outpath, blflag, latest, no_aircraft, skip_all_transition

def setup_paths(datestr):
    """
    Set up input and figure output paths.
    
    Args:
        datestr: Date string in YYYYMMDD format
        
    Returns:
        tuple: (inpath, figpath, cobalt_cmd_path)
    """
    homepath = Path.home()
    
    # Input path for processed radar data
    #inpath = os.path.join(NCAS_OBS_PROC_PATH, CAMPAIGN, 'L1_v1.0.1')
    inpath = os.path.join(NCAS_OBS_PROC_PATH, CAMPAIGN, 'L1c')

    # Output path for figures
    figpath = os.path.join(inpath, 'quicklooks')
    
    # Path to COBALT command files (aircraft tracking info)
    cobalt_cmd_path = os.path.join(NCAS_OBS_RAW_PATH, 'data', 'campaign', CAMPAIGN, 'cobalt-command')

    return inpath, figpath, cobalt_cmd_path

def read_cobalt_cmd_file(file_path):
    """
    Read COBALT command file containing aircraft information.
    
    Args:
        file_path: Path to command file
        
    Returns:
        dict: Dictionary of command data (flight_id, icao_id, position, etc.)
    """
    try:
        with open(file_path, "r") as file:
            line = file.readline()
            # Clean up the line format
            newline = line.rstrip().replace(': ', ':')
            pairs = newline.split()
            data = {}

            for pair in pairs:
                if ':' in pair:
                    key, value = pair.split(':', 1)
                    data[key.strip()] = value.strip()

        return data
    except Exception as e:
        print(f"Error reading command file {file_path}: {e}")
        return {}

def extract_azimuth_from_file(file_path):
    """
    Extract azimuth angle from COBALT command file.
    
    Args:
        file_path: Path to command file
        
    Returns:
        float: Azimuth angle in degrees
    """
    data = read_cobalt_cmd_file(file_path)
    return float(data.get('azimuth', 0))

def find_nearest_file(directory, netcdf_time, target_azimuth):
    """
    Find the nearest COBALT command file to a given time and azimuth.
    
    Args:
        directory: Directory containing command files
        netcdf_time: Target datetime
        target_azimuth: Target azimuth angle in degrees
        
    Returns:
        tuple: (nearest_file, file_time) or (None, None) if not found
    """
    log_file_pattern = re.compile(r"kepler-(\d{8})-(\d{6})\.log")
    
    nearest_file = None
    nearest_file_time = None
    nearest_azimuth_diff = float('inf')
    candidates = []

    # First pass: find files with the closest azimuth within 10 minutes
    for file_name in os.listdir(directory):
        match = log_file_pattern.match(file_name)
        if match:
            # Extract datetime from filename
            date_part = match.group(1)
            time_part = match.group(2)
            file_time_str = f"{date_part}{time_part}"
            file_time = datetime.datetime.strptime(file_time_str, "%Y%m%d%H%M%S")

            # Only consider files within 10 minutes
            time_diff = abs(netcdf_time - file_time)
            if time_diff < datetime.timedelta(seconds=600):
                # Extract azimuth and find closest match
                azimuth = extract_azimuth_from_file(os.path.join(directory, file_name))
                azimuth_diff = abs(target_azimuth - azimuth)

                if azimuth_diff < nearest_azimuth_diff:
                    nearest_azimuth_diff = azimuth_diff
                    candidates = [(file_name, file_time, time_diff)]
                elif azimuth_diff == nearest_azimuth_diff:
                    candidates.append((file_name, file_time, time_diff))

    # Second pass: among azimuth candidates, find closest time
    nearest_time_diff = datetime.timedelta.max
    for file_name, file_time, time_diff in candidates:
        if time_diff < nearest_time_diff:
            nearest_file = file_name
            nearest_file_time = file_time
            nearest_time_diff = time_diff

    return nearest_file, nearest_file_time

def get_time_coverage_start(netcdf_file_path):
    """
    Extract time_coverage_start from NetCDF file.
    
    Args:
        netcdf_file_path: Path to NetCDF file
        
    Returns:
        str: Time coverage start string or None if not found
    """
    try:
        with nc4.Dataset(netcdf_file_path, 'r') as ds:
            if 'time_coverage_start' in ds.ncattrs():
                return ds.getncattr('time_coverage_start')
            else:
                raise ValueError("Attribute 'time_coverage_start' not found")
    except Exception as e:
        print(f"Error reading NetCDF file: {e}")
        return None

def great_circle_distance(earth_radius, target_range, zenith_angle):
    """
    Calculate great-circle distance from radar to target.
    
    Args:
        earth_radius: Earth radius in km
        target_range: Range to target in km
        zenith_angle: Zenith angle in radians
        
    Returns:
        float: Great-circle distance in km
    """
    target_altitude = get_target_altitude(earth_radius, target_range, zenith_angle)
    
    # Calculate angular distance
    delta = np.arcsin(target_range * np.sin(zenith_angle) / (earth_radius + target_altitude))
    
    # Convert to great-circle distance
    great_circle_dist = earth_radius * delta
    return great_circle_dist

def get_target_altitude(earth_radius, target_range, zenith_angle):
    """
    Calculate altitude of radar target above Earth's surface.
    
    Args:
        earth_radius: Earth radius in km
        target_range: Range to target in km
        zenith_angle: Zenith angle in radians
        
    Returns:
        float: Target altitude in km
    """
    # Law of cosines to get distance from Earth center to target
    r_t = np.sqrt(
        earth_radius**2 + target_range**2 + 
        2 * earth_radius * target_range * np.cos(zenith_angle)
    )
    
    # Altitude is distance from Earth surface
    altitude = np.abs(r_t - earth_radius)
    return altitude

def setup_plot_limits(blflag, wind_scan=False):
    """
    Set up plot limits based on scan type and boundary layer flag.
    
    Args:
        blflag: Boundary layer flag (limits height to 4km)
        wind_scan: Whether this is a wind scan (wider elevation range)
        
    Returns:
        tuple: (hmax, xmin, xmax) - height max and horizontal limits
    """
    if blflag:
        return 4, 0, 20
    elif wind_scan:
        return 14, -14, 14
    else:
        return 14, 0, 40

def make_cobalt_rhi_plot(ncfile, figpath, cobalt_cmd_path, blflag=False, no_aircraft=False, skip_all_transition=False):
    """
    Create RHI (Range Height Indicator) quicklook plot for COBALT data.
    
    Args:
        ncfile: Path to CF-Radial NetCDF file
        figpath: Output directory for figure
        cobalt_cmd_path: Path to COBALT command files
        blflag: Boundary layer flag for plot limits
        no_aircraft: If True, skip aircraft markers and related processing
        skip_all_transition: If True, skip sweeps where 100% of rays are antenna_transition=1
    """
    print(f"Creating RHI plot for: {ncfile}")
    
    # Get scan information
    time_coverage_start = get_time_coverage_start(ncfile)
    if not time_coverage_start:
        print(f"Could not get time coverage from {ncfile}")
        return
        
    netcdf_time = datetime.datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%SZ")
    
    with nc4.Dataset(ncfile) as ds:
        scan_azimuth = ds['fixed_angle'][0]
        product_version = ds.product_version
        elevation_min = np.min(ds['elevation'][:])
        elevation_max = np.max(ds['elevation'][:])
        scan_span = np.abs(elevation_max - elevation_min)
    
    print(f"Scan azimuth: {scan_azimuth:.1f}deg, elevation range: {elevation_min:.1f}° to {elevation_max:.1f}° (span: {scan_span:.1f}°)")
    
    # Determine if this is a wind scan (elevation passes through 90 degrees AND span > 60 degrees)
    passes_through_zenith = elevation_min <= 90.0 <= elevation_max
    large_span = scan_span > 60
    wind_scan = passes_through_zenith and large_span
    
    if wind_scan:
        print("Detected as wind scan (elevation passes through 90° AND span > 60°)")
    elif passes_through_zenith and not large_span:
        print(f"Detected as regular RHI scan (passes through 90° but span only {scan_span:.1f}°)")
    elif large_span and not passes_through_zenith:
        print(f"Detected as regular RHI scan (large span {scan_span:.1f}° but doesn't pass through 90°)")
    else:
        print("Detected as regular RHI scan")
    
    # Set up plot limits
    hmax, xmin, xmax = setup_plot_limits(blflag, wind_scan)
    
    # Find corresponding aircraft information if not a wind scan and aircraft plotting is enabled
    cmd_data = {}
    logfile_time_str = "N/A"
    
    if not wind_scan and not no_aircraft:
        nearest_file, logfile_time = find_nearest_file(cobalt_cmd_path, netcdf_time, scan_azimuth)
        
        if nearest_file and logfile_time:
            logfile_time_str = logfile_time.strftime('%Y-%m-%dT%H:%M:%SZ')
            cmd_data = read_cobalt_cmd_file(os.path.join(cobalt_cmd_path, nearest_file))
            print(f"Found matching command file: {nearest_file}")
        else:
            print("No matching command file found")
    elif no_aircraft:
        print("Aircraft marker plotting disabled")
    elif wind_scan:
        print("Skipping aircraft data lookup for wind scan")
    
    # Read radar data
    try:
        radar_ds = pyart.io.read_cfradial(ncfile)
    except Exception as e:
        print(f"Error reading radar file: {e}")
        return
    
    # Set up gate filter
    gatefilter = pyart.correct.GateFilter(radar_ds)
    gatefilter.exclude_below('SNR', -20)
    
    # Create display object
    display = pyart.graph.RadarDisplay(radar_ds)
    
    # Set up output directory
    rhi_figpath = os.path.join(figpath, 'rhi', netcdf_time.strftime('%Y%m%d'))
    os.makedirs(rhi_figpath, exist_ok=True)
    
    # Get velocity limits from data
    vel_field = radar_ds.fields['VEL']
    vel_limit_lower = vel_field['field_limit_lower']
    vel_limit_upper = vel_field['field_limit_upper']
    
    # Create plots for each sweep
    nsweeps = radar_ds.nsweeps
    dtime0 = cftime.num2pydate(radar_ds.time['data'][0], radar_ds.time['units'])
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S")
    
    # Get valid sweep indices (optionally filtering out all-transition sweeps)
    valid_sweep_indices = get_valid_sweep_indices(radar_ds, skip_all_transition=skip_all_transition)
    
    if skip_all_transition and len(valid_sweep_indices) < nsweeps:
        print(f"Skipping {nsweeps - len(valid_sweep_indices)} sweep(s) that are 100% antenna transitions")
    
    if nsweeps > 1:
        # Multiple sweeps - create separate plots
        for sweep_idx in valid_sweep_indices:
            print(f"Processing sweep {sweep_idx}/{nsweeps}")
            
            fig, axes = plt.subplots(2, 2, figsize=(15, 15), constrained_layout=False)
            
            # Reserve space for logos BEFORE plotting
            fig.subplots_adjust(left=0.08, right=0.95, top=0.88, bottom=0.08, 
                              wspace=0.3, hspace=0.3)
            
            rhi_az = radar_ds.get_azimuth(sweep_idx)[0]
            
            # Plot each field
            _plot_rhi_fields(display, axes, sweep_idx, gatefilter, vel_limit_lower, 
                           vel_limit_upper, hmax, xmin, xmax)
            
            # Add aircraft markers if enabled
            if cmd_data and not no_aircraft and not wind_scan:
                _add_aircraft_markers(axes, cmd_data)
            
            # Update title with aircraft information
            _update_plot_title(fig, axes, cmd_data, logfile_time_str, wind_scan, no_aircraft)
            
            # Add logos to figure (in reserved space)
            add_logos_to_figure(fig)
            
            # Adjust layout with padding for logos before adding them
            plt.tight_layout(pad=3.0, rect=[0.08, 0.08, 0.95, 0.90])
            
            # Save individual sweep plot
            sweep_start_index = radar_ds.get_start(sweep_idx)
            dtime_sweep = cftime.num2pydate(radar_ds.time['data'][sweep_start_index], 
                                          radar_ds.time['units'])
            dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S")
            
            # Format azimuth to 2 decimal places
            figname = (f'ncas-radar-mobile-ka-band-1_cao_cobalt_{dtime_sweep_str}'
                      f'_rhi_az{rhi_az:.2f}_l1_{product_version}.png')
            
            plt.savefig(os.path.join(rhi_figpath, figname), dpi=300)
            plt.close()
            
    else:
        # Single sweep - create combined plot
        fig, axes = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=False)
        
        # Reserve space for logos BEFORE plotting
        fig.subplots_adjust(left=0.08, right=0.95, top=0.88, bottom=0.08, 
                          wspace=0.3, hspace=0.3)
        
        # Get azimuth for filename
        rhi_az = radar_ds.get_azimuth(0)[0]
        
        # Plot fields
        _plot_rhi_fields(display, axes, 0, gatefilter, vel_limit_lower, 
                        vel_limit_upper, hmax, xmin, xmax)
        
        # Add aircraft position markers if available and enabled
        if cmd_data and not wind_scan and not no_aircraft:
            _add_aircraft_markers(axes, cmd_data)
        
        # Update title with aircraft information
        _update_plot_title(fig, axes, cmd_data, logfile_time_str, wind_scan, no_aircraft)
        
        # Add logos to figure (in reserved space)
        add_logos_to_figure(fig)
        
        # Save plot with azimuth formatted to 2 decimal places
        figname = f'ncas-mobile-ka-band-radar-1_cao_cobalt_{dtime0_str}_rhi_az{rhi_az:.2f}_l1_{product_version}.png'
        plt.savefig(os.path.join(rhi_figpath, figname), dpi=300)
        plt.close()
    
    print(f"RHI plot saved to: {rhi_figpath}")

def _plot_rhi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max, 
                    hmax, xmin, xmax):
    """Helper function to plot RHI fields on axes."""
    
    # Reflectivity (DBZ)
    display.plot_rhi("DBZ", ax=axes[0,0], sweep=sweep_idx, vmin=-60, vmax=40,
                     cmap=COLORMAPS['dbz'], colorbar_orient='horizontal',
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[0,0], hmax, xmin, xmax)
    
    # Velocity (VEL)  
    display.plot_rhi("VEL", ax=axes[1,0], sweep=sweep_idx, vmin=vel_min, vmax=vel_max,
                     cmap=COLORMAPS['vel'], colorbar_orient='horizontal',
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[1,0], hmax, xmin, xmax)
    
    # Spectrum width (WIDTH)
    display.plot_rhi("WIDTH", ax=axes[1,1], sweep=sweep_idx,
                     norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                     cmap=COLORMAPS['width'], colorbar_orient='horizontal',
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[1,1], hmax, xmin, xmax)
    
    # Linear depolarization ratio (LDR)
    display.plot_rhi("LDR", ax=axes[0,1], sweep=sweep_idx, vmin=-35, vmax=5,
                     cmap=COLORMAPS['ldr'], colorbar_orient='horizontal',
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[0,1], hmax, xmin, xmax)

def _setup_rhi_axes(ax, hmax, xmin, xmax):
    """Helper function to set up RHI plot axes."""
    ax.set_ylim(0, hmax)
    ax.set_xlim(xmin, xmax)
    ax.grid(True)
    ax.set_aspect('equal')

def _add_aircraft_markers(axes, cmd_data):
    """Add aircraft position markers to RHI plots."""
    try:
        target_range = float(cmd_data['range']) / 1000.0  # Convert to km
        target_elev = np.deg2rad(float(cmd_data['elevation']))
        
        # Calculate positions
        target_height = target_range * np.sin(target_elev)
        target_zenith_angle = np.pi/2.0 - target_elev
        target_horiz_dist = target_range * np.cos(target_elev)
        
        # Great circle distance and altitude
        target_arc_dist = great_circle_distance(EARTH_RADIUS_KM, target_range, target_zenith_angle)
        target_altitude = get_target_altitude(EARTH_RADIUS_KM, target_range, np.pi - target_zenith_angle)
        
        print(f'Target position - Range: {target_range:.1f}km, Height: {target_height:.1f}km')
        
        # Add markers to all subplots
        for ax in axes.flat:
            ax.plot(target_arc_dist, target_altitude, 'k+', markersize=8)  # Black cross
            ax.plot(target_horiz_dist, target_height, 'b+', markersize=8)  # Blue cross
            
    except (KeyError, ValueError) as e:
        print(f"Could not add aircraft markers: {e}")

def _update_plot_title(fig, axes, cmd_data, logfile_time_str, wind_scan, no_aircraft=False):
    """Update plot title with aircraft information."""
    
    # Get original title from first subplot
    orig_title = axes[0,0].get_title()
    title_lines = orig_title.split("\n")
    
    # Clear individual subplot titles
    for ax in axes.flat:
        ax.set_title("", fontsize=8)
    
    if cmd_data and not wind_scan and not no_aircraft:
        # Create aircraft information string
        aircraft_info = (f"{logfile_time_str} {cmd_data.get('flight_id', 'N/A')} "
                        f"{cmd_data.get('icao_id', 'N/A')} {cmd_data.get('aircraft_type', 'N/A')} "
                        f"lean_burn={cmd_data.get('lean_burn', 'N/A')} "
                        f"Az={float(cmd_data.get('azimuth', 0)):.1f}deg "
                        f"El={float(cmd_data.get('elevation', 0)):.1f}deg")
        
        if 'advected' in cmd_data:
            aircraft_info += f" advected {cmd_data['advected']}"
            
        title_lines[1] = aircraft_info
    else:
        if wind_scan:
            title_lines[1] = "Wind Scan"
        elif no_aircraft:
            title_lines[1] = "Aircraft markers disabled"
        else:
            title_lines[1] = "No aircraft data"
    
    # Set figure title
    new_title = "\n".join(title_lines)
    fig.suptitle(new_title, fontsize=11)

def make_cobalt_vpt_plot_day(datestr, inpath, figpath, blflag=False):
    """
    Create VPT (Vertical Pointing) quicklook plot for a full day.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Output directory path
        blflag: Boundary layer flag for plot limits
    """
    print(f"Creating VPT plot for date: {datestr}")
    
    # Set height limits
    hmax = 4 if blflag else 12
    
    # Find VPT file
    inpath_date = os.path.join(inpath, datestr)
    os.chdir(inpath_date)
    
    try:
        vpt_files = glob.glob(f'*{datestr}*vpt*.nc')
        if not vpt_files:
            print(f"No VPT files found for {datestr}")
            return
            
        vpt_file = os.path.join(inpath_date, vpt_files[0])
        print(f"Processing VPT file: {vpt_file}")
        
    except Exception as e:
        print(f"Error finding VPT file: {e}")
        return
    
    # Read data
    try:
        with nc4.Dataset(vpt_file) as ds:
            product_version = ds.product_version
            
        radar_ds = pyart.io.read_cfradial(vpt_file)
        
    except Exception as e:
        print(f"Error reading VPT file: {e}")
        return
    
    # Get velocity limits
    vel_field = radar_ds.fields['VEL']
    vel_limit_lower = vel_field['field_limit_lower'] 
    vel_limit_upper = vel_field['field_limit_upper']
    
    # Try to load previous day's data for continuity
    nsweeps_prev = 0
    radar_ds_prev = None
    
    try:
        prev_date = datetime.datetime.strptime(datestr, '%Y%m%d') - datetime.timedelta(days=1)
        prevstr = prev_date.strftime('%Y%m%d')
        inpath_prev = os.path.join(inpath, prevstr)
        
        os.chdir(inpath_prev)
        vpt_files_prev = glob.glob(f'*{prevstr}*vpt*.nc')
        
        if vpt_files_prev:
            vpt_file_prev = os.path.join(inpath_prev, vpt_files_prev[0])
            radar_ds_prev = pyart.io.read_cfradial(vpt_file_prev)
            nsweeps_prev = radar_ds_prev.nsweeps
            print(f"Loaded previous day data: {nsweeps_prev} sweeps")
            
    except Exception as e:
        print(f"Could not load previous day data: {e}")
    
    # Create figure with more space
    fig, axes = plt.subplots(4, 1, figsize=(12, 18), constrained_layout=False)
    
    # Reserve space for logos - more bottom space for AMOF, less top space for NCAS
    fig.subplots_adjust(left=0.10, right=0.95, top=0.88, bottom=0.12, hspace=0.3)
    
    # Set up time limits for full day
    dtime = cftime.num2pydate(radar_ds.time['data'], radar_ds.time['units'])
    dt_min = dtime[0].replace(hour=0, minute=0, second=0)
    dt_max = dt_min + datetime.timedelta(days=1)
    time_str = dtime[0].strftime("%Y-%m-%d")
    
    # Plot first sweep to establish colorbars and titles
    radar_sweep_ds = radar_ds.extract_sweeps([0])
    gatefilter = pyart.correct.GateFilter(radar_sweep_ds)
    gatefilter.exclude_below('SNR', -20)
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
        gatefilter.exclude_below('SNR', -20)
        display = pyart.graph.RadarDisplay(radar_sweep_ds)
        
        _plot_vpt_fields_overlay(display, axes, gatefilter, vel_limit_lower, vel_limit_upper)
    
    # Final plot setup
    for ax in axes:
        ax.set_xlim(dt_min, dt_max)
        ax.grid(True)
        ax.set_xlabel('Time (UTC)')
    
    # Add logos to figure
    add_logos_to_figure(fig)
    
    # Save figure
    figname = f'ncas-mobile-ka-band-radar-1_cao_cobalt_{datestr}_vpt_l1_{product_version}.png'
    if blflag:
        figname = figname.replace('.png', '_bl.png')
    
    vpt_figpath = os.path.join(figpath, 'vpt')
    os.makedirs(vpt_figpath, exist_ok=True)
    
    plt.savefig(os.path.join(vpt_figpath, figname), dpi=300)
    plt.close()
    
    print(f"VPT plot saved to: {vpt_figpath}")

def _plot_vpt_fields(display, axes, radar_ds, gatefilter, time_str, hmax, 
                    vel_min, vel_max):
    """Plot VPT fields with titles and colorbars."""
    
    # DBZ
    field_name = pyart.graph.common.generate_field_name(radar_ds, "DBZ")
    title = f"{pyart.graph.common.generate_radar_name(radar_ds)} {time_str}\n{field_name}"
    display.plot_vpt("DBZ", ax=axes[0], time_axis_flag=True, title=title, edges=False,
                     gatefilter=gatefilter, vmin=-60, vmax=40, cmap=COLORMAPS['dbz'],
                     colorbar_orient='horizontal', filter_transitions=True)
    axes[0].set_ylim(0, hmax)
    
    # VEL  
    field_name = pyart.graph.common.generate_field_name(radar_ds, "VEL")
    title = f"{pyart.graph.common.generate_radar_name(radar_ds)} {time_str}\n{field_name}"
    display.plot_vpt("VEL", ax=axes[1], time_axis_flag=True, title=title, edges=False,
                     gatefilter=gatefilter, vmin=vel_min, vmax=vel_max, 
                     cmap=COLORMAPS['vel'], colorbar_orient='horizontal')
    axes[1].set_ylim(0, hmax)
    
    # WIDTH
    field_name = pyart.graph.common.generate_field_name(radar_ds, "WIDTH")
    title = f"{pyart.graph.common.generate_radar_name(radar_ds)} {time_str}\n{field_name}"
    display.plot_vpt("WIDTH", ax=axes[2], time_axis_flag=True, title=title, edges=False,
                     gatefilter=gatefilter, norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                     cmap=COLORMAPS['width'], colorbar_orient='horizontal')
    axes[2].set_ylim(0, hmax)
    
    # LDR
    field_name = pyart.graph.common.generate_field_name(radar_ds, "LDR")
    title = f"{pyart.graph.common.generate_radar_name(radar_ds)} {time_str}\n{field_name}"
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

def make_cobalt_rhi_plots_day(datestr, inpath, figpath, blflag=False, no_aircraft=False, skip_all_transition=False):
    """
    Create RHI plots for all files from a given day.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Output directory path  
        blflag: Boundary layer flag
        no_aircraft: If True, skip aircraft markers and related processing
        skip_all_transition: If True, skip sweeps where 100% of rays are antenna_transition=1
    """
    print(f"Creating RHI plots for date: {datestr}")
    
    inpath_date = os.path.join(inpath, datestr)
    print(f"Looking for RHI files in: {inpath_date}")
    
    # Check if directory exists
    if not os.path.exists(inpath_date):
        print(f"ERROR: Directory does not exist: {inpath_date}")
        return
    
    os.chdir(inpath_date)
    
    # Find all RHI files for the date - be more specific about the pattern
    rhi_pattern = f'*{datestr}*rhi*.nc'
    print(f"Searching for files matching pattern: {rhi_pattern}")
    
    all_files = glob.glob('*.nc')
    print(f"All .nc files in directory: {len(all_files)}")
    if all_files:
        print(f"Sample files: {all_files[:5]}")
    
    rhi_files = [os.path.join(inpath_date, f) for f in glob.glob(rhi_pattern)]
    print(f"Found {len(rhi_files)} RHI files matching pattern")
    
    if not rhi_files:
        # Try alternative patterns
        alt_patterns = [
            f'*{datestr}*RHI*.nc',
            f'*rhi*{datestr}*.nc',
            f'*RHI*{datestr}*.nc',
            '*rhi*.nc',
            '*RHI*.nc'
        ]
        
        for pattern in alt_patterns:
            alt_files = glob.glob(pattern)
            if alt_files:
                print(f"Found {len(alt_files)} files with alternative pattern '{pattern}': {alt_files[:3]}")
                rhi_files = [os.path.join(inpath_date, f) for f in alt_files]
                break
    
    if not rhi_files:
        print("No RHI files found with any pattern")
        return
    
    # Set up COBALT command path
    _, _, cobalt_cmd_path = setup_paths(datestr)
    
    # Process each file
    for i, rhi_file in enumerate(rhi_files):
        print(f"Processing RHI file {i+1}/{len(rhi_files)}: {os.path.basename(rhi_file)}")
        try:
            make_cobalt_rhi_plot(rhi_file, figpath, cobalt_cmd_path, blflag=blflag, no_aircraft=no_aircraft, skip_all_transition=skip_all_transition)
        except Exception as e:
            print(f"Error processing {rhi_file}: {e}")
            import traceback
            traceback.print_exc()

def add_logos_to_figure(fig, logo_dir=None):
    """
    Add AMOF and NCAS logos to figure using inset axes.
    
    Args:
        fig: Matplotlib figure object
        logo_dir: Directory containing logo files (optional, not used in this implementation)
    """
    import matplotlib.image as mpimg
    
    # Default logo directory (adjust path as needed)
    if logo_dir is None:
        logo_dir = "/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/logos"
    
    # Logo file paths
    ncas_logo_path = os.path.join(logo_dir, "ncas_logo.png")
    amof_logo_path = os.path.join(logo_dir, "amof_logo.png")
    
    try:
        # Add AMOF logo as an inset axis in the bottom left outside the plot
        if os.path.exists(amof_logo_path):
            amof_logo_img = mpimg.imread(amof_logo_path)
            amof_logo_h = 0.04  # Smaller height to fit in reserved space
            amof_logo_w = amof_logo_h * (amof_logo_img.shape[1] / amof_logo_img.shape[0])
            # Position in the reserved bottom margin
            amof_logo_ax = fig.add_axes([0.02, 0.01, amof_logo_w, amof_logo_h])
            amof_logo_ax.imshow(amof_logo_img)
            amof_logo_ax.axis('off')
            print("Added AMOF logo to figure")
        else:
            print(f"AMOF logo not found at: {amof_logo_path}")

        # Add NCAS logo in the reserved top-right space
        if os.path.exists(ncas_logo_path):
            ncas_logo_img = mpimg.imread(ncas_logo_path)
            ncas_logo_h = 0.05  # Smaller height to fit in reserved space
            ncas_logo_w = ncas_logo_h * (ncas_logo_img.shape[1] / ncas_logo_img.shape[0])
            # Position in the reserved top margin
            ncas_logo_ax = fig.add_axes([1.0 - ncas_logo_w - 0.02, 0.93, ncas_logo_w, ncas_logo_h])
            ncas_logo_ax.imshow(ncas_logo_img)
            ncas_logo_ax.axis('off')
            print("Added NCAS logo to figure")
        else:
            print(f"NCAS logo not found at: {ncas_logo_path}")
            
    except Exception as e:
        print(f"Error adding logos: {e}")

def main():
    """Main function."""
    print(f"COBALT Quicklooks Generator v{VERSION}")
    print("=" * 50)
    
    # Parse command line arguments
    datestr, inpath_arg, outpath_arg, blflag, latest, no_aircraft, skip_all_transition = parse_command_line()
    
    # Set up paths
    inpath, figpath, _ = setup_paths(datestr)
    
    # Override paths if provided via command line
    if inpath_arg:
        inpath = inpath_arg
    if outpath_arg:
        figpath = outpath_arg
    
    print(f"Processing date: {datestr}")
    print(f"Input path: {inpath}")
    print(f"Output path: {figpath}")
    print(f"Boundary layer mode: {blflag}")
    print(f"Aircraft markers: {'disabled' if no_aircraft else 'enabled'}")
    print(f"Skip all-transition sweeps: {'enabled' if skip_all_transition else 'disabled'}")
    
    # Create quicklook plots
    try:
        # RHI plots - ADD no_aircraft parameter
        make_cobalt_rhi_plots_day(datestr, inpath, figpath, blflag=blflag, no_aircraft=no_aircraft, skip_all_transition=skip_all_transition)
        
        # VPT plots  
        make_cobalt_vpt_plot_day(datestr, inpath, figpath, blflag=blflag)
        
        print("Quicklook generation completed successfully!")
        
    except Exception as e:
        print(f"Error in main processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

