#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
make_picasso_quicklooks.py

Generate quicklook plots for PICASSO campaign radar data.

This script creates visualization plots for NCAS Mobile Ka-band Radar data
from the PICASSO (Profiling Ice Clouds by Airborne Scanning from Shetland 
to Arctic Ocean) campaign. It generates RHI, PPI, VPT, and special MAN 
(manual tracking) plots with phase identification.

Usage:
    python make_picasso_quicklooks.py -d YYYYMMDD [-i input_path] [-o output_path] [-b] [--skip-all-transition] [--azimuth-offset DEGREES] [--zoom-km KM] [--zoom-level LEVEL]

Arguments:
    -d, --date:              Date string in YYYYMMDD format
    -i, --inpath:            Input directory containing CF-Radial files (optional)
    -o, --outpath:           Output directory for quicklook images (optional)
    -b:                      Boundary layer flag (limits plots to 4km height)
    --skip-all-transition:   Skip sweeps where 100% of rays are antenna_transition=1
                             (by default, all sweeps are plotted including antenna transitions)
    --azimuth-offset:        Apply azimuth offset in degrees (e.g., --azimuth-offset 5.0)
    --zoom-km:               Distance from radar for map zoom in km (e.g., --zoom-km 5.0)
    --zoom-level:            Basemap tile zoom level 8-17 (higher = more detail, overrides auto-calculation)

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
Last modified: 25-02-2026
Version: 1.0.0
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

# Configure matplotlib to use non-interactive backend (for headless systems)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
from matplotlib.cm import ScalarMappable
import cmocean

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import contextily as ctx
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyproj import Geod

# Import kepler utilities
from kepler_utils import get_valid_sweep_indices

# Configuration
VERSION = 'v1.0.0'
TRACKING_TAG = 'CFARR_0002'
CAMPAIGN = 'picasso'

# Default paths
NCAS_OBS_PROC_PATH = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1'

# Plot configuration
COLORMAPS = {
    'dbz': 'HomeyerRainbow',
    'vel': 'balance', 
    'ldr': 'viridis',
    'width': 'SpectralExtended'
}

# Earth radius in km for distance calculations
EARTH_RADIUS_KM = 6371.0

# Phase colors for MAN scan visualization
PHASE_COLORS = {
    'upward_rhi': '#FF6B6B',         # Red
    'dwell': '#4ECDC4',              # Cyan
    'vertical_pointing': '#00CED1',  # Dark turquoise
    'downward_rhi': '#95E1D3',       # Light green
    'turning': '#FFE66D',            # Yellow
    'other': '#C7CEEA',              # Light purple
    'tracking': '#95A3B3'            # Light gray
}

def parse_command_line():
    """
    Parse command line arguments.
    
    Returns:
        tuple: (datestr, inpath, outpath, blflag, skip_all_transition, custom_bounds, ppi_maps_only, azimuth_offset, zoom_km, zoom_level)
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:i:o:b", 
                                   ["date=", "inpath=", "outpath=", "skip-all-transition",
                                    "lat-min=", "lat-max=", "lon-min=", "lon-max=", "ppi-maps-only",
                                    "azimuth-offset=", "zoom-km=", "zoom-level="])
    except getopt.GetoptError as err:
        print(f"Error: {err}")
        print("Usage: python make_picasso_quicklooks.py -d YYYYMMDD [-i input_path] [-o output_path] [-b] [--skip-all-transition]")
        print("       [--lat-min LAT] [--lat-max LAT] [--lon-min LON] [--lon-max LON] [--ppi-maps-only]")
        sys.exit(2)

    # Default values
    data_date = datetime.datetime.now()
    datestr = data_date.strftime('%Y%m%d')
    inpath = None
    outpath = None
    blflag = False
    skip_all_transition = False
    lat_min = None
    lat_max = None
    lon_min = None
    lon_max = None
    ppi_maps_only = False
    azimuth_offset = 0.0
    zoom_km = None
    zoom_level = None

    for option, argument in opts:
        if option in ("-d", "--date"):
            datestr = argument
        elif option in ("-i", "--inpath"):
            inpath = argument
        elif option in ("-o", "--outpath"):
            outpath = argument
        elif option == "-b":
            blflag = True
        elif option == "--skip-all-transition":
            skip_all_transition = True
        elif option == "--lat-min":
            try:
                lat_min = float(argument)
            except ValueError:
                print(f"Error: lat-min must be a number, got '{argument}'")
                sys.exit(2)
        elif option == "--lat-max":
            try:
                lat_max = float(argument)
            except ValueError:
                print(f"Error: lat-max must be a number, got '{argument}'")
                sys.exit(2)
        elif option == "--lon-min":
            try:
                lon_min = float(argument)
            except ValueError:
                print(f"Error: lon-min must be a number, got '{argument}'")
                sys.exit(2)
        elif option == "--lon-max":
            try:
                lon_max = float(argument)
            except ValueError:
                print(f"Error: lon-max must be a number, got '{argument}'")
                sys.exit(2)
        elif option == "--ppi-maps-only":
            ppi_maps_only = True
        elif option == "--azimuth-offset":
            try:
                azimuth_offset = float(argument)
                print(f"Applying azimuth offset: {azimuth_offset}°")
            except ValueError:
                print(f"Error: azimuth-offset must be a number, got '{argument}'")
                sys.exit(2)
        elif option == "--zoom-km":
            try:
                zoom_km = float(argument)
                if zoom_km <= 0:
                    print(f"Error: zoom-km must be positive, got {zoom_km}")
                    sys.exit(2)
            except ValueError:
                print(f"Error: Invalid zoom-km '{argument}'. Must be a positive number.")
                sys.exit(2)
        elif option == "--zoom-level":
            try:
                zoom_level = int(argument)
                if not (8 <= zoom_level <= 17):
                    print(f"Error: zoom-level must be between 8 and 17 (OpenTopo limit), got {zoom_level}")
                    sys.exit(2)
            except ValueError:
                print(f"Error: Invalid zoom-level '{argument}'. Must be an integer.")
                sys.exit(2)
        else:
            print(f"Unhandled option: {option}")
            sys.exit(2)

    # Validate bounds if any are specified
    custom_bounds = None
    if any(b is not None for b in [lat_min, lat_max, lon_min, lon_max]):
        if not all(b is not None for b in [lat_min, lat_max, lon_min, lon_max]):
            print("Error: If specifying custom bounds, all four values must be provided: --lat-min, --lat-max, --lon-min, --lon-max")
            sys.exit(2)
        if lat_min >= lat_max:
            print(f"Error: lat-min ({lat_min}) must be less than lat-max ({lat_max})")
            sys.exit(2)
        if lon_min >= lon_max:
            print(f"Error: lon-min ({lon_min}) must be less than lon-max ({lon_max})")
            sys.exit(2)
        custom_bounds = (lat_min, lat_max, lon_min, lon_max)
        print(f"Using custom map bounds: lat=[{lat_min}, {lat_max}], lon=[{lon_min}, {lon_max}]")

    return datestr, inpath, outpath, blflag, skip_all_transition, custom_bounds, ppi_maps_only, azimuth_offset, zoom_km, zoom_level

def setup_paths(datestr, inpath=None, outpath=None):
    """
    Set up input and figure output paths.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Optional input path override
        outpath: Optional output path override
        
    Returns:
        tuple: (inpath, figpath)
    """
    # Input path for processed radar data
    if inpath is None:
        inpath = os.path.join(NCAS_OBS_PROC_PATH, CAMPAIGN, 'L1_v1.0.0')
    
    # Output path for figures
    if outpath is None:
        figpath = os.path.join(inpath, 'quicklooks')
    else:
        figpath = outpath

    return inpath, figpath

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
                # Try to construct from time variable
                time_units = ds.variables['time'].units
                time_ref = time_units.split('since ')[-1]
                return time_ref
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

def setup_plot_limits(blflag, scan_type='rhi'):
    """
    Set up plot limits based on scan type and boundary layer flag.
    
    Args:
        blflag: Boundary layer flag (limits height to 4km)
        scan_type: Type of scan ('rhi', 'ppi', 'vpt', 'man')
        
    Returns:
        tuple: (hmax, xmin, xmax) - height max and horizontal limits
    """
    if blflag:
        return 4, 0, 40  # Boundary layer mode: lower height
    else:
        return 14, 0, 40  # Full troposphere mode

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

def calculate_zoom_level(lat_min, lat_max, lon_min, lon_max, map_width_pixels=1000, override_zoom=None):
    """
    Calculate appropriate zoom level based on bounding box.
    
    Args:
        lat_min, lat_max, lon_min, lon_max: Bounding box coordinates
        map_width_pixels: Expected map width in pixels
        override_zoom: Optional zoom level override (8-17)
        
    Returns:
        int: Appropriate zoom level (1-20)
    """
    import math
    
    # If override provided, use it
    if override_zoom is not None:
        print(f"Using manual zoom level: {override_zoom}")
        return override_zoom
    
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
    
    # Cap zoom level to match OpenTopo tile provider limit
    zoom = min(zoom, 17)
    zoom = max(zoom, 8)

    print(f"Bounding box: {lat_span:.4f}° lat x {lon_span:.4f}° lon")
    print(f"Max range: {max_span:.4f}°, calculated zoom level: {zoom}")

    return zoom

def split_radar_at_time_gaps(radar, max_gap_seconds=60):
    """
    Split a radar object into multiple segments at large time gaps.
    
    This prevents pixels from stretching across periods when radar was not collecting data.
    
    Args:
        radar: PyART Radar object
        max_gap_seconds: Maximum allowed time gap between consecutive rays (default 60s)
        
    Returns:
        list: List of PyART Radar objects, one per continuous segment
    """
    # Get time data
    times = cftime.num2pydate(radar.time['data'], radar.time['units'])
    
    # Calculate time differences between consecutive rays in seconds
    time_diffs = np.diff([t.timestamp() for t in times])
    
    # Find gaps larger than threshold
    gap_indices = np.where(time_diffs > max_gap_seconds)[0]
    
    # If no large gaps, return original radar
    if len(gap_indices) == 0:
        return [radar]
    
    # Split radar at gap locations
    segments = []
    start_ray = 0
    
    for gap_idx in gap_indices:
        # Create segment from start_ray to gap_idx (inclusive)
        end_ray = gap_idx
        if end_ray >= start_ray:
            ray_indices = np.arange(start_ray, end_ray + 1)
            segment = radar.extract_sweeps([0])  # Get structure
            
            # Extract the rays for this segment
            segment.time['data'] = radar.time['data'][ray_indices]
            segment.range['data'] = radar.range['data']
            
            # Extract field data for these rays
            for field_name in radar.fields.keys():
                field_data = radar.fields[field_name]['data'][ray_indices, :]
                segment.fields[field_name]['data'] = field_data
            
            # Update metadata
            segment.azimuth['data'] = radar.azimuth['data'][ray_indices]
            segment.elevation['data'] = radar.elevation['data'][ray_indices]
            segment.nsweeps = 1
            segment.sweep_number['data'] = np.array([0])
            segment.sweep_mode['data'] = np.array(['vpt'])
            segment.fixed_angle['data'] = np.array([90.0])
            segment.sweep_start_ray_index['data'] = np.array([0])
            segment.sweep_end_ray_index['data'] = np.array([len(ray_indices) - 1])
            segment.nrays = len(ray_indices)
            segment.ngates = len(segment.range['data'])
            
            segments.append(segment)
        
        start_ray = gap_idx + 1
    
    # Add final segment
    if start_ray < len(times):
        ray_indices = np.arange(start_ray, len(times))
        segment = radar.extract_sweeps([0])
        
        segment.time['data'] = radar.time['data'][ray_indices]
        segment.range['data'] = radar.range['data']
        
        for field_name in radar.fields.keys():
            field_data = radar.fields[field_name]['data'][ray_indices, :]
            segment.fields[field_name]['data'] = field_data
        
        segment.azimuth['data'] = radar.azimuth['data'][ray_indices]
        segment.elevation['data'] = radar.elevation['data'][ray_indices]
        segment.nsweeps = 1
        segment.sweep_number['data'] = np.array([0])
        segment.sweep_mode['data'] = np.array(['vpt'])
        segment.fixed_angle['data'] = np.array([90.0])
        segment.sweep_start_ray_index['data'] = np.array([0])
        segment.sweep_end_ray_index['data'] = np.array([len(ray_indices) - 1])
        segment.nrays = len(ray_indices)
        segment.ngates = len(segment.range['data'])
        
        segments.append(segment)
    
    print(f"  Split radar data into {len(segments)} segments at time gaps > {max_gap_seconds}s")
    return segments

def make_picasso_rhi_plot(ncfile, figpath, blflag=False, skip_all_transition=False):
    """
    Create RHI (Range Height Indicator) quicklook plot for PICASSO data.
    
    Args:
        ncfile: Path to CF-Radial NetCDF file
        figpath: Output directory for figure
        blflag: Boundary layer flag for plot limits
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
        product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
        instrument_name = ds.source if 'source' in ds.ncattrs() else 'ncas-mobile-ka-band-radar-1'
    
    # Set up plot limits
    hmax, xmin, xmax = setup_plot_limits(blflag, 'rhi')
    
    # Read radar data
    try:
        radar_ds = pyart.io.read_cfradial(ncfile)
    except Exception as e:
        print(f"Error reading radar file: {e}")
        return
    
    # Set up gate filter
    gatefilter = pyart.correct.GateFilter(radar_ds)
    gatefilter.exclude_below('SNR', -20)
    
    # Exclude rays flagged as antenna_transition (repositioning between sweeps)
    if hasattr(radar_ds, 'antenna_transition') and radar_ds.antenna_transition is not None:
        if 'data' in radar_ds.antenna_transition:
            transition_flags = radar_ds.antenna_transition['data']
            transition_indices = np.where(transition_flags == 1)[0]
            if len(transition_indices) > 0:
                # Exclude all gates in rays flagged as antenna_transition
                # Properly access the gate_excluded attribute (2D boolean array)
                for ray_idx in transition_indices:
                    if hasattr(gatefilter, 'gate_excluded'):
                        gatefilter.gate_excluded[ray_idx, :] = True
                    elif hasattr(gatefilter, '_gate_excluded'):
                        gatefilter._gate_excluded[ray_idx, :] = True
    
    # Create display object
    display = pyart.graph.RadarDisplay(radar_ds)
    
    # Set up output directory
    rhi_figpath = os.path.join(figpath, 'rhi', netcdf_time.strftime('%Y%m%d'))
    os.makedirs(rhi_figpath, exist_ok=True)
    
    # Get velocity limits from data
    vel_field = radar_ds.fields['VEL']
    vel_limit_lower = vel_field.get('field_limit_lower', -10)
    vel_limit_upper = vel_field.get('field_limit_upper', 10)
    
    # Create plots for each sweep
    nsweeps = radar_ds.nsweeps
    dtime0 = cftime.num2pydate(radar_ds.time['data'][0], radar_ds.time['units'])
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S")
    
    # Get valid sweep indices (optionally filtering out all-transition sweeps)
    valid_sweep_indices = get_valid_sweep_indices(radar_ds, skip_all_transition=skip_all_transition)
    
    if skip_all_transition and len(valid_sweep_indices) < nsweeps:
        print(f"Skipping {nsweeps - len(valid_sweep_indices)} sweep(s) that are 100% antenna transitions")
    
    for sweep_idx in valid_sweep_indices:
        print(f"Processing sweep {sweep_idx+1}/{nsweeps}")
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Get sweep start and end times
        start_ray_idx = radar_ds.sweep_start_ray_index['data'][sweep_idx]
        end_ray_idx = radar_ds.sweep_end_ray_index['data'][sweep_idx]
        
        # Calculate actual mean azimuth for this sweep (handles wrapping)
        azimuth_data = radar_ds.azimuth['data'][start_ray_idx:end_ray_idx+1]
        rhi_az = np.degrees(np.arctan2(
            np.mean(np.sin(np.radians(azimuth_data))),
            np.mean(np.cos(np.radians(azimuth_data)))
        )) % 360
        
        start_time = cftime.num2pydate(radar_ds.time['data'][start_ray_idx], radar_ds.time['units'])
        end_time = cftime.num2pydate(radar_ds.time['data'][end_ray_idx], radar_ds.time['units'])
        
        time_range_str = f"{start_time.strftime('%H:%M:%S')} - {end_time.strftime('%H:%M:%S')} UTC"
        
        # Plot each field
        _plot_rhi_fields(display, axes, sweep_idx, gatefilter, vel_limit_lower, 
                       vel_limit_upper, hmax, xmin, xmax)
        
        # Add title
        fig.suptitle(f'PICASSO RHI Scan - Azimuth: {rhi_az:.1f}°\n'
                    f'{instrument_name}\n'
                    f'{start_time.strftime("%Y-%m-%d")} | {time_range_str}',
                    fontsize=12, fontweight='bold')
        
        # Save plot
        if nsweeps > 1:
            sweep_start_index = radar_ds.get_start(sweep_idx)
            dtime_sweep = cftime.num2pydate(radar_ds.time['data'][sweep_start_index], 
                                          radar_ds.time['units'])
            dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S")
            figname = (f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime_sweep_str}'
                      f'_rhi_az{rhi_az:.2f}_l1_{product_version}.png')
        else:
            figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime0_str}_rhi_az{rhi_az:.2f}_l1_{product_version}.png'
        
        plt.tight_layout()
        plt.savefig(os.path.join(rhi_figpath, figname), dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"RHI plot saved to: {rhi_figpath}")

def _plot_rhi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max, 
                    hmax, xmin, xmax):
    """Helper function to plot RHI fields on axes."""
    
    radar = display._radar
    
    # Reflectivity (DBZ)
    display.plot_rhi("DBZ", ax=axes[0,0], sweep=sweep_idx, vmin=-60, vmax=40,
                     cmap=COLORMAPS['dbz'], colorbar_flag=False,
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[0,0], hmax, xmin, xmax)
    dbz_title = radar.fields['DBZ'].get('long_name', 'Reflectivity')
    axes[0,0].set_title(dbz_title)
    if len(axes[0,0].collections) > 0:
        plt.colorbar(axes[0,0].collections[-1], ax=axes[0,0], orientation='horizontal', 
                     label='Reflectivity (dBZ)', pad=0.25, aspect=30)
    
    # Velocity (VEL)  
    display.plot_rhi("VEL", ax=axes[1,0], sweep=sweep_idx, vmin=vel_min, vmax=vel_max,
                     cmap=COLORMAPS['vel'], colorbar_flag=False,
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[1,0], hmax, xmin, xmax)
    vel_title = radar.fields['VEL'].get('long_name', 'Doppler Velocity')
    axes[1,0].set_title(vel_title)
    if len(axes[1,0].collections) > 0:
        plt.colorbar(axes[1,0].collections[-1], ax=axes[1,0], orientation='horizontal',
                     label='Velocity (m/s)', pad=0.25, aspect=30)
    
    # Spectrum width (WIDTH)
    display.plot_rhi("WIDTH", ax=axes[1,1], sweep=sweep_idx,
                     norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                     cmap=COLORMAPS['width'], colorbar_flag=False,
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[1,1], hmax, xmin, xmax)
    width_title = radar.fields['WIDTH'].get('long_name', 'Spectrum Width')
    axes[1,1].set_title(width_title)
    if len(axes[1,1].collections) > 0:
        plt.colorbar(axes[1,1].collections[-1], ax=axes[1,1], orientation='horizontal',
                     label='Spectrum Width (m/s)', pad=0.25, aspect=30)
    
    # Linear depolarization ratio (LDR)
    display.plot_rhi("LDR", ax=axes[0,1], sweep=sweep_idx, vmin=-35, vmax=5,
                     cmap=COLORMAPS['ldr'], colorbar_flag=False,
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[0,1], hmax, xmin, xmax)
    ldr_title = radar.fields['LDR'].get('long_name', 'Linear Depolarization Ratio')
    axes[0,1].set_title(ldr_title)
    if len(axes[0,1].collections) > 0:
        plt.colorbar(axes[0,1].collections[-1], ax=axes[0,1], orientation='horizontal',
                     label='LDR (dB)', pad=0.25, aspect=30)

def _setup_rhi_axes(ax, hmax, xmin, xmax):
    """Helper function to set up RHI plot axes."""
    ax.set_ylim(0, hmax)
    ax.set_xlim(xmin, xmax)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

def make_picasso_ppi_plot(ncfile, figpath, blflag=False, skip_all_transition=False, azimuth_offset=0.0):
    """
    Create PPI (Plan Position Indicator) quicklook plot for PICASSO data.
    
    Args:
        ncfile: Path to CF-Radial NetCDF file
        figpath: Output directory for figure
        blflag: Boundary layer flag for plot limits
        skip_all_transition: If True, skip sweeps where 100% of rays are antenna_transition=1
        azimuth_offset: Additional azimuth offset to apply in degrees (default: 0.0)
    
    Note:
        Azimuth data in the NetCDF file is already corrected during processing.
        Additional azimuth offset can be applied via azimuth_offset parameter.
    """
    print(f"Creating PPI plot for: {ncfile}")
    
    # Get scan information
    time_coverage_start = get_time_coverage_start(ncfile)
    if not time_coverage_start:
        print(f"Could not get time coverage from {ncfile}")
        return
        
    netcdf_time = datetime.datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%SZ")
    
    with nc4.Dataset(ncfile) as ds:
        scan_elevation = ds['fixed_angle'][0]
        product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
        instrument_name = ds.source if 'source' in ds.ncattrs() else 'ncas-mobile-ka-band-radar-1'
    
    # Set maximum range for PPI plots (km)
    max_range = 40.0
    
    # Read radar data
    try:
        radar_ds = pyart.io.read_cfradial(ncfile)
    except Exception as e:
        print(f"Error reading radar file: {e}")
        return
    
    # Apply azimuth offset if specified
    if azimuth_offset != 0.0:
        az_before = radar_ds.azimuth['data'].copy()
        radar_ds.azimuth['data'] = (radar_ds.azimuth['data'] + azimuth_offset) % 360
        print(f"  Applying azimuth offset: {azimuth_offset}°")
        print(f"    Original azimuth range: [{az_before.min():.1f}, {az_before.max():.1f}]°")
        print(f"    Modified azimuth range: [{radar_ds.azimuth['data'].min():.1f}, {radar_ds.azimuth['data'].max():.1f}]°")
        print(f"    Example: {az_before[0]:.2f}° -> {radar_ds.azimuth['data'][0]:.2f}°")
        print(f"    Visual effect: Data will rotate {azimuth_offset:.1f}° {'clockwise' if azimuth_offset > 0 else 'counter-clockwise'}")
        print(f"                   (e.g., North data moves to {azimuth_offset:.1f}° azimuth)")
        
        # Clear any cached gate locations that PyART might have computed
        # PyART caches gate_x, gate_y, gate_z, gate_altitude, gate_latitude, gate_longitude
        for attr in ['gate_x', 'gate_y', 'gate_z', 'gate_altitude', 'gate_latitude', 'gate_longitude']:
            if hasattr(radar_ds, attr):
                delattr(radar_ds, attr)
                print(f"    Cleared cached attribute: {attr}")
    
    # Set up gate filter
    gatefilter = pyart.correct.GateFilter(radar_ds)
    gatefilter.exclude_below('SNR', -20)
    
    # Exclude rays flagged as antenna_transition (repositioning between sweeps)
    if hasattr(radar_ds, 'antenna_transition') and radar_ds.antenna_transition is not None:
        if 'data' in radar_ds.antenna_transition:
            transition_flags = radar_ds.antenna_transition['data']
            transition_indices = np.where(transition_flags == 1)[0]
            if len(transition_indices) > 0:
                # Exclude all gates in rays flagged as antenna_transition
                # Properly access the gate_excluded attribute (2D boolean array)
                for ray_idx in transition_indices:
                    if hasattr(gatefilter, 'gate_excluded'):
                        gatefilter.gate_excluded[ray_idx, :] = True
                    elif hasattr(gatefilter, '_gate_excluded'):
                        gatefilter._gate_excluded[ray_idx, :] = True
    
    # Create display object
    display = pyart.graph.RadarDisplay(radar_ds)
    
    # Set up output directory
    ppi_figpath = os.path.join(figpath, 'ppi', netcdf_time.strftime('%Y%m%d'))
    os.makedirs(ppi_figpath, exist_ok=True)
    
    # Get velocity limits from data
    vel_field = radar_ds.fields['VEL']
    vel_limit_lower = vel_field.get('field_limit_lower', -10)
    vel_limit_upper = vel_field.get('field_limit_upper', 10)
    
    # Create plots for each sweep
    nsweeps = radar_ds.nsweeps
    dtime0 = cftime.num2pydate(radar_ds.time['data'][0], radar_ds.time['units'])
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S")
    
    # Get valid sweep indices (optionally filtering out all-transition sweeps)
    valid_sweep_indices = get_valid_sweep_indices(radar_ds, skip_all_transition=skip_all_transition)
    
    if skip_all_transition and len(valid_sweep_indices) < nsweeps:
        print(f"Skipping {nsweeps - len(valid_sweep_indices)} sweep(s) that are 100% antenna transitions")
    
    for sweep_idx in valid_sweep_indices:
        print(f"Processing sweep {sweep_idx+1}/{nsweeps}")
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 14))
        
        # Get sweep start and end times
        start_ray_idx = radar_ds.sweep_start_ray_index['data'][sweep_idx]
        end_ray_idx = radar_ds.sweep_end_ray_index['data'][sweep_idx]
        
        # Calculate actual mean elevation for this sweep
        elevation_data = radar_ds.elevation['data'][start_ray_idx:end_ray_idx+1]
        ppi_el = np.mean(elevation_data)
        
        start_time = cftime.num2pydate(radar_ds.time['data'][start_ray_idx], radar_ds.time['units'])
        end_time = cftime.num2pydate(radar_ds.time['data'][end_ray_idx], radar_ds.time['units'])
        
        time_range_str = f"{start_time.strftime('%H:%M:%S')} - {end_time.strftime('%H:%M:%S')} UTC"
        
        # Plot each field
        _plot_ppi_fields(display, axes, sweep_idx, gatefilter, vel_limit_lower, 
                        vel_limit_upper, max_range)
        
        # Add title
        fig.suptitle(f'PICASSO PPI Scan - Elevation: {ppi_el:.1f}°\n'
                    f'{instrument_name}\n'
                    f'{start_time.strftime("%Y-%m-%d")} | {time_range_str}',
                    fontsize=12, fontweight='bold')
        
        # Save plot
        if nsweeps > 1:
            sweep_start_index = radar_ds.get_start(sweep_idx)
            dtime_sweep = cftime.num2pydate(radar_ds.time['data'][sweep_start_index], 
                                          radar_ds.time['units'])
            dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S")
            figname = (f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime_sweep_str}'
                      f'_ppi_el{ppi_el:.2f}_l1_{product_version}.png')
        else:
            figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime0_str}_ppi_el{ppi_el:.2f}_l1_{product_version}.png'
        
        plt.tight_layout()
        plt.savefig(os.path.join(ppi_figpath, figname), dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"PPI plot saved to: {ppi_figpath}")

def parse_geospatial_bounds(bounds_str):
    """
    Parse geospatial_bounds string from NetCDF file.
    
    Expected format: "Bounding box: 51.00N -1.43E, 51.14N -0.78E"
    where the first pair is (lat_min, lon_min) and second is (lat_max, lon_max).
    
    Args:
        bounds_str: String containing bounding box information
        
    Returns:
        tuple: (lat_min, lat_max, lon_min, lon_max) or None if parsing fails
    """
    try:
        # Remove "Bounding box: " prefix if present
        if "Bounding box:" in bounds_str:
            bounds_str = bounds_str.split("Bounding box:")[-1].strip()
        
        # Split into two coordinate pairs
        parts = bounds_str.split(',')
        if len(parts) != 2:
            return None
        
        # Parse first pair (lat_min, lon_min)
        first_pair = parts[0].strip().split()
        if len(first_pair) != 2:
            return None
        lat_min = float(first_pair[0].replace('N', '').replace('S', ''))
        lon_min = float(first_pair[1].replace('E', '').replace('W', ''))
        
        # Parse second pair (lat_max, lon_max)
        second_pair = parts[1].strip().split()
        if len(second_pair) != 2:
            return None
        lat_max = float(second_pair[0].replace('N', '').replace('S', ''))
        lon_max = float(second_pair[1].replace('E', '').replace('W', ''))
        
        return (lat_min, lat_max, lon_min, lon_max)
        
    except Exception as e:
        print(f"Could not parse geospatial_bounds: {e}")
        return None

def make_picasso_ppi_map_plot(ncfile, figpath, blflag=False, skip_all_transition=False, 
                               basemap_type="opentopo", lat_bounds=None, lon_bounds=None, 
                               zoom_km=None, azimuth_offset=0.0, zoom_level=None):
    """
    Create PPI map plots with geographic basemap for PICASSO data.
    
    Args:
        ncfile: Path to CF-Radial NetCDF file
        figpath: Output directory for figure
        blflag: Boundary layer flag (not used for map plots)
        skip_all_transition: If True, skip sweeps where 100% of rays are antenna_transition=1
        basemap_type: Type of basemap ("opentopo", "satellite", "cartodb")
        lat_bounds: Optional tuple (lat_min, lat_max) for map extent
        lon_bounds: Optional tuple (lon_min, lon_max) for map extent
        zoom_km: Optional zoom distance in km from radar (overrides lat_bounds/lon_bounds)
        azimuth_offset: Additional azimuth offset to apply in degrees (default: 0.0)
        zoom_level: Optional basemap tile zoom level (8-17, overrides auto-calculation)
    
    Note:
        Azimuth data in the NetCDF file is already corrected during processing.
        Additional azimuth offset can be applied via azimuth_offset parameter.
    """
    print(f"Creating PPI map plot for: {ncfile}")
    
    # Get scan information
    time_coverage_start = get_time_coverage_start(ncfile)
    if not time_coverage_start:
        print(f"Could not get time coverage from {ncfile}")
        return
        
    netcdf_time = datetime.datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%SZ")
    
    # Read metadata from NetCDF file
    geospatial_bounds_str = None
    with nc4.Dataset(ncfile) as ds:
        scan_elevation = ds['fixed_angle'][0]
        product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
        instrument_name = ds.source if 'source' in ds.ncattrs() else 'ncas-mobile-ka-band-radar-1'
        if 'geospatial_bounds' in ds.ncattrs():
            geospatial_bounds_str = ds.getncattr('geospatial_bounds')
    
    # Read radar data
    try:
        radar = pyart.io.read_cfradial(ncfile)
    except Exception as e:
        print(f"Error reading radar file: {e}")
        return
    
    # Apply azimuth offset if specified
    if azimuth_offset != 0.0:
        az_before = radar.azimuth['data'].copy()
        radar.azimuth['data'] = (radar.azimuth['data'] + azimuth_offset) % 360
        print(f"  Applying azimuth offset: {azimuth_offset}°")
        print(f"    Original azimuth range: [{az_before.min():.1f}, {az_before.max():.1f}]°")
        print(f"    Modified azimuth range: [{radar.azimuth['data'].min():.1f}, {radar.azimuth['data'].max():.1f}]°")
        print(f"    Example: {az_before[0]:.2f}° -> {radar.azimuth['data'][0]:.2f}°")
        print(f"    Visual effect: Data will rotate {azimuth_offset:.1f}° {'clockwise' if azimuth_offset > 0 else 'counter-clockwise'}")
        print(f"                   (e.g., North data moves to {azimuth_offset:.1f}° azimuth)")
        
        # Clear any cached gate locations that PyART might have computed
        # PyART caches gate_x, gate_y, gate_z, gate_altitude, gate_latitude, gate_longitude
        for attr in ['gate_x', 'gate_y', 'gate_z', 'gate_altitude', 'gate_latitude', 'gate_longitude']:
            if hasattr(radar, attr):
                delattr(radar, attr)
                print(f"    Cleared cached attribute: {attr}")
    
    # Get radar location
    radar_lat = radar.latitude['data'][0]
    radar_lon = radar.longitude['data'][0]
    
    # Set up projection
    projection = ccrs.AzimuthalEquidistant(central_longitude=radar_lon, central_latitude=radar_lat)
    
    # Set up gate filter
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_below('SNR', -20)
    
    # Exclude rays flagged as antenna_transition
    if hasattr(radar, 'antenna_transition') and radar.antenna_transition is not None:
        if 'data' in radar.antenna_transition:
            transition_flags = radar.antenna_transition['data']
            transition_indices = np.where(transition_flags == 1)[0]
            if len(transition_indices) > 0:
                for ray_idx in transition_indices:
                    if hasattr(gatefilter, 'gate_excluded'):
                        gatefilter.gate_excluded[ray_idx, :] = True
                    elif hasattr(gatefilter, '_gate_excluded'):
                        gatefilter._gate_excluded[ray_idx, :] = True
    
    # Create display object
    display = pyart.graph.RadarMapDisplay(radar)
    
    # Set up output directory
    if zoom_km is not None:
        ppi_map_figpath = os.path.join(figpath, 'ppi_map_zoom', netcdf_time.strftime('%Y%m%d'))
    else:
        ppi_map_figpath = os.path.join(figpath, 'ppi_map', netcdf_time.strftime('%Y%m%d'))
    os.makedirs(ppi_map_figpath, exist_ok=True)
    
    # Calculate average ray duration
    avg_ray_duration = get_average_ray_duration(radar)
    if avg_ray_duration is not None:
        ray_duration_str = f"{avg_ray_duration:.1f}ms"
    else:
        ray_duration_str = "unknown"
    
    # Get velocity limits from data
    vel_field = radar.fields['VEL']
    vel_limit_lower = vel_field.get('field_limit_lower', -10)
    vel_limit_upper = vel_field.get('field_limit_upper', 10)
    
    # Define colormaps and limits
    colormaps = {
        "DBZ": "HomeyerRainbow",
        "VEL": "balance",
        "WIDTH": "viridis",
        "LDR": "viridis"
    }
    
    default_vmin_vmax = {
        "DBZ": (-20, 50),
        "VEL": (vel_limit_lower, vel_limit_upper),
        "WIDTH": (1e-1 * np.sqrt(1e-1), np.sqrt(1e1)),
        "LDR": (-35, 5)
    }
    
    # Set map bounds - priority order: zoom_km, explicit bounds, NetCDF bounds, calculated
    if zoom_km is not None:
        # Calculate bounds from zoom distance in km
        # 1 km in latitude ≈ 1/111.32 degrees
        lat_offset = zoom_km / 111.32
        # 1 km in longitude depends on latitude: 1/(111.32 * cos(lat_radians))
        lon_offset = zoom_km / (111.32 * np.cos(np.radians(radar_lat)))
        
        lat_min = radar_lat - lat_offset
        lat_max = radar_lat + lat_offset
        lon_min = radar_lon - lon_offset
        lon_max = radar_lon + lon_offset
        print(f"Using zoom bounds: +/- {zoom_km} km from radar")
    elif lat_bounds is not None and lon_bounds is not None:
        # Use explicitly provided bounds
        lat_min, lat_max = lat_bounds
        lon_min, lon_max = lon_bounds
        print(f"Using user-provided bounds")
    elif geospatial_bounds_str is not None:
        # Try to parse geospatial_bounds from NetCDF file
        bounds = parse_geospatial_bounds(geospatial_bounds_str)
        if bounds is not None:
            lat_min, lat_max, lon_min, lon_max = bounds
            print(f"Using geospatial_bounds from NetCDF file: {geospatial_bounds_str}")
        else:
            # Fall back to calculated bounds
            print(f"Could not parse geospatial_bounds, using calculated bounds")
            lat_min, lat_max = radar_lat - 0.4, radar_lat + 0.4
            lon_factor = 1.0 / np.cos(np.radians(radar_lat))
            lon_offset = 0.4 * lon_factor
            lon_min, lon_max = radar_lon - lon_offset, radar_lon + lon_offset
    else:
        # Fall back to calculated bounds based on radar range (40km)
        print(f"No geospatial_bounds in NetCDF file, using calculated bounds")
        # Roughly 40km in latitude is ~0.36 degrees
        lat_min, lat_max = radar_lat - 0.4, radar_lat + 0.4
        # Adjust for latitude when converting longitude range
        # At 60°N, 1° longitude ≈ 55.8 km, so 40km ≈ 0.72°
        lon_factor = 1.0 / np.cos(np.radians(radar_lat))
        lon_offset = 0.4 * lon_factor
        lon_min, lon_max = radar_lon - lon_offset, radar_lon + lon_offset
    
    print(f"Using map bounds: lat [{lat_min:.3f}, {lat_max:.3f}], lon [{lon_min:.3f}, {lon_max:.3f}]")
    
    # Calculate zoom level (or use provided override)
    zoom_level = calculate_zoom_level(lat_min, lat_max, lon_min, lon_max, override_zoom=zoom_level)
    
    # Maximum radar range in km
    max_range_km = 40.0
    
    # Get valid sweep indices
    valid_sweep_indices = get_valid_sweep_indices(radar, skip_all_transition=skip_all_transition)
    
    if skip_all_transition and len(valid_sweep_indices) < radar.nsweeps:
        print(f"Skipping {radar.nsweeps - len(valid_sweep_indices)} sweep(s) that are 100% antenna transitions")
    
    # Process each sweep
    for sweep_idx in valid_sweep_indices:
        print(f"Processing map for sweep {sweep_idx+1}/{radar.nsweeps}")
        
        # Get sweep information
        ppi_el = radar.get_elevation(sweep_idx)[0]
        dtime_sweep = cftime.num2pydate(radar.time['data'][radar.get_start(sweep_idx)], 
                                        radar.time['units'])
        dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S")
        
        # Create figure with 2x2 grid
        variables = ["DBZ", "VEL", "WIDTH", "LDR"]
        positions = [(0, 0), (1, 0), (0, 1), (1, 1)]
        fig, ax = plt.subplots(2, 2, figsize=(30, 30), subplot_kw=dict(projection=projection))
        fig.subplots_adjust(wspace=0.1, hspace=0.1)
        
        # Plot each variable
        for var, pos in zip(variables, positions):
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
                var_name = var
            
            # Plot the PPI map
            display.plot_ppi_map(
                var, sweep_idx,
                vmin=default_vmin_vmax[var][0], vmax=default_vmin_vmax[var][1],
                min_lat=lat_min - 5, max_lat=lat_max + 5,  # Keep PyART's internal limits wide
                min_lon=lon_min - 5, max_lon=lon_max + 5,
                lon_lines=None,
                lat_lines=None,
                projection=projection,
                fig=fig,
                ax=ax[pos],
                lat_0=radar_lat,
                lon_0=radar_lon,
                cmap=cmap_orig,
                norm=norm,
                alpha=0.5,
                colorbar_flag=False,
                gatefilter=gatefilter,
                edges=False,
                embellish=False,
                resolution="10m",
                title=f"{var_name}\nElevation: {ppi_el:.2f}° | Date: {dtime_sweep.strftime('%Y-%m-%d')} | Time: {dtime_sweep.strftime('%H:%M:%S')} | Ray duration: {ray_duration_str}",
            )
            
            # Create colorbar
            colorbar_mappable = ScalarMappable(norm=norm, cmap=cmap_orig)
            colorbar_mappable.set_array([])
            cbar = fig.colorbar(
                colorbar_mappable,
                ax=ax[pos],
                orientation='horizontal',
                shrink=0.6,
                pad=0.12
            )
            cbar.set_alpha(1.0)
            
            # Add label to colorbar
            if var in radar.fields:
                long_name = radar.fields[var].get('long_name', var)
                units = radar.fields[var].get('units', '')
                cbar.set_label(f"{long_name} ({units})", fontsize=12)
            else:
                cbar.set_label(f"{var}", fontsize=12)
            
            # Plot range rings
            display.plot_range_rings([10, 20, 30, 40], ax=ax[pos], lw=0.5, col='0.5')
        
        # Set map extents and add basemap for all panels
        for pos in positions:
            ax[pos].set_extent([lon_min, lon_max, lat_min, lat_max], ccrs.PlateCarree())
            
            # Add basemap
            if basemap_type == "satellite":
                try:
                    ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.Esri.WorldImagery, 
                                  crs=projection)
                    print(f"Added satellite basemap with zoom level {zoom_level}")
                except Exception as e:
                    print(f"Failed to load satellite imagery, falling back to OpenTopoMap: {e}")
                    try:
                        ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.OpenTopoMap, 
                                      crs=projection)
                    except Exception as e2:
                        print(f"OpenTopoMap also failed, using lower zoom: {e2}")
                        ctx.add_basemap(ax[pos], zoom=max(1, zoom_level-2), 
                                      source=ctx.providers.OpenTopoMap, crs=projection)
            elif basemap_type == "cartodb":
                try:
                    ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.CartoDB.Positron, 
                                  crs=projection)
                    print(f"Added CartoDB basemap with zoom level {zoom_level}")
                except Exception as e:
                    print(f"Failed to load CartoDB, trying lower zoom: {e}")
                    ctx.add_basemap(ax[pos], zoom=max(1, zoom_level-2), 
                                  source=ctx.providers.CartoDB.Positron, crs=projection)
            else:  # Default to opentopo
                try:
                    ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.OpenTopoMap, 
                                  crs=projection)
                    print(f"Added OpenTopoMap basemap with zoom level {zoom_level}")
                except Exception as e:
                    print(f"Failed to load OpenTopoMap, falling back to CartoDB: {e}")
                    try:
                        ctx.add_basemap(ax[pos], zoom=zoom_level, source=ctx.providers.CartoDB.Positron, 
                                      crs=projection)
                    except Exception as e2:
                        print(f"CartoDB also failed, using lower zoom: {e2}")
                        ctx.add_basemap(ax[pos], zoom=max(1, zoom_level-2), 
                                      source=ctx.providers.CartoDB.Positron, crs=projection)
            
            # Add gridlines with labels only on left and bottom to avoid overlap
            gridliner = ax[pos].gridlines(draw_labels=True, linestyle="--", linewidth=0.5, color='gray')
            gridliner.top_labels = False
            gridliner.right_labels = False
            gridliner.xformatter = LONGITUDE_FORMATTER
            gridliner.yformatter = LATITUDE_FORMATTER
            gridliner.xlabel_style = {'size': 10, 'color': 'black'}
            gridliner.ylabel_style = {'size': 10, 'color': 'black'}
        
        # Add overall figure title
        fig.suptitle(f'PICASSO PPI Map - Elevation: {ppi_el:.1f}°\n{instrument_name}', 
                     fontsize=16, fontweight='bold', y=0.98)
        
        # Save figure
        if zoom_km is not None:
            figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime_sweep_str}_ppi_el{ppi_el:.2f}_map_zoom{zoom_km:.0f}km_l1_{product_version}.png'
        else:
            figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime_sweep_str}_ppi_el{ppi_el:.2f}_map_l1_{product_version}.png'
        plt.savefig(os.path.join(ppi_map_figpath, figname), dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Saved: {figname}")
    
    print(f"PPI map plots saved to: {ppi_map_figpath}")

def _plot_ppi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max, max_range):
    """Helper function to plot PPI fields on axes for turning scans."""
    
    radar = display._radar
    
    # Reflectivity (DBZ)
    display.plot_ppi("DBZ", ax=axes[0,0], sweep=sweep_idx, vmin=-60, vmax=40,
                     cmap=COLORMAPS['dbz'], colorbar_flag=False,
                     gatefilter=gatefilter)
    _setup_ppi_axes(axes[0,0], max_range)
    dbz_title = radar.fields['DBZ'].get('long_name', 'Reflectivity')
    axes[0,0].set_title(dbz_title)
    if len(axes[0,0].collections) > 0:
        plt.colorbar(axes[0,0].collections[-1], ax=axes[0,0], orientation='vertical',
                     label='Reflectivity (dBZ)', pad=0.02, aspect=20)
    
    # Velocity (VEL)  
    display.plot_ppi("VEL", ax=axes[1,0], sweep=sweep_idx, vmin=vel_min, vmax=vel_max,
                     cmap=COLORMAPS['vel'], colorbar_flag=False,
                     gatefilter=gatefilter)
    _setup_ppi_axes(axes[1,0], max_range)
    vel_title = radar.fields['VEL'].get('long_name', 'Doppler Velocity')
    axes[1,0].set_title(vel_title)
    if len(axes[1,0].collections) > 0:
        plt.colorbar(axes[1,0].collections[-1], ax=axes[1,0], orientation='vertical',
                     label='Velocity (m/s)', pad=0.02, aspect=20)
    
    # Spectrum width (WIDTH)
    display.plot_ppi("WIDTH", ax=axes[1,1], sweep=sweep_idx,
                     norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                     cmap=COLORMAPS['width'], colorbar_flag=False,
                     gatefilter=gatefilter)
    _setup_ppi_axes(axes[1,1], max_range)
    width_title = radar.fields['WIDTH'].get('long_name', 'Spectrum Width')
    axes[1,1].set_title(width_title)
    if len(axes[1,1].collections) > 0:
        plt.colorbar(axes[1,1].collections[-1], ax=axes[1,1], orientation='vertical',
                     label='Spectrum Width (m/s)', pad=0.02, aspect=20)
    
    # Linear depolarization ratio (LDR)
    display.plot_ppi("LDR", ax=axes[0,1], sweep=sweep_idx, vmin=-35, vmax=5,
                     cmap=COLORMAPS['ldr'], colorbar_flag=False,
                     gatefilter=gatefilter)
    _setup_ppi_axes(axes[0,1], max_range)
    ldr_title = radar.fields['LDR'].get('long_name', 'Linear Depolarization Ratio')
    axes[0,1].set_title(ldr_title)
    if len(axes[0,1].collections) > 0:
        plt.colorbar(axes[0,1].collections[-1], ax=axes[0,1], orientation='vertical',
                     label='LDR (dB)', pad=0.02, aspect=20)

def _setup_ppi_axes(ax, max_range):
    """Helper function to set up PPI plot axes."""
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

def _plot_sweep_angles(axes, sweep_times, elevation_data, azimuth_data, mode, 
                       radar_ds, start_ray_idx, end_ray_idx, mean_elevation=None):
    """Plot elevation and azimuth angles vs time for a sweep.
    
    Points marked as antenna_transition are plotted in grey.
    
    Args:
        mean_elevation: Actual mean elevation angle (used for auto-scaling y-axis)
    """
    
    # Get antenna transition flags for this sweep
    transition_flags = None
    if hasattr(radar_ds, 'antenna_transition') and radar_ds.antenna_transition is not None:
        if 'data' in radar_ds.antenna_transition:
            transition_flags = radar_ds.antenna_transition['data'][start_ray_idx:end_ray_idx+1]
    
    # Helper function to plot contiguous segments
    def plot_segments(ax, times, data, flags, normal_color, normal_label):
        """Plot data with different colors for normal vs transition segments."""
        if flags is None:
            ax.plot(times, data, color=normal_color, linewidth=1.5)
            return False
        
        is_transition = flags == 1
        has_transitions = np.any(is_transition)
        
        # Find contiguous segments
        # Add sentinels to detect all transitions
        padded = np.concatenate([[False], is_transition, [False]])
        diff = np.diff(padded.astype(int))
        
        # Find starts and ends of transition segments
        trans_starts = np.where(diff == 1)[0]
        trans_ends = np.where(diff == -1)[0]
        
        # Plot segments
        last_end = 0
        for start, end in zip(trans_starts, trans_ends):
            # Plot normal segment before transition (if any)
            if last_end < start:
                label = normal_label if last_end == 0 else None
                ax.plot(times[last_end:start+1], data[last_end:start+1], 
                       color=normal_color, linewidth=1.5, label=label)
            
            # Plot transition segment
            label = 'Antenna transition' if start == trans_starts[0] else None
            ax.plot(times[start:end], data[start:end], 
                   color='grey', linewidth=1.5, alpha=0.7, label=label)
            
            last_end = end
        
        # Plot final normal segment (if any)
        if last_end < len(times):
            label = normal_label if len(trans_starts) == 0 else None
            ax.plot(times[last_end:], data[last_end:], 
                   color=normal_color, linewidth=1.5, label=label)
        
        return has_transitions
    
    # Elevation angle panel (bottom left)
    has_transition = plot_segments(axes[2,0], sweep_times, elevation_data, 
                                     transition_flags, 'b', 'Normal')
    
    if has_transition:
        axes[2,0].legend(loc='best', fontsize=8)
    
    axes[2,0].set_xlabel('Time (UTC)')
    axes[2,0].set_ylabel('Elevation (°)')
    axes[2,0].set_title('Elevation Angle')
    axes[2,0].grid(True, alpha=0.3)
    axes[2,0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.setp(axes[2,0].xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Set aspect ratio to match radar field panels better (make it shorter)
    axes[2,0].set_box_aspect(0.4)
    
    # Set elevation limits based on actual mean elevation
    if mean_elevation is not None and mean_elevation > 85.0:
        # Near-vertical - zoom in on high elevations
        axes[2,0].set_ylim(80, 95)
    else:
        # RHI or turning - show full range
        axes[2,0].set_ylim(0, 95)
    
    # Azimuth angle panel (bottom right)
    has_transition = plot_segments(axes[2,1], sweep_times, azimuth_data, 
                                     transition_flags, 'r', 'Normal')
    
    if has_transition:
        axes[2,1].legend(loc='best', fontsize=8)
    
    axes[2,1].set_xlabel('Time (UTC)')
    axes[2,1].set_ylabel('Azimuth (°)')
    axes[2,1].set_title('Azimuth Angle')
    axes[2,1].grid(True, alpha=0.3)
    axes[2,1].set_ylim(0, 360)
    axes[2,1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.setp(axes[2,1].xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Set aspect ratio to match radar field panels better (make it shorter)
    axes[2,1].set_box_aspect(0.4)

def _export_antenna_transitions_to_csv(radar_ds, ncfile, output_dir):
    """
    Export antenna_transition indices to CSV file with scan rate information.
    
    For MAN scans, includes both elevation and azimuth scan rates.
    
    Args:
        radar_ds: PyART radar object
        ncfile: Path to the source NetCDF file
        output_dir: Directory to save CSV file
    """
    # Check if antenna_transition exists
    if not hasattr(radar_ds, 'antenna_transition') or radar_ds.antenna_transition is None:
        print(f"  No antenna_transition data found in {ncfile}")
        return
    
    if 'data' not in radar_ds.antenna_transition:
        print(f"  antenna_transition has no data field in {ncfile}")
        return
    
    transition_flags = radar_ds.antenna_transition['data']
    
    # Get time data for reference
    time_data = radar_ds.time['data']
    time_units = radar_ds.time['units']
    times = cftime.num2pydate(time_data, time_units)
    
    # Get elevation and azimuth angles
    elevation = radar_ds.elevation['data']
    azimuth = radar_ds.azimuth['data']
    
    # Check if this is a MAN scan with separate elevation and azimuth scan rates
    has_separate_rates = (hasattr(radar_ds, 'elevation_scan_rate') and 
                         radar_ds.elevation_scan_rate is not None and
                         hasattr(radar_ds, 'azimuth_scan_rate') and 
                         radar_ds.azimuth_scan_rate is not None)
    
    if has_separate_rates:
        # MAN scan: use both elevation and azimuth scan rates
        elevation_scan_rate = radar_ds.elevation_scan_rate['data']
        azimuth_scan_rate = radar_ds.azimuth_scan_rate['data']
        
        # Create DataFrame with separate rate columns
        df = pd.DataFrame({
            'ray_index': np.arange(len(transition_flags)),
            'antenna_transition': transition_flags,
            'timestamp': [t.strftime('%Y-%m-%dT%H:%M:%S.%fZ') for t in times],
            'elevation_deg': elevation,
            'azimuth_deg': azimuth,
            'elevation_scan_rate_deg_per_sec': elevation_scan_rate,
            'azimuth_scan_rate_deg_per_sec': azimuth_scan_rate
        })
        
        rate_columns = ['elevation_scan_rate_deg_per_sec', 'azimuth_scan_rate_deg_per_sec']
        print(f"  Using separate elevation and azimuth scan rates for MAN scan")
        
    else:
        # Standard RHI/PPI: use single scan rate
        if hasattr(radar_ds, 'scan_rate') and radar_ds.scan_rate is not None:
            if 'data' in radar_ds.scan_rate:
                scan_rate = radar_ds.scan_rate['data']
            else:
                print(f"  Warning: scan_rate has no data field, will calculate from angles")
                scan_rate = _calculate_scan_rate(elevation, times)
        else:
            print(f"  Warning: No scan_rate found in file, will calculate from angles")
            scan_rate = _calculate_scan_rate(elevation, times)
        
        # Create DataFrame with single rate column
        df = pd.DataFrame({
            'ray_index': np.arange(len(transition_flags)),
            'antenna_transition': transition_flags,
            'timestamp': [t.strftime('%Y-%m-%dT%H:%M:%S.%fZ') for t in times],
            'elevation_deg': elevation,
            'azimuth_deg': azimuth,
            'scan_rate_deg_per_sec': scan_rate
        })
        
        rate_columns = ['scan_rate_deg_per_sec']
    
    # Generate output filename based on input NetCDF file
    basename = os.path.splitext(os.path.basename(ncfile))[0]
    csv_filename = f"{basename}_antenna_transitions.csv"
    csv_filepath = os.path.join(output_dir, csv_filename)
    
    # Export to CSV
    df.to_csv(csv_filepath, index=False)
    print(f"  Antenna transitions exported to: {csv_filepath}")
    
    # Also create a summary file with just transition indices
    transition_indices = np.where(transition_flags == 1)[0]
    if len(transition_indices) > 0:
        summary_filename = f"{basename}_transition_indices.csv"
        summary_filepath = os.path.join(output_dir, summary_filename)
        
        # Build summary dataframe
        summary_data = {
            'ray_index': transition_indices,
            'timestamp': [times[i].strftime('%Y-%m-%dT%H:%M:%S.%fZ') for i in transition_indices],
            'elevation_deg': elevation[transition_indices],
            'azimuth_deg': azimuth[transition_indices]
        }
        
        # Add rate columns
        for col in rate_columns:
            summary_data[col] = df[col].values[transition_indices]
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(summary_filepath, index=False)
        print(f"  Transition indices summary: {summary_filepath}")
        print(f"    Total rays: {len(transition_flags)}, Transition rays: {len(transition_indices)} ({100*len(transition_indices)/len(transition_flags):.1f}%)")
    else:
        print(f"  No antenna transitions found (all flags are 0)")

def _calculate_scan_rate(elevation, times):
    """Calculate scan rate from elevation angles and timestamps (fallback)."""
    scan_rate = np.zeros(len(elevation))
    for i in range(len(elevation) - 1):
        dt = (times[i+1] - times[i]).total_seconds()
        if dt > 0:
            scan_rate[i] = (elevation[i+1] - elevation[i]) / dt
        else:
            scan_rate[i] = 0.0
    # Last point uses backward difference
    if len(elevation) > 1:
        dt = (times[-1] - times[-2]).total_seconds()
        if dt > 0:
            scan_rate[-1] = (elevation[-1] - elevation[-2]) / dt
        else:
            scan_rate[-1] = 0.0
    return scan_rate

def make_picasso_man_plot(ncfile, figpath, blflag=False, skip_all_transition=False):
    """
    Create MAN (Manual tracking) quicklook plots for PICASSO aircraft tracking data.
    
    This creates special visualizations for MAN scans including:
    - Overview plot with angle timeseries and full-track RHI
    - Individual plots for each phase/sweep (upward, dwell, downward, turning)
    
    Args:
        ncfile: Path to CF-Radial NetCDF file
        figpath: Output directory for figure
        blflag: Boundary layer flag for plot limits
        skip_all_transition: If True, skip sweeps where 100% of rays are antenna_transition=1
    """
    print(f"Creating MAN plot for: {ncfile}")
    
    # Get scan information
    time_coverage_start = get_time_coverage_start(ncfile)
    if not time_coverage_start:
        print(f"Could not get time coverage from {ncfile}")
        return
        
    netcdf_time = datetime.datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%SZ")
    
    # Read radar data
    try:
        radar_ds = pyart.io.read_cfradial(ncfile)
    except KeyError as e:
        print(f"Error reading radar file - missing field: {e}")
        print(f"  Checking file structure...")
        try:
            with nc4.Dataset(ncfile, 'r') as ds:
                print(f"  Available variables: {list(ds.variables.keys())}")
                print(f"  Global attributes: {ds.ncattrs()}")
                if 'latitude' not in ds.variables and 'latitude' not in ds.ncattrs():
                    print(f"  ERROR: latitude field is missing from file!")
                if 'longitude' not in ds.variables and 'longitude' not in ds.ncattrs():
                    print(f"  ERROR: longitude field is missing from file!")
        except Exception as check_error:
            print(f"  Could not diagnose file structure: {check_error}")
        return
    except Exception as e:
        print(f"Error reading radar file: {e}")
        return
    
    # Get product version, instrument name, and phase_sequence from NetCDF file
    try:
        with nc4.Dataset(ncfile) as ds:
            product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
            instrument_name = ds.source if 'source' in ds.ncattrs() else 'ncas-mobile-ka-band-radar-1'
            # Read phase_sequence if available (for phase-split MAN files)
            if 'phase_sequence' in ds.ncattrs():
                phase_sequence_str = ds.phase_sequence
                print(f"  Found phase_sequence in NetCDF: {phase_sequence_str}")
                # Store it in radar_ds.metadata so it's accessible throughout the function
                radar_ds.metadata['phase_sequence'] = phase_sequence_str
    except:
        product_version = VERSION
        instrument_name = 'ncas-mobile-ka-band-radar-1'
    
    # Set up gate filter
    gatefilter = pyart.correct.GateFilter(radar_ds)
    gatefilter.exclude_below('SNR', -20)
    
    # Exclude rays flagged as antenna_transition (repositioning between sweeps)
    if hasattr(radar_ds, 'antenna_transition') and radar_ds.antenna_transition is not None:
        if 'data' in radar_ds.antenna_transition:
            transition_flags = radar_ds.antenna_transition['data']
            transition_indices = np.where(transition_flags == 1)[0]
            if len(transition_indices) > 0:
                # Exclude all gates in rays flagged as antenna_transition
                # Properly access the gate_excluded attribute (2D boolean array)
                for ray_idx in transition_indices:
                    if hasattr(gatefilter, 'gate_excluded'):
                        gatefilter.gate_excluded[ray_idx, :] = True
                    elif hasattr(gatefilter, '_gate_excluded'):
                        gatefilter._gate_excluded[ray_idx, :] = True
    
    # Set up output directory
    man_figpath = os.path.join(figpath, 'man', netcdf_time.strftime('%Y%m%d'))
    os.makedirs(man_figpath, exist_ok=True)
    
    # Export antenna transition indices to CSV (disabled for quicklook generation)
    # _export_antenna_transitions_to_csv(radar_ds, ncfile, man_figpath)
    
    dtime0 = cftime.num2pydate(radar_ds.time['data'][0], radar_ds.time['units'])
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S")
    
    # Create overview plot with angle timeseries
    _make_man_overview_plot(radar_ds, gatefilter, man_figpath, dtime0_str, product_version, instrument_name, blflag)
    
    # Create individual sweep/phase plots
    _make_man_sweep_plots(radar_ds, gatefilter, man_figpath, dtime0_str, product_version, instrument_name, blflag)
    
    print(f"MAN plots saved to: {man_figpath}")

def _make_man_overview_plot(radar_ds, gatefilter, figpath, dtime_str, version, instrument_name, blflag, skip_all_transition=False):
    """Create MAN scan overview plot with angle timeseries."""
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # Extract data
    time_data = radar_ds.time['data']
    time_units = radar_ds.time['units']
    times = cftime.num2pydate(time_data, time_units)
    
    azimuth = radar_ds.azimuth['data']
    elevation = radar_ds.elevation['data']
    
    # Get sweep information
    nsweeps = radar_ds.nsweeps
    
    # Get valid sweep indices (optionally filtering out all-transition sweeps)
    valid_sweep_indices = get_valid_sweep_indices(radar_ds, skip_all_transition=skip_all_transition)
    
    sweep_modes = []
    sweep_starts = []
    sweep_ends = []
    
    # Check if phase_sequence metadata exists (from phase-split MAN scans)
    phase_sequence = None
    if 'phase_sequence' in radar_ds.metadata:
        phase_list = radar_ds.metadata['phase_sequence'].split(', ')
        if len(phase_list) == nsweeps:
            phase_sequence = phase_list
            print(f"  Overview plot using phase_sequence from NetCDF: {phase_list}")
    
    for i in valid_sweep_indices:
        # Use phase_sequence if available, otherwise use sweep_mode
        if phase_sequence:
            mode = phase_sequence[i]
        else:
            # Safely access sweep_mode - it might be stored differently in PyART
            try:
                if hasattr(radar_ds, 'sweep_mode') and isinstance(radar_ds.sweep_mode, dict):
                    mode = radar_ds.sweep_mode['data'][i]
                else:
                    # Fallback: try to read from sweep metadata or use default
                    mode = 'manual_rhi'  # Default for MAN scans
                    
                # Decode byte string if necessary
                if isinstance(mode, bytes):
                    mode = mode.decode('utf-8').strip()
                elif isinstance(mode, np.ndarray):
                    mode = mode.tobytes().decode('utf-8').strip()
                else:
                    mode = str(mode).strip()
                
                # Remove null characters
                mode = mode.replace('\x00', '')
            except (KeyError, TypeError, AttributeError) as e:
                print(f"  Warning: Could not read sweep_mode for sweep {i}, using 'manual_rhi'")
                mode = 'manual_rhi'
        
        sweep_modes.append(mode)
        sweep_starts.append(radar_ds.sweep_start_ray_index['data'][i])
        sweep_ends.append(radar_ds.sweep_end_ray_index['data'][i])
    
    # Top left: Elevation vs Time
    ax1 = fig.add_subplot(gs[0, 0])
    for i in range(nsweeps):
        start_idx = sweep_starts[i]
        end_idx = sweep_ends[i]
        phase = sweep_modes[i]
        color = PHASE_COLORS.get(phase, PHASE_COLORS['other'])
        
        ax1.plot(times[start_idx:end_idx+1], elevation[start_idx:end_idx+1],
                color=color, linewidth=2, label=phase if i == 0 or phase != sweep_modes[i-1] else '')
    
    ax1.axhline(88, color='red', linestyle='--', alpha=0.5, label='88° threshold')
    ax1.set_xlabel('Time (UTC)')
    ax1.set_ylabel('Elevation (°)')
    ax1.set_title('Elevation Angle vs Time')
    ax1.grid(True, alpha=0.3)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Top right: Azimuth vs Time
    ax2 = fig.add_subplot(gs[0, 1])
    for i in range(nsweeps):
        start_idx = sweep_starts[i]
        end_idx = sweep_ends[i]
        phase = sweep_modes[i]
        color = PHASE_COLORS.get(phase, PHASE_COLORS['other'])
        
        ax2.plot(times[start_idx:end_idx+1], azimuth[start_idx:end_idx+1],
                color=color, linewidth=2, label=phase if i == 0 or phase != sweep_modes[i-1] else '')
    
    ax2.set_xlabel('Time (UTC)')
    ax2.set_ylabel('Azimuth (°)')
    ax2.set_title('Azimuth Angle vs Time')
    ax2.grid(True, alpha=0.3)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Create legend for phases
    handles = [plt.Line2D([0], [0], color=color, linewidth=2, label=phase) 
               for phase, color in PHASE_COLORS.items() if phase in sweep_modes]
    ax2.legend(handles=handles, loc='best', fontsize=8)
    
    # Middle: Reflectivity (all sweeps)
    display = pyart.graph.RadarDisplay(radar_ds)
    hmax, xmin, xmax = setup_plot_limits(blflag, 'man')
    
    ax3 = fig.add_subplot(gs[1, :])
    display.plot_rhi("DBZ", ax=ax3, sweep=0, vmin=-60, vmax=40,
                     cmap=COLORMAPS['dbz'], colorbar_label='Reflectivity (dBZ)',
                     gatefilter=gatefilter, filter_transitions=True)
    ax3.set_ylim(0, hmax)
    ax3.set_xlim(xmin, xmax)
    ax3.grid(True, alpha=0.3)
    ax3.set_title(f'Reflectivity - Full Track ({nsweeps} phase{"s" if nsweeps > 1 else ""})')
    
    # Bottom: Velocity (all sweeps)
    ax4 = fig.add_subplot(gs[2, :])
    vel_field = radar_ds.fields['VEL']
    vel_min = vel_field.get('field_limit_lower', -10)
    vel_max = vel_field.get('field_limit_upper', 10)
    
    display.plot_rhi("VEL", ax=ax4, sweep=0, vmin=vel_min, vmax=vel_max,
                     cmap=COLORMAPS['vel'], colorbar_label='Velocity (m/s)',
                     gatefilter=gatefilter, filter_transitions=True)
    ax4.set_ylim(0, hmax)
    ax4.set_xlim(xmin, xmax)
    ax4.grid(True, alpha=0.3)
    ax4.set_title('Doppler Velocity - Full Track')
    
    # Add main title
    fig.suptitle(f'PICASSO MAN Scan - Aircraft Tracking\n{instrument_name}\n{times[0].strftime("%Y-%m-%d %H:%M:%S")} UTC', 
                 fontsize=14, fontweight='bold')
    
    # Save
    figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime_str}_man_overview_l1_{version}.png'
    plt.savefig(os.path.join(figpath, figname), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Overview plot created: {figname}")

def _make_man_sweep_plots(radar_ds, gatefilter, figpath, dtime_str, version, instrument_name, blflag, skip_all_transition=False):
    """Create individual plots for each MAN scan phase/sweep."""
    
    display = pyart.graph.RadarDisplay(radar_ds)
    hmax, xmin, xmax = setup_plot_limits(blflag, 'man')
    
    # For PPI plots, use max_range based on boundary layer flag
    max_range = 20 if blflag else 30
    
    vel_field = radar_ds.fields['VEL']
    vel_min = vel_field.get('field_limit_lower', -10)
    vel_max = vel_field.get('field_limit_upper', 10)
    
    nsweeps = radar_ds.nsweeps
    
    # Get valid sweep indices (optionally filtering out all-transition sweeps)
    valid_sweep_indices = get_valid_sweep_indices(radar_ds, skip_all_transition=skip_all_transition)
    
    if skip_all_transition and len(valid_sweep_indices) < nsweeps:
        print(f"  Skipping {nsweeps - len(valid_sweep_indices)} sweep(s) that are 100% antenna transitions")
    
    # Check if phase_sequence metadata exists (from phase-split MAN scans)
    phase_sequence = None
    if 'phase_sequence' in radar_ds.metadata:
        phase_list = radar_ds.metadata['phase_sequence'].split(', ')
        if len(phase_list) == nsweeps:
            phase_sequence = phase_list
            print(f"  Using phase_sequence from NetCDF: {phase_list}")
        else:
            print(f"  Warning: phase_sequence length ({len(phase_list)}) doesn't match nsweeps ({nsweeps})")
    
    for sweep_idx in valid_sweep_indices:
        # Get sweep phase info - use phase_sequence if available
        if phase_sequence:
            mode = phase_sequence[sweep_idx]
            print(f"  Sweep {sweep_idx+1}: Using phase from phase_sequence: '{mode}'")
        else:
            # Safely access sweep_mode - it might be stored differently in PyART
            try:
                if hasattr(radar_ds, 'sweep_mode') and isinstance(radar_ds.sweep_mode, dict):
                    mode = radar_ds.sweep_mode['data'][sweep_idx]
                else:
                    # Fallback: try to read from sweep metadata or use default
                    mode = 'manual_rhi'  # Default for MAN scans
                    
                if isinstance(mode, bytes):
                    mode = mode.decode('utf-8').strip()
                elif isinstance(mode, np.ndarray):
                    mode = mode.tobytes().decode('utf-8').strip()
                else:
                    mode = str(mode).strip()
                mode = mode.replace('\x00', '')
            except (KeyError, TypeError, AttributeError) as e:
                print(f"  Warning: Could not read sweep_mode for sweep {sweep_idx}, using 'manual_rhi'")
                mode = 'manual_rhi'
        
        # Clean mode string (remove null bytes, strip whitespace)
        mode = str(mode).replace('\x00', '').strip().lower()
        
        # Get sweep start and end times
        start_ray_idx = radar_ds.sweep_start_ray_index['data'][sweep_idx]
        end_ray_idx = radar_ds.sweep_end_ray_index['data'][sweep_idx]
        
        # Check if sweep has enough rays to plot
        n_rays = end_ray_idx - start_ray_idx + 1
        if n_rays < 2:
            print(f"  Sweep {sweep_idx+1}: Skipping - insufficient rays ({n_rays})")
            continue
        
        start_time = cftime.num2pydate(radar_ds.time['data'][start_ray_idx], radar_ds.time['units'])
        end_time = cftime.num2pydate(radar_ds.time['data'][end_ray_idx], radar_ds.time['units'])
        
        time_range_str = f"{start_time.strftime('%H:%M:%S')} - {end_time.strftime('%H:%M:%S')} UTC"
        
        # Calculate actual mean azimuth/elevation for this sweep
        azimuth_data = radar_ds.azimuth['data'][start_ray_idx:end_ray_idx+1]
        elevation_data = radar_ds.elevation['data'][start_ray_idx:end_ray_idx+1]
        
        # Calculate circular mean for azimuth (handles 0/360 wrapping)
        mean_azimuth = np.degrees(np.arctan2(
            np.mean(np.sin(np.radians(azimuth_data))),
            np.mean(np.cos(np.radians(azimuth_data)))
        )) % 360
        
        # Regular mean for elevation
        mean_elevation = np.mean(elevation_data)
        
        # Create figure with 3x2 layout (4 data fields + 2 angle panels)
        fig, axes = plt.subplots(3, 2, figsize=(14, 14))
        
        # Get time array for this sweep
        sweep_times = cftime.num2pydate(
            radar_ds.time['data'][start_ray_idx:end_ray_idx+1], 
            radar_ds.time['units']
        )
        
        # Plot fields - use appropriate visualization based on phase type
        print(f"  Sweep {sweep_idx+1}: mode='{mode}', mean_elevation={mean_elevation:.1f}°")
        
        try:
            if mode == 'turning':
                _plot_ppi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max, max_range)
                subtitle = f'Mean Elevation: {mean_elevation:.1f}° | {time_range_str}'
            elif mode in ['dwell', 'vertical_pointing'] or mean_elevation > 85.0:
                # Use VPT-style plotting for:
                # 1. Dwell and vertical_pointing modes
                # 2. ANY sweep with mean elevation > 85° (nearly vertical)
                # These show nothing useful in RHI format and should be time-height plots
                if mean_elevation > 85.0 and mode not in ['dwell', 'vertical_pointing']:
                    print(f"    High elevation sweep ({mean_elevation:.1f}°), using VPT format")
                else:
                    print(f"    Using VPT format for {mode} phase")
                _plot_vpt_fields_single_sweep(display, axes, sweep_idx, gatefilter, vel_min, vel_max, hmax, skip_all_transition)
                subtitle = f'Mean Elevation: {mean_elevation:.1f}° | {time_range_str}'
            else:
                _plot_rhi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max, hmax, xmin, xmax)
                subtitle = f'Mean Azimuth: {mean_azimuth:.1f}° | {time_range_str}'
            
            # Add elevation and azimuth panels
            _plot_sweep_angles(axes, sweep_times, elevation_data, azimuth_data, mode,
                              radar_ds, start_ray_idx, end_ray_idx, mean_elevation)
            
            # Add phase information to title (use original case for display)
            mode_display = mode.replace('_', ' ').title().replace(' ', '_')
            phase_color = PHASE_COLORS.get(mode, PHASE_COLORS.get('other', '#808080'))
            fig.suptitle(f'MAN Scan - Phase: {mode_display.upper()} (Sweep {sweep_idx+1}/{nsweeps})\n{instrument_name}\n{subtitle}', 
                        fontsize=12, fontweight='bold',
                        color=phase_color)
            
            # Save
            figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime_str}_man_sweep{sweep_idx:02d}_{mode.replace("_", "-")}_l1_{version}.png'
            plt.tight_layout()
            plt.savefig(os.path.join(figpath, figname), dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"  Phase plot created: sweep {sweep_idx+1}/{nsweeps} ({mode})")
        except Exception as e:
            print(f"  Warning: Failed to plot sweep {sweep_idx+1}: {e}")
            plt.close()
            continue


def _plot_vpt_fields_single_sweep(display, axes, sweep_idx, gatefilter, vel_min, vel_max, hmax, skip_all_transition=False):
    """Helper function to plot a single sweep in VPT (time-height) format.
    
    Used for dwell and vertical_pointing phases where RHI format shows nothing useful.
    Manually extracts and plots data to avoid PyART dimension issues with single sweeps.
    For near-vertical pointing (>85°), y-axis is "Height (km)". Otherwise "Range (km)".
    
    Args:
        skip_all_transition: If True, completely exclude antenna_transition rays from the plot.
                           If False, include all rays and use gatefilter for masking.
    """
    radar = display._radar
    
    # Get sweep ray indices
    start_ray = radar.sweep_start_ray_index['data'][sweep_idx]
    end_ray = radar.sweep_end_ray_index['data'][sweep_idx]
    
    # Filter out antenna_transition rays for VPT plots only if skip_all_transition is True
    sweep_ray_indices = np.arange(start_ray, end_ray + 1)
    if skip_all_transition and hasattr(radar, 'antenna_transition') and radar.antenna_transition is not None:
        if 'data' in radar.antenna_transition:
            transition_flags = radar.antenna_transition['data'][start_ray:end_ray+1]
            # Keep only non-transition rays (antenna_transition == 0)
            valid_mask = (transition_flags == 0)
            sweep_ray_indices = sweep_ray_indices[valid_mask]
            
            if len(sweep_ray_indices) == 0:
                print(f"    Warning: All rays in sweep {sweep_idx} are antenna_transition, cannot plot")
                return
    
    # Check if we have enough rays to plot
    if len(sweep_ray_indices) < 2:
        print(f"    Warning: Not enough rays in sweep {sweep_idx} to plot ({len(sweep_ray_indices)} rays)")
        return
    
    # Calculate mean elevation to determine if truly vertical
    mean_elevation = np.mean(radar.elevation['data'][sweep_ray_indices])
    is_vertical = mean_elevation > 85.0  # Consider >85° as vertical
    y_label = 'Height (km)' if is_vertical else 'Range (km)'
    
    # Extract time and height data for this sweep (only valid rays)
    sweep_times = cftime.num2pydate(radar.time['data'][sweep_ray_indices], radar.time['units'])
    sweep_range = radar.range['data'] / 1000.0  # Convert to km
    
    # Check if we have gates to plot
    if len(sweep_range) == 0:
        print(f"    Warning: No range gates in sweep {sweep_idx}, cannot plot")
        return
    
    # Create meshgrid for plotting
    time_nums = mdates.date2num(sweep_times)
    X, Y = np.meshgrid(time_nums, sweep_range)
    
    # Extract field long names from metadata for titles
    dbz_long = radar.fields['DBZ'].get('long_name', 'Reflectivity')
    vel_long = radar.fields['VEL'].get('long_name', 'Velocity')
    width_long = radar.fields['WIDTH'].get('long_name', 'Spectrum Width')
    ldr_long = radar.fields['LDR'].get('long_name', 'LDR')
    
    # Plot Reflectivity (DBZ)
    field_data = radar.fields['DBZ']['data'][sweep_ray_indices, :].T
    if gatefilter is not None:
        field_mask = gatefilter.gate_excluded[sweep_ray_indices, :].T
        field_data = np.ma.masked_where(field_mask, field_data)
    
    pcm = axes[0,0].pcolormesh(X, Y, field_data, vmin=-60, vmax=40, 
                               cmap=COLORMAPS['dbz'], shading='nearest')
    plt.colorbar(pcm, ax=axes[0,0], label='Reflectivity (dBZ)', 
                 orientation='horizontal', pad=0.18, aspect=30)
    axes[0,0].set_ylim(0, hmax)
    axes[0,0].set_title(dbz_long)
    axes[0,0].set_xlabel('Time (UTC)')
    axes[0,0].set_ylabel(y_label)
    axes[0,0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Plot Velocity (VEL)
    field_data = radar.fields['VEL']['data'][sweep_ray_indices, :].T
    if gatefilter is not None:
        field_mask = gatefilter.gate_excluded[sweep_ray_indices, :].T
        field_data = np.ma.masked_where(field_mask, field_data)
    
    pcm = axes[1,0].pcolormesh(X, Y, field_data, vmin=vel_min, vmax=vel_max,
                               cmap=COLORMAPS['vel'], shading='nearest')
    plt.colorbar(pcm, ax=axes[1,0], label='Velocity (m/s)',
                 orientation='horizontal', pad=0.18, aspect=30)
    axes[1,0].set_ylim(0, hmax)
    axes[1,0].set_title(vel_long)
    axes[1,0].set_xlabel('Time (UTC)')
    axes[1,0].set_ylabel(y_label)
    axes[1,0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Plot Spectrum width (WIDTH)
    field_data = radar.fields['WIDTH']['data'][sweep_ray_indices, :].T
    if gatefilter is not None:
        field_mask = gatefilter.gate_excluded[sweep_ray_indices, :].T
        field_data = np.ma.masked_where(field_mask, field_data)
    
    pcm = axes[1,1].pcolormesh(X, Y, field_data, 
                               norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                               cmap=COLORMAPS['width'], shading='nearest')
    plt.colorbar(pcm, ax=axes[1,1], label='Spectrum Width (m/s)',
                 orientation='horizontal', pad=0.18, aspect=30)
    axes[1,1].set_ylim(0, hmax)
    axes[1,1].set_title(width_long)
    axes[1,1].set_xlabel('Time (UTC)')
    axes[1,1].set_ylabel(y_label)
    axes[1,1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Plot Linear depolarization ratio (LDR)
    field_data = radar.fields['LDR']['data'][sweep_ray_indices, :].T
    if gatefilter is not None:
        field_mask = gatefilter.gate_excluded[sweep_ray_indices, :].T
        field_data = np.ma.masked_where(field_mask, field_data)
    
    pcm = axes[0,1].pcolormesh(X, Y, field_data, vmin=-35, vmax=5,
                               cmap=COLORMAPS['ldr'], shading='nearest')
    plt.colorbar(pcm, ax=axes[0,1], label='LDR (dB)',
                 orientation='horizontal', pad=0.18, aspect=30)
    axes[0,1].set_ylim(0, hmax)
    axes[0,1].set_title(ldr_long)
    axes[0,1].set_xlabel('Time (UTC)')
    axes[0,1].set_ylabel(y_label)
    axes[0,1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))


def make_picasso_vpt_plot(ncfile, figpath, blflag=False, skip_all_transition=False):
    """
    Create VPT (Vertical Pointing) quicklook plot for PICASSO data.
    
    Args:
        ncfile: Path to CF-Radial NetCDF file
        figpath: Output directory for figure
        blflag: Boundary layer flag for plot limits
        skip_all_transition: If True, skip sweeps where 100% of rays are antenna_transition=1
    """
    print(f"Creating VPT plot for: {ncfile}")
    
    # Get scan information
    time_coverage_start = get_time_coverage_start(ncfile)
    if not time_coverage_start:
        print(f"Could not get time coverage from {ncfile}")
        return
    
    # Read radar data
    try:
        radar_ds = pyart.io.read_cfradial(ncfile)
    except Exception as e:
        print(f"Error reading radar file: {e}")
        return
    
    # Check if we should skip this file based on antenna_transition filtering
    if skip_all_transition:
        valid_sweep_indices = get_valid_sweep_indices(radar_ds, skip_all_transition=True)
        if len(valid_sweep_indices) == 0:
            print(f"  Skipping VPT file - all sweeps are 100% antenna transitions")
            return
    
    # Get product version and phase_sequence from NetCDF file
    try:
        with nc4.Dataset(ncfile) as ds:
            product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
            instrument_name = ds.source if 'source' in ds.ncattrs() else 'ncas-mobile-ka-band-radar-1'
            # Read phase_sequence if available (for phase-split files)
            if 'phase_sequence' in ds.ncattrs():
                phase_sequence_str = ds.phase_sequence
                print(f"  Found phase_sequence in NetCDF: {phase_sequence_str}")
                # Store it in radar_ds.metadata so it's accessible throughout the function
                radar_ds.metadata['phase_sequence'] = phase_sequence_str
    except:
        product_version = VERSION
        instrument_name = 'ncas-mobile-ka-band-radar-1'
    
    # Set up output directory
    netcdf_time = datetime.datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%SZ")
    vpt_figpath = os.path.join(figpath, 'vpt')
    os.makedirs(vpt_figpath, exist_ok=True)
    
    # Set height limits
    hmax = 4 if blflag else 12
    
    # Create time-height plots
    dtime0 = cftime.num2pydate(radar_ds.time['data'][0], radar_ds.time['units'])
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S")
    
    # Split radar data at large time gaps to prevent pixel stretching
    radar_segments = split_radar_at_time_gaps(radar_ds, max_gap_seconds=600)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    
    # Get velocity limits
    vel_field = radar_ds.fields['VEL']
    vel_min = vel_field.get('field_limit_lower', -10)
    vel_max = vel_field.get('field_limit_upper', 10)
    
    # Plot each segment
    for seg_idx, radar_segment in enumerate(radar_segments):
        display = pyart.graph.RadarDisplay(radar_segment)
        times = cftime.num2pydate(radar_segment.time['data'], radar_segment.time['units'])
        
        # Extract field labels from metadata
        if seg_idx == 0:
            dbz_long = radar_segment.fields['DBZ'].get('long_name', 'Reflectivity')
            dbz_units = radar_segment.fields['DBZ'].get('units', 'dBZ')
            dbz_title = f"{dbz_long} ({dbz_units})"
            
            vel_long = radar_segment.fields['VEL'].get('long_name', 'Velocity')
            vel_units = radar_segment.fields['VEL'].get('units', 'm/s')
            vel_title = f"{vel_long} ({vel_units})"
            
            width_long = radar_segment.fields['WIDTH'].get('long_name', 'Spectrum Width')
            width_units = radar_segment.fields['WIDTH'].get('units', 'm/s')
            width_title = f"{width_long} ({width_units})"
            
            ldr_long = radar_segment.fields['LDR'].get('long_name', 'LDR')
            ldr_units = radar_segment.fields['LDR'].get('units', 'dB')
            ldr_title = f"{ldr_long} ({ldr_units})"
        
        # First segment - create plots with colorbars
        if seg_idx == 0:
            _plot_vpt_field(display, axes[0,0], 'DBZ', dbz_title, 'Reflectivity (dBZ)',
                           COLORMAPS['dbz'], -60, 40, hmax, times)
            _plot_vpt_field(display, axes[0,1], 'VEL', vel_title, 'Velocity (m/s)',
                           COLORMAPS['vel'], vel_min, vel_max, hmax, times)
            _plot_vpt_field(display, axes[1,0], 'WIDTH', width_title, 'Spectrum Width (m/s)',
                           COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True)
            _plot_vpt_field(display, axes[1,1], 'LDR', ldr_title, 'LDR (dB)',
                           COLORMAPS['ldr'], -35, 5, hmax, times)
        else:
            # Subsequent segments - overlay without colorbars
            _plot_vpt_field_overlay(display, axes[0,0], 'DBZ', COLORMAPS['dbz'], -60, 40, hmax, times)
            _plot_vpt_field_overlay(display, axes[0,1], 'VEL', COLORMAPS['vel'], vel_min, vel_max, hmax, times)
            _plot_vpt_field_overlay(display, axes[1,0], 'WIDTH', COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True)
            _plot_vpt_field_overlay(display, axes[1,1], 'LDR', COLORMAPS['ldr'], -35, 5, hmax, times)
    
    fig.suptitle(f'PICASSO VPT - {dtime0.strftime("%Y-%m-%d")}\n{instrument_name}', fontsize=14, fontweight='bold')
    
    # Save
    figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime0_str}_vpt_l1_{product_version}.png'
    plt.tight_layout()
    plt.savefig(os.path.join(vpt_figpath, figname), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"VPT plot saved to: {vpt_figpath}")

def _plot_vpt_field(display, ax, field_name, title, colorbar_label, cmap, vmin, vmax, hmax, times, log_scale=False):
    """Helper function to plot VPT field."""
    if log_scale:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    
    # Plot without built-in colorbar so we can add horizontal one
    pcm = display.plot_vpt(field_name, ax=ax, time_axis_flag=True, edges=True,
                           vmin=vmin, vmax=vmax, cmap=cmap, 
                           colorbar_flag=False, norm=norm, shading='auto')
    
    ax.set_ylim(0, hmax)
    ax.set_xlabel('Time (UTC)')
    ax.set_ylabel('Height (km)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    
    # Add horizontal colorbar to match composite VPT style
    plt.colorbar(pcm, ax=ax, label=colorbar_label, orientation='horizontal', pad=0.15, aspect=30)

def _plot_vpt_field_overlay(display, ax, field_name, cmap, vmin, vmax, hmax, times, log_scale=False):
    """Helper function to plot VPT field overlay (no colorbar)."""
    if log_scale:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    
    display.plot_vpt(field_name, ax=ax, time_axis_flag=True, edges=True,
                     vmin=vmin, vmax=vmax, cmap=cmap, 
                     colorbar_flag=False, norm=norm, shading='auto')

def make_picasso_vpt_composite(datestr, inpath, figpath, blflag=False, skip_all_transition=False):
    """
    Create composite VPT plot combining VPT sweeps and vertical pointing MAN dwell segments for a day.
    
    This function merges all vertical pointing data from:
    1. Dedicated VPT scan files
    2. Vertical pointing dwell segments (elevation > 85°) from MAN (manual tracking) scans
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Output directory path
        blflag: Boundary layer flag for plot limits
        skip_all_transition: If True, skip sweeps where 100% of rays are antenna_transition=1
    """
    print(f"\nCreating composite VPT plot (VPT + vertical pointing MAN dwells) for {datestr}")
    
    date_path = os.path.join(inpath, datestr)
    if not os.path.exists(date_path):
        print(f"Date directory does not exist: {date_path}")
        return
    
    # Collect metadata about VPT segments (not the data itself to save memory)
    vpt_segment_info = []
    
    # 1. Get VPT file segment metadata
    vpt_files = sorted(glob.glob(os.path.join(date_path, '*_vpt_*.nc')))
    print(f"Found {len(vpt_files)} VPT files")
    
    for ncfile in vpt_files:
        try:
            # Just read times to identify segments, don't keep full radar data
            radar_ds = pyart.io.read_cfradial(ncfile)
            
            # Check if we should skip this file based on antenna_transition filtering
            if skip_all_transition:
                valid_sweep_indices = get_valid_sweep_indices(radar_ds, skip_all_transition=True)
                if len(valid_sweep_indices) == 0:
                    print(f"  Skipping VPT file (all sweeps are 100% antenna transitions): {os.path.basename(ncfile)}")
                    del radar_ds
                    continue
            
            times = cftime.num2pydate(radar_ds.time['data'], radar_ds.time['units'])
            
            # Find time gaps
            time_diffs = np.diff([t.timestamp() for t in times])
            gap_indices = np.where(time_diffs > 600)[0]
            
            # Store segment info
            if len(gap_indices) == 0:
                vpt_segment_info.append({
                    'file': ncfile,
                    'type': 'vpt_full',
                    'start_time': times[0],
                    'source': 'VPT'
                })
                print(f"  Indexed VPT file: {os.path.basename(ncfile)} (1 continuous segment)")
            else:
                # Multiple segments - split at time gaps
                num_segments = len(gap_indices) + 1
                print(f"  Indexed VPT file: {os.path.basename(ncfile)} (split into {num_segments} segments at time gaps > 10 min)")
                start_ray = 0
                for gap_idx in gap_indices:
                    vpt_segment_info.append({
                        'file': ncfile,
                        'type': 'vpt_partial',
                        'ray_start': start_ray,
                        'ray_end': gap_idx,
                        'start_time': times[start_ray],
                        'source': 'VPT'
                    })
                    start_ray = gap_idx + 1
                
                # Final segment
                if start_ray < len(times):
                    vpt_segment_info.append({
                        'file': ncfile,
                        'type': 'vpt_partial',
                        'ray_start': start_ray,
                        'ray_end': len(times) - 1,
                        'start_time': times[start_ray],
                        'source': 'VPT'
                    })
            
            del radar_ds  # Free memory immediately
            
        except Exception as e:
            print(f"  Error indexing VPT file {ncfile}: {e}")
            continue
    
    # 2. Get MAN vertical pointing dwell segment metadata
    man_files = sorted(glob.glob(os.path.join(date_path, '*_man_*.nc')))
    print(f"Found {len(man_files)} MAN files to search for vertical pointing dwells")
    
    for ncfile in man_files:
        try:
            radar_ds = pyart.io.read_cfradial(ncfile)
            
            # Check if this file has phase_sequence metadata
            if 'phase_sequence' not in radar_ds.metadata:
                del radar_ds
                continue
            
            phase_list = radar_ds.metadata['phase_sequence'].split(', ')
            nsweeps = radar_ds.nsweeps
            
            # Get valid sweep indices (filter out all-transition sweeps if requested)
            valid_sweep_indices = get_valid_sweep_indices(radar_ds, skip_all_transition=skip_all_transition)
            
            # Find dwell sweeps that are vertical pointing
            for sweep_idx in range(nsweeps):
                if sweep_idx >= len(phase_list):
                    continue
                
                # Skip if this sweep is filtered out
                if sweep_idx not in valid_sweep_indices:
                    continue
                    
                phase = phase_list[sweep_idx]
                
                # Only include dwell or vertical_pointing phases
                if phase in ['dwell', 'vertical_pointing']:
                    start_ray_idx = radar_ds.sweep_start_ray_index['data'][sweep_idx]
                    end_ray_idx = radar_ds.sweep_end_ray_index['data'][sweep_idx]
                    
                    # Check mean elevation - only include if vertical pointing (>85°)
                    elevation_data = radar_ds.elevation['data'][start_ray_idx:end_ray_idx+1]
                    mean_elevation = np.mean(elevation_data)
                    
                    if mean_elevation <= 85.0:
                        print(f"  Skipping {phase} in {os.path.basename(ncfile)}, sweep {sweep_idx+1}: mean elevation {mean_elevation:.1f}° (not vertical pointing)")
                        continue
                    
                    times = cftime.num2pydate(radar_ds.time['data'][start_ray_idx:end_ray_idx+1], 
                                             radar_ds.time['units'])
                    
                    vpt_segment_info.append({
                        'file': ncfile,
                        'type': 'man_dwell',
                        'sweep_idx': sweep_idx,
                        'start_time': times[0],
                        'source': f'MAN_{phase}'
                    })
                    
                    print(f"  Found vertical pointing {phase} in {os.path.basename(ncfile)}, sweep {sweep_idx+1} (elevation {mean_elevation:.1f}°)")
            
            del radar_ds  # Free memory
                    
        except KeyError as e:
            # Likely missing latitude/longitude - these are old files that need reprocessing
            print(f"  Skipping MAN file {os.path.basename(ncfile)}: missing required field {e} (needs reprocessing)")
            continue
        except Exception as e:
            print(f"  Error indexing MAN file {ncfile}: {e}")
            continue
    
    if not vpt_segment_info:
        print(f"No VPT data found for {datestr}")
        return
    
    print(f"\nTotal VPT segments: {len(vpt_segment_info)}")
    
    # Sort all segments by start time
    vpt_segment_info.sort(key=lambda x: x['start_time'])
    
    # Get product version and instrument name
    try:
        first_file = vpt_files[0] if vpt_files else man_files[0]
        with nc4.Dataset(first_file) as ds:
            product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
            instrument_name = ds.source if 'source' in ds.ncattrs() else 'ncas-mobile-ka-band-radar-1'
    except:
        product_version = VERSION
        instrument_name = 'ncas-mobile-ka-band-radar-1'
    
    # Set up output directory
    vpt_figpath = os.path.join(figpath, 'vpt')
    os.makedirs(vpt_figpath, exist_ok=True)
    
    # Set height limits
    hmax = 4 if blflag else 12
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    
    # Track whether we've plotted the first valid segment
    first_plotted = False
    
    # Plot each segment - load data on demand
    for idx, seg_info in enumerate(vpt_segment_info):
        try:
            # Load only this segment's data
            if seg_info['type'] == 'vpt_full':
                radar_ds = pyart.io.read_cfradial(seg_info['file'])
                radar_segment = split_radar_at_time_gaps(radar_ds, max_gap_seconds=600)[0]
                del radar_ds
                
            elif seg_info['type'] == 'vpt_partial':
                radar_ds = pyart.io.read_cfradial(seg_info['file'])
                # Extract rays for this segment
                ray_indices = np.arange(seg_info['ray_start'], seg_info['ray_end'] + 1)
                radar_segment = radar_ds.extract_sweeps([0])
                radar_segment.time['data'] = radar_ds.time['data'][ray_indices]
                for field_name in radar_ds.fields.keys():
                    radar_segment.fields[field_name]['data'] = radar_ds.fields[field_name]['data'][ray_indices, :]
                radar_segment.azimuth['data'] = radar_ds.azimuth['data'][ray_indices]
                radar_segment.elevation['data'] = radar_ds.elevation['data'][ray_indices]
                # Update metadata to reflect actual extracted dimensions
                radar_segment.nrays = len(ray_indices)
                radar_segment.sweep_end_ray_index['data'][0] = len(ray_indices) - 1
                if hasattr(radar_segment, 'antenna_transition') and radar_segment.antenna_transition is not None:
                    if 'data' in radar_segment.antenna_transition:
                        radar_segment.antenna_transition['data'] = radar_ds.antenna_transition['data'][ray_indices]
                del radar_ds
                
            elif seg_info['type'] == 'man_dwell':
                radar_ds = pyart.io.read_cfradial(seg_info['file'])
                dwell_sweep = radar_ds.extract_sweeps([seg_info['sweep_idx']])
                radar_segment = split_radar_at_time_gaps(dwell_sweep, max_gap_seconds=600)[0]
                del radar_ds, dwell_sweep
            
            # Validate segment has sufficient dimensions for plotting
            if radar_segment.nrays < 2 or radar_segment.ngates < 2:
                print(f"  Skipping segment {idx} with insufficient dimensions: nrays={radar_segment.nrays}, ngates={radar_segment.ngates}")
                del radar_segment
                continue
            
            # Validate that data dimensions match coordinate dimensions
            try:
                field_data = radar_segment.fields['DBZ']['data']
                time_len = len(radar_segment.time['data'])
                range_len = len(radar_segment.range['data'])
                
                if field_data.shape[0] != time_len or field_data.shape[1] != range_len:
                    print(f"  Skipping segment {idx} with mismatched field dimensions: field shape={field_data.shape}, time={time_len}, range={range_len}")
                    del radar_segment
                    continue
                
                if time_len < 2 or range_len < 2:
                    print(f"  Skipping segment {idx} with insufficient coordinate length: time={time_len}, range={range_len}")
                    del radar_segment
                    continue
            except Exception as e:
                print(f"  Skipping segment {idx}: could not validate dimensions - {e}")
                del radar_segment
                continue
            
            times = cftime.num2pydate(radar_segment.time['data'], radar_segment.time['units'])
            display = pyart.graph.RadarDisplay(radar_segment)
            
            # Extract field labels from metadata
            dbz_long = radar_segment.fields['DBZ'].get('long_name', 'Reflectivity')
            dbz_units = radar_segment.fields['DBZ'].get('units', 'dBZ')
            dbz_title = f"{dbz_long} ({dbz_units})"
            
            vel_long = radar_segment.fields['VEL'].get('long_name', 'Velocity')
            vel_units = radar_segment.fields['VEL'].get('units', 'm/s')
            vel_title = f"{vel_long} ({vel_units})"
            
            width_long = radar_segment.fields['WIDTH'].get('long_name', 'Spectrum Width')
            width_units = radar_segment.fields['WIDTH'].get('units', 'm/s')
            width_title = f"{width_long} ({width_units})"
            
            ldr_long = radar_segment.fields['LDR'].get('long_name', 'LDR')
            ldr_units = radar_segment.fields['LDR'].get('units', 'dB')
            ldr_title = f"{ldr_long} ({ldr_units})"
            
            if not first_plotted:
                # First valid segment - create plots with colorbars
                _plot_vpt_field_composite(display, axes[0,0], 'DBZ', dbz_title, 'Reflectivity (dBZ)',
                               COLORMAPS['dbz'], -60, 40, hmax, times, first=True)
                
                vel_field = radar_segment.fields['VEL']
                vel_min = vel_field.get('field_limit_lower', -10)
                vel_max = vel_field.get('field_limit_upper', 10)
                _plot_vpt_field_composite(display, axes[0,1], 'VEL', vel_title, 'Velocity (m/s)',
                               COLORMAPS['vel'], vel_min, vel_max, hmax, times, first=True)
                
                _plot_vpt_field_composite(display, axes[1,0], 'WIDTH', width_title, 'Spectrum Width (m/s)',
                               COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True, first=True)
                
                _plot_vpt_field_composite(display, axes[1,1], 'LDR', ldr_title, 'LDR (dB)',
                               COLORMAPS['ldr'], -35, 5, hmax, times, first=True)
                
                first_plotted = True
            else:
                # Subsequent segments - overlay without colorbars
                _plot_vpt_field_composite(display, axes[0,0], 'DBZ', dbz_title, 'Reflectivity (dBZ)',
                               COLORMAPS['dbz'], -60, 40, hmax, times, first=False)
                
                vel_field = radar_segment.fields['VEL']
                vel_min = vel_field.get('field_limit_lower', -10)
                vel_max = vel_field.get('field_limit_upper', 10)
                _plot_vpt_field_composite(display, axes[0,1], 'VEL', vel_title, 'Velocity (m/s)',
                               COLORMAPS['vel'], vel_min, vel_max, hmax, times, first=False)
                
                _plot_vpt_field_composite(display, axes[1,0], 'WIDTH', width_title, 'Spectrum Width (m/s)',
                               COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True, first=False)
                
                _plot_vpt_field_composite(display, axes[1,1], 'LDR', ldr_title, 'LDR (dB)',
                               COLORMAPS['ldr'], -35, 5, hmax, times, first=False)
            
            # Free memory immediately after plotting
            del radar_segment, display, times
            
        except Exception as e:
            print(f"  Error plotting segment {idx}: {e}")
            continue
    
    # Count sources
    vpt_count = sum(1 for s in vpt_segment_info if s['source'] == 'VPT')
    man_count = sum(1 for s in vpt_segment_info if s['source'].startswith('MAN_'))
    
    # Set title
    date_obj = datetime.datetime.strptime(datestr, '%Y%m%d')
    fig.suptitle(f'PICASSO Composite VPT - {date_obj.strftime("%Y-%m-%d")}\n'
                 f'{instrument_name}\n'
                 f'{len(vpt_segment_info)} segments (VPT: {vpt_count}, MAN vertical pointing: {man_count})', 
                 fontsize=14, fontweight='bold')
    
    # Save
    figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{datestr}_vpt_composite_l1_{product_version}.png'
    plt.tight_layout()
    plt.savefig(os.path.join(vpt_figpath, figname), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Composite VPT plot saved to: {vpt_figpath}")
    print(f"  Combined {len(vpt_segment_info)} segments: {vpt_count} VPT, {man_count} MAN vertical pointing")

def make_picasso_man_dwell_composite(datestr, inpath, figpath, blflag=False, skip_all_transition=False):
    """
    Create composite VPT-style plot from all vertical pointing dwell segments in MAN scans for a day.
    
    This function:
    1. Finds all MAN files for the specified date
    2. Extracts dwell and vertical_pointing segments with elevation > 85° from each file
    3. Combines them chronologically into a single time-height plot
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Output directory path
        blflag: Boundary layer flag for plot limits
        skip_all_transition: If True, skip sweeps where 100% of rays are antenna_transition=1
    """
    print(f"\nCreating MAN vertical pointing dwell composite VPT plot for {datestr}")
    
    # Find all MAN files for this date
    date_path = os.path.join(inpath, datestr)
    if not os.path.exists(date_path):
        print(f"Date directory does not exist: {date_path}")
        return
    
    man_files = sorted(glob.glob(os.path.join(date_path, '*_man_*.nc')))
    if not man_files:
        print(f"No MAN files found for {datestr}")
        return
    
    print(f"Found {len(man_files)} MAN files to search for vertical pointing dwells")
    
    # Collect all vertical pointing dwell segments
    dwell_data = []
    
    for ncfile in man_files:
        try:
            # Read radar data
            radar_ds = pyart.io.read_cfradial(ncfile)
            
            # Check if this file has phase_sequence metadata
            if 'phase_sequence' not in radar_ds.metadata:
                continue
            
            phase_list = radar_ds.metadata['phase_sequence'].split(', ')
            nsweeps = radar_ds.nsweeps
            
            # Get valid sweep indices (filter out all-transition sweeps if requested)
            valid_sweep_indices = get_valid_sweep_indices(radar_ds, skip_all_transition=skip_all_transition)
            
            # Find dwell sweeps that are vertical pointing
            for sweep_idx in range(nsweeps):
                if sweep_idx >= len(phase_list):
                    continue
                
                # Skip if this sweep is filtered out
                if sweep_idx not in valid_sweep_indices:
                    continue
                    
                phase = phase_list[sweep_idx]
                
                # Only include dwell or vertical_pointing phases
                if phase in ['dwell', 'vertical_pointing']:
                    # Extract this sweep to check elevation
                    start_idx = radar_ds.sweep_start_ray_index['data'][sweep_idx]
                    end_idx = radar_ds.sweep_end_ray_index['data'][sweep_idx]
                    
                    # Check mean elevation - only include if vertical pointing (>85°)
                    elevation_data = radar_ds.elevation['data'][start_idx:end_idx+1]
                    mean_elevation = np.mean(elevation_data)
                    
                    if mean_elevation <= 85.0:
                        print(f"  Skipping {phase} in {os.path.basename(ncfile)}, sweep {sweep_idx+1}: mean elevation {mean_elevation:.1f}° (not vertical pointing)")
                        continue
                    
                    # Extract this dwell sweep
                    dwell_sweep = radar_ds.extract_sweeps([sweep_idx])
                    
                    # Split dwell at time gaps to prevent pixel stretching
                    dwell_segments = split_radar_at_time_gaps(dwell_sweep, max_gap_seconds=600)
                    
                    # Add each segment separately
                    for segment in dwell_segments:
                        # Validate segment has sufficient dimensions for plotting
                        if segment.nrays < 2 or segment.ngates < 2:
                            print(f"    Skipping segment with insufficient dimensions: nrays={segment.nrays}, ngates={segment.ngates}")
                            continue
                        
                        # Validate field data shape matches coordinates
                        try:
                            field_shape = segment.fields['DBZ']['data'].shape
                            time_len = len(segment.time['data'])
                            range_len = len(segment.range['data'])
                            
                            if field_shape[0] != time_len or field_shape[1] != range_len:
                                print(f"    Skipping segment with mismatched field dimensions: field shape={field_shape}, time={time_len}, range={range_len}")
                                continue
                            
                            if time_len < 2 or range_len < 2:
                                print(f"    Skipping segment with insufficient coordinate length: time={time_len}, range={range_len}")
                                continue
                        except Exception as e:
                            print(f"    Skipping segment: could not validate field dimensions - {e}")
                            continue
                        
                        times = cftime.num2pydate(segment.time['data'], segment.time['units'])
                        
                        dwell_data.append({
                            'radar': segment,
                            'times': times,
                            'start_time': times[0]
                        })
                    
                    print(f"  Found vertical pointing {phase} in {os.path.basename(ncfile)}, sweep {sweep_idx+1} (elevation {mean_elevation:.1f}°), split into {len(dwell_segments)} segment(s)")
                    
        except Exception as e:
            print(f"  Error processing {ncfile}: {e}")
            continue
    
    if not dwell_data:
        print(f"No vertical pointing dwell segments found for {datestr}")
        return
    
    print(f"\nFound {len(dwell_data)} total vertical pointing dwell segments")
    
    # Sort by start time
    dwell_data.sort(key=lambda x: x['start_time'])
    
    # Get product version and instrument name from first file
    try:
        with nc4.Dataset(man_files[0]) as ds:
            product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
            instrument_name = ds.source if 'source' in ds.ncattrs() else 'ncas-mobile-ka-band-radar-1'
    except:
        product_version = VERSION
        instrument_name = 'ncas-mobile-ka-band-radar-1'
    
    # Set up output directory
    vpt_figpath = os.path.join(figpath, 'vpt')
    os.makedirs(vpt_figpath, exist_ok=True)
    
    # Set height limits
    hmax = 4 if blflag else 12
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    
    # Process each dwell segment and plot
    for idx, dwell_info in enumerate(dwell_data):
        radar = dwell_info['radar']
        times = dwell_info['times']
        
        # Re-validate dimensions (may have changed after processing)
        if radar.nrays < 2 or radar.ngates < 2:
            print(f"  Skipping segment {idx} with insufficient dimensions: nrays={radar.nrays}, ngates={radar.ngates}")
            continue
            
        # Validate that data dimensions match coordinate dimensions
        try:
            field_data = radar.fields['DBZ']['data']
            time_len = len(radar.time['data'])
            range_len = len(radar.range['data'])
            
            if field_data.shape[0] != time_len or field_data.shape[1] != range_len:
                print(f"  Skipping segment {idx} with mismatched field dimensions: field shape={field_data.shape}, time={time_len}, range={range_len}")
                continue
            
            if time_len < 2 or range_len < 2:
                print(f"  Skipping segment {idx} with insufficient coordinate length: time={time_len}, range={range_len}")
                continue
        except Exception as e:
            print(f"  Skipping segment {idx}: could not validate dimensions - {e}")
            continue
        
        try:
            display = pyart.graph.RadarDisplay(radar)
            
            # Extract field labels from metadata (once per dwell, they should be the same)
            if idx == 0:
                dbz_long = radar.fields['DBZ'].get('long_name', 'Reflectivity')
                dbz_units = radar.fields['DBZ'].get('units', 'dBZ')
                dbz_title = f"{dbz_long} ({dbz_units})"
                
                vel_long = radar.fields['VEL'].get('long_name', 'Velocity')
                vel_units = radar.fields['VEL'].get('units', 'm/s')
                vel_title = f"{vel_long} ({vel_units})"
                
                width_long = radar.fields['WIDTH'].get('long_name', 'Spectrum Width')
                width_units = radar.fields['WIDTH'].get('units', 'm/s')
                width_title = f"{width_long} ({width_units})"
                
                ldr_long = radar.fields['LDR'].get('long_name', 'LDR')
                ldr_units = radar.fields['LDR'].get('units', 'dB')
                ldr_title = f"{ldr_long} ({ldr_units})"
            
            # Plot each field (will overlay all dwells)
            if idx == 0:
                # First dwell - create the plots
                _plot_vpt_field_composite(display, axes[0,0], 'DBZ', dbz_title, 'Reflectivity (dBZ)',
                               COLORMAPS['dbz'], -60, 40, hmax, times, first=True)
                
                vel_field = radar.fields['VEL']
                vel_min = vel_field.get('field_limit_lower', -10)
                vel_max = vel_field.get('field_limit_upper', 10)
                _plot_vpt_field_composite(display, axes[0,1], 'VEL', vel_title, 'Velocity (m/s)',
                               COLORMAPS['vel'], vel_min, vel_max, hmax, times, first=True)
                
                _plot_vpt_field_composite(display, axes[1,0], 'WIDTH', width_title, 'Spectrum Width (m/s)',
                               COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True, first=True)
                
                _plot_vpt_field_composite(display, axes[1,1], 'LDR', ldr_title, 'LDR (dB)',
                               COLORMAPS['ldr'], -35, 5, hmax, times, first=True)
            else:
                # Subsequent dwells - add to existing plots
                _plot_vpt_field_composite(display, axes[0,0], 'DBZ', dbz_title, 'Reflectivity (dBZ)',
                               COLORMAPS['dbz'], -60, 40, hmax, times, first=False)
                
                vel_field = radar.fields['VEL']
                vel_min = vel_field.get('field_limit_lower', -10)
                vel_max = vel_field.get('field_limit_upper', 10)
                _plot_vpt_field_composite(display, axes[0,1], 'VEL', vel_title, 'Velocity (m/s)',
                               COLORMAPS['vel'], vel_min, vel_max, hmax, times, first=False)
                
                _plot_vpt_field_composite(display, axes[1,0], 'WIDTH', width_title, 'Spectrum Width (m/s)',
                               COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True, first=False)
                
                _plot_vpt_field_composite(display, axes[1,1], 'LDR', ldr_title, 'LDR (dB)',
                               COLORMAPS['ldr'], -35, 5, hmax, times, first=False)
        except Exception as e:
            print(f"  Error plotting segment {idx}: {e}")
            continue
    
    # Set title
    date_obj = datetime.datetime.strptime(datestr, '%Y%m%d')
    fig.suptitle(f'PICASSO MAN Vertical Pointing Dwell Composite VPT - {date_obj.strftime("%Y-%m-%d")}\n'
                 f'{instrument_name}\n'
                 f'{len(dwell_data)} vertical pointing dwell segments', 
                 fontsize=14, fontweight='bold')
    
    # Save
    figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{datestr}_man-dwell-vpt_l1_{product_version}.png'
    plt.tight_layout()
    plt.savefig(os.path.join(vpt_figpath, figname), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"MAN vertical pointing dwell composite VPT saved to: {vpt_figpath}")
    print(f"  Combined {len(dwell_data)} dwell segments spanning "
          f"{dwell_data[0]['start_time'].strftime('%H:%M')} to "
          f"{dwell_data[-1]['times'][-1].strftime('%H:%M')} UTC")

def _plot_vpt_field_composite(display, ax, field_name, title, colorbar_label, cmap, vmin, vmax, hmax, times, log_scale=False, first=True):
    """
    Helper function to plot VPT field for composite plots using manual pcolormesh.
    
    Args:
        display: PyART RadarDisplay object
        ax: Matplotlib axis
        field_name: Name of field to plot
        title: Title for the panel
        colorbar_label: Label for colorbar
        cmap: Colormap name
        vmin, vmax: Color limits
        hmax: Maximum height
        times: Array of datetime objects
        log_scale: Whether to use log scale
        first: If True, set up axes with colorbar; if False, only add data
    """
    radar = display._radar
    
    # Get data
    field_data = radar.fields[field_name]['data']
    sweep_range = radar.range['data'] / 1000.0  # Convert to km
    
    # Create meshgrid for plotting
    time_nums = mdates.date2num(times)
    X, Y = np.meshgrid(time_nums, sweep_range)
    
    # Transpose field data for pcolormesh (expects [range, time])
    field_data_plot = field_data.T
    
    # Plot data with appropriate normalization
    if log_scale:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        pcm = ax.pcolormesh(X, Y, field_data_plot, cmap=cmap, norm=norm, shading='nearest')
    else:
        pcm = ax.pcolormesh(X, Y, field_data_plot, vmin=vmin, vmax=vmax, 
                            cmap=cmap, shading='nearest')
    
    if first:
        # First segment - set up axes and colorbar
        plt.colorbar(pcm, ax=ax, label=colorbar_label, orientation='horizontal', pad=0.15, aspect=30)
        ax.set_ylim(0, hmax)
        ax.set_xlabel('Time (UTC)')
        ax.set_ylabel('Height (km)')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

def collect_ppi_bounds_for_day(ppi_files):
    """
    Collect and merge geospatial bounds from all PPI files for a day.
    
    Args:
        ppi_files: List of PPI NetCDF file paths
        
    Returns:
        tuple: (lat_min, lat_max, lon_min, lon_max) or None if no bounds found
    """
    all_bounds = []
    
    for ncfile in ppi_files:
        try:
            with nc4.Dataset(ncfile, 'r') as ds:
                if 'geospatial_bounds' in ds.ncattrs():
                    geospatial_bounds_str = ds.getncattr('geospatial_bounds')
                    bounds = parse_geospatial_bounds(geospatial_bounds_str)
                    if bounds is not None:
                        all_bounds.append(bounds)
        except Exception as e:
            print(f"  Warning: Could not read bounds from {os.path.basename(ncfile)}: {e}")
            continue
    
    if not all_bounds:
        return None
    
    # Calculate union of all bounds (min of mins, max of maxs)
    lat_mins = [b[0] for b in all_bounds]
    lat_maxs = [b[1] for b in all_bounds]
    lon_mins = [b[2] for b in all_bounds]
    lon_maxs = [b[3] for b in all_bounds]
    
    unified_bounds = (
        min(lat_mins),  # lat_min
        max(lat_maxs),  # lat_max
        min(lon_mins),  # lon_min
        max(lon_maxs)   # lon_max
    )
    
    return unified_bounds

def process_day(datestr, inpath, figpath, blflag=False, skip_all_transition=False, custom_bounds=None, ppi_maps_only=False, azimuth_offset=0.0, zoom_km=None, zoom_level=None):
    """
    Process all radar files for a given day.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Output directory path
        blflag: Boundary layer flag
        skip_all_transition: If True, skip sweeps where 100% of rays are antenna_transition=1
        custom_bounds: Optional tuple (lat_min, lat_max, lon_min, lon_max) for PPI maps
        ppi_maps_only: If True, only generate PPI map plots
        azimuth_offset: Additional azimuth offset to apply in degrees (default: 0.0)
    """
    print(f"\n{'='*60}")
    print(f"Processing PICASSO quicklooks for {datestr}")
    print(f"{'='*60}\n")
    
    # Find all files for this date
    date_path = os.path.join(inpath, datestr)
    if not os.path.exists(date_path):
        print(f"Date directory does not exist: {date_path}")
        return
    
    if not ppi_maps_only:
        # Process RHI files
        rhi_files = sorted(glob.glob(os.path.join(date_path, '*_rhi_*.nc')))
        print(f"\nFound {len(rhi_files)} RHI files")
        for ncfile in rhi_files:
            try:
                make_picasso_rhi_plot(ncfile, figpath, blflag, skip_all_transition)
            except Exception as e:
                print(f"Error processing RHI file {ncfile}: {e}")
        
        # Process MAN files
        man_files = sorted(glob.glob(os.path.join(date_path, '*_man_*.nc')))
        print(f"\nFound {len(man_files)} MAN files")
        for ncfile in man_files:
            try:
                make_picasso_man_plot(ncfile, figpath, blflag, skip_all_transition)
            except Exception as e:
                print(f"Error processing MAN file {ncfile}: {e}")
        
        # Process VPT files
        vpt_files = sorted(glob.glob(os.path.join(date_path, '*_vpt_*.nc')))
        print(f"\nFound {len(vpt_files)} VPT files")
        for ncfile in vpt_files:
            try:
                make_picasso_vpt_plot(ncfile, figpath, blflag, skip_all_transition)
            except Exception as e:
                print(f"Error processing VPT file {ncfile}: {e}")
        
        # Create composite VPT plot from MAN vertical pointing dwell segments
        # Disabled - not needed
        # try:
        #     make_picasso_man_dwell_composite(datestr, inpath, figpath, blflag, skip_all_transition)
        # except Exception as e:
        #     print(f"Error creating MAN vertical pointing dwell composite: {e}")
        
        # Create merged composite VPT plot (VPT files + MAN vertical pointing dwells)
        try:
            make_picasso_vpt_composite(datestr, inpath, figpath, blflag, skip_all_transition)
        except Exception as e:
            print(f"Error creating composite VPT plot: {e}")
        
        # Process PPI files (standard polar plots)
        ppi_files = sorted(glob.glob(os.path.join(date_path, '*_ppi_*.nc')))
        print(f"\nFound {len(ppi_files)} PPI files")
        for ncfile in ppi_files:
            try:
                make_picasso_ppi_plot(ncfile, figpath, blflag, skip_all_transition, azimuth_offset)
            except Exception as e:
                print(f"Error processing PPI file {ncfile}: {e}")
    else:
        # PPI maps only mode - just find PPI files
        ppi_files = sorted(glob.glob(os.path.join(date_path, '*_ppi_*.nc')))
        print(f"\nFound {len(ppi_files)} PPI files (maps only mode)")
    
    # Process PPI map plots with unified bounds for the day
    if ppi_files:
        # Use custom bounds if provided, otherwise collect from files
        if custom_bounds is not None:
            lat_min, lat_max, lon_min, lon_max = custom_bounds
            print(f"\nGenerating PPI map plots with custom bounds")
            print(f"  Custom bounds:")
            print(f"    Latitude:  [{lat_min:.3f}, {lat_max:.3f}]")
            print(f"    Longitude: [{lon_min:.3f}, {lon_max:.3f}]")
            
            # Generate map plots with custom bounds
            for ncfile in ppi_files:
                try:
                    make_picasso_ppi_map_plot(ncfile, figpath, blflag, skip_all_transition, 
                                             basemap_type="opentopo",
                                             lat_bounds=(lat_min, lat_max),
                                             lon_bounds=(lon_min, lon_max),
                                             azimuth_offset=azimuth_offset,
                                             zoom_level=zoom_level)
                except Exception as e:
                    print(f"  Error processing PPI map for {ncfile}: {e}")
        else:
            print(f"\nGenerating PPI map plots with unified bounds")
            
            # Collect bounds from all PPI files for the day
            print(f"  Collecting geospatial bounds from {len(ppi_files)} PPI files...")
            unified_bounds = collect_ppi_bounds_for_day(ppi_files)
            
            if unified_bounds is not None:
                lat_min, lat_max, lon_min, lon_max = unified_bounds
                print(f"  Unified bounds for {datestr}:")
                print(f"    Latitude:  [{lat_min:.3f}, {lat_max:.3f}]")
                print(f"    Longitude: [{lon_min:.3f}, {lon_max:.3f}]")
                
                # Generate map plots with unified bounds
                for ncfile in ppi_files:
                    try:
                        make_picasso_ppi_map_plot(ncfile, figpath, blflag, skip_all_transition, 
                                                 basemap_type="opentopo",
                                                 lat_bounds=(lat_min, lat_max),
                                                 lon_bounds=(lon_min, lon_max),
                                                 azimuth_offset=azimuth_offset,
                                                 zoom_level=zoom_level)
                    except Exception as e:
                        print(f"  Error processing PPI map for {ncfile}: {e}")
            else:
                print(f"  Warning: Could not determine unified bounds, using individual bounds")
                # Fall back to individual bounds
                for ncfile in ppi_files:
                    try:
                        make_picasso_ppi_map_plot(ncfile, figpath, blflag, skip_all_transition, 
                                                 basemap_type="opentopo",
                                                 azimuth_offset=azimuth_offset,
                                                 zoom_level=zoom_level)
                    except Exception as e:
                        print(f"  Error processing PPI map for {ncfile}: {e}")
    
    # Generate zoomed-in PPI map plots (use --zoom-km or default to +/- 2km from radar)
    if ppi_files:
        zoom_distance = zoom_km if zoom_km is not None else 2.0
        print(f"\nGenerating zoomed PPI map plots (+/- {zoom_distance}km)")
        for ncfile in ppi_files:
            try:
                make_picasso_ppi_map_plot(ncfile, figpath, blflag, skip_all_transition, 
                                         basemap_type="opentopo",
                                         zoom_km=zoom_distance,
                                         azimuth_offset=azimuth_offset,
                                         zoom_level=zoom_level)
            except Exception as e:
                print(f"  Error processing zoomed PPI map for {ncfile}: {e}")
    
    print(f"\n{'='*60}")
    print(f"Completed PICASSO quicklooks for {datestr}")
    print(f"{'='*60}\n")


def compute_aircraft_tracking_angles(kml_file, radar_lat, radar_lon, radar_alt, 
                                     output_plot=None):
    """
    Compute azimuth and elevation angles from radar to aircraft flight track.
    
    This function reads a KML file containing aircraft GPS coordinates and computes
    the azimuth and elevation angles that the radar would need to track the aircraft.
    This is useful for:
    - Validating MAN (manual tracking) scan performance
    - Planning future tracking scans
    - Comparing actual vs. ideal tracking angles
    
    Args:
        kml_file: Path to KML file with aircraft flight track
        radar_lat: Radar latitude (degrees North)
        radar_lon: Radar longitude (degrees East)
        radar_alt: Radar altitude (meters above sea level)
        output_plot: Optional path to save figure (if None, displays interactively)
        
    Returns:
        Dictionary containing:
            - 'times': Array of datetime objects
            - 'azimuth': Array of azimuth angles (degrees, 0=N, 90=E)
            - 'elevation': Array of elevation angles (degrees above horizon)
            - 'range': Array of slant ranges (km)
            - 'lat': Array of aircraft latitudes
            - 'lon': Array of aircraft longitudes
            - 'alt': Array of aircraft altitudes (m)
    
    Example:
        >>> angles = compute_aircraft_tracking_angles(
        ...     'flight_track.kml', 
        ...     radar_lat=51.144, 
        ...     radar_lon=-1.438, 
        ...     radar_alt=83.0
        ... )
        >>> print(f"Max elevation: {np.max(angles['elevation']):.1f}°")
    """
    import xml.etree.ElementTree as ET
    from math import radians, degrees, sin, cos, atan2, sqrt, asin
    
    print(f"Reading flight track from: {kml_file}")
    print(f"Radar location: lat={radar_lat:.4f}°, lon={radar_lon:.4f}°, alt={radar_alt:.1f}m")
    
    # Parse KML file
    tree = ET.parse(kml_file)
    root = tree.getroot()
    
    # Handle KML namespace
    ns = {'kml': 'http://www.opengis.net/kml/2.2'}
    
    # Extract coordinates and times
    track_data = []
    
    # Try to find gx:Track elements (Google Earth GPS track)
    for track in root.findall('.//gx:Track', {'gx': 'http://www.google.com/kml/ext/2.2', **ns}):
        when_elements = track.findall('.//kml:when', ns)
        coord_elements = track.findall('.//gx:coord', {'gx': 'http://www.google.com/kml/ext/2.2'})
        
        for when, coord in zip(when_elements, coord_elements):
            time_str = when.text.strip()
            coords = coord.text.strip().split()
            if len(coords) >= 2:
                lon, lat = float(coords[0]), float(coords[1])
                alt = float(coords[2]) if len(coords) > 2 else 0.0
                
                # Parse timestamp
                try:
                    timestamp = datetime.datetime.fromisoformat(time_str.replace('Z', '+00:00'))
                except:
                    timestamp = datetime.datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ')
                
                track_data.append({'time': timestamp, 'lat': lat, 'lon': lon, 'alt': alt})
    
    # If no gx:Track found, try standard Placemark/Point coordinates
    if not track_data:
        for placemark in root.findall('.//kml:Placemark', ns):
            time_elem = placemark.find('.//kml:TimeStamp/kml:when', ns)
            coord_elem = placemark.find('.//kml:Point/kml:coordinates', ns)
            
            if time_elem is not None and coord_elem is not None:
                time_str = time_elem.text.strip()
                coords = coord_elem.text.strip().split(',')
                
                if len(coords) >= 2:
                    lon, lat = float(coords[0]), float(coords[1])
                    alt = float(coords[2]) if len(coords) > 2 else 0.0
                    
                    try:
                        timestamp = datetime.datetime.fromisoformat(time_str.replace('Z', '+00:00'))
                    except:
                        timestamp = datetime.datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ')
                    
                    track_data.append({'time': timestamp, 'lat': lat, 'lon': lon, 'alt': alt})
    
    if not track_data:
        raise ValueError(f"No track data found in KML file: {kml_file}")
    
    print(f"Found {len(track_data)} track points")
    
    # Convert to arrays
    times = np.array([d['time'] for d in track_data])
    aircraft_lat = np.array([d['lat'] for d in track_data])
    aircraft_lon = np.array([d['lon'] for d in track_data])
    aircraft_alt = np.array([d['alt'] for d in track_data])
    
    # Calculate azimuth and elevation for each point
    n_points = len(track_data)
    azimuth = np.zeros(n_points)
    elevation = np.zeros(n_points)
    slant_range = np.zeros(n_points)
    
    for i in range(n_points):
        # Convert to radians
        lat1, lon1 = radians(radar_lat), radians(radar_lon)
        lat2, lon2 = radians(aircraft_lat[i]), radians(aircraft_lon[i])
        
        # Calculate azimuth (bearing)
        dlon = lon2 - lon1
        x = sin(dlon) * cos(lat2)
        y = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)
        azimuth[i] = (degrees(atan2(x, y)) + 360) % 360
        
        # Calculate horizontal distance using Haversine formula
        dlat = lat2 - lat1
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a))
        earth_radius_km = 6371.0
        horizontal_dist_km = earth_radius_km * c
        
        # Calculate elevation angle
        altitude_diff_m = aircraft_alt[i] - radar_alt
        altitude_diff_km = altitude_diff_m / 1000.0
        
        # Slant range
        slant_range[i] = sqrt(horizontal_dist_km**2 + altitude_diff_km**2)
        
        # Elevation angle (degrees above horizon)
        if horizontal_dist_km > 0:
            elevation[i] = degrees(atan2(altitude_diff_km, horizontal_dist_km))
        else:
            elevation[i] = 90.0 if altitude_diff_m > 0 else -90.0
    
    print(f"Azimuth range: {np.min(azimuth):.1f}° to {np.max(azimuth):.1f}°")
    print(f"Elevation range: {np.min(elevation):.1f}° to {np.max(elevation):.1f}°")
    print(f"Slant range: {np.min(slant_range):.2f} to {np.max(slant_range):.2f} km")
    
    # Create plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
    
    # Plot azimuth vs time
    ax1.plot(times, azimuth, 'b-', linewidth=1.5)
    ax1.set_ylabel('Azimuth (°)', fontsize=12)
    ax1.set_title(f'Aircraft Tracking Angles from Radar\nRadar: {radar_lat:.4f}°N, {radar_lon:.4f}°E, {radar_alt:.0f}m', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Plot elevation vs time
    ax2.plot(times, elevation, 'r-', linewidth=1.5)
    ax2.axhline(88, color='orange', linestyle='--', alpha=0.5, label='88° (dwell threshold)')
    ax2.axhline(89.5, color='red', linestyle='--', alpha=0.5, label='89.5° (vertical threshold)')
    ax2.set_ylabel('Elevation (°)', fontsize=12)
    ax2.legend(loc='best', fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Plot slant range vs time
    ax3.plot(times, slant_range, 'g-', linewidth=1.5)
    ax3.set_ylabel('Slant Range (km)', fontsize=12)
    ax3.set_xlabel('Time (UTC)', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Rotate x-axis labels
    for ax in [ax1, ax2, ax3]:
        ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    
    if output_plot:
        plt.savefig(output_plot, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_plot}")
        plt.close()
    else:
        plt.show()
    
    # Return results
    return {
        'times': times,
        'azimuth': azimuth,
        'elevation': elevation,
        'range': slant_range,
        'lat': aircraft_lat,
        'lon': aircraft_lon,
        'alt': aircraft_alt
    }


def main():
    """Main entry point."""
    # Parse command line
    datestr, inpath, outpath, blflag, skip_all_transition, custom_bounds, ppi_maps_only, azimuth_offset, zoom_km, zoom_level = parse_command_line()
    
    # Setup paths
    inpath, figpath = setup_paths(datestr, inpath, outpath)
    
    print(f"PICASSO Quicklook Generator v{VERSION}")
    print(f"Date: {datestr}")
    print(f"Input path: {inpath}")
    print(f"Output path: {figpath}")
    print(f"Boundary layer mode: {blflag}")
    print(f"Skip all-transition sweeps: {'enabled' if skip_all_transition else 'disabled'}")
    if azimuth_offset != 0.0:
        print(f"Azimuth offset: {azimuth_offset}°")
    if zoom_km is not None:
        print(f"Map zoom distance: +/- {zoom_km} km")
    if zoom_level is not None:
        print(f"Basemap tile zoom level: {zoom_level}")
    if ppi_maps_only:
        print(f"Mode: PPI maps only (ppi_map + ppi_map_zoom)")
    
    # Process the day
    process_day(datestr, inpath, figpath, blflag, skip_all_transition, custom_bounds, ppi_maps_only, azimuth_offset, zoom_km, zoom_level)

if __name__ == "__main__":
    main()
