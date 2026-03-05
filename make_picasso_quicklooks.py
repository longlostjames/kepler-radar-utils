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
    python make_picasso_quicklooks.py -d YYYYMMDD [-i input_path] [-o output_path] [-b]

Arguments:
    -d, --date:     Date string in YYYYMMDD format
    -i, --inpath:   Input directory containing CF-Radial files (optional)
    -o, --outpath:  Output directory for quicklook images (optional)
    -b:             Boundary layer flag (limits plots to 4km height)

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

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
import cmocean

# Configuration
VERSION = '1.0.0'
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
        tuple: (datestr, inpath, outpath, blflag)
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:i:o:b", 
                                   ["date=", "inpath=", "outpath="])
    except getopt.GetoptError as err:
        print(f"Error: {err}")
        print("Usage: python make_picasso_quicklooks.py -d YYYYMMDD [-i input_path] [-o output_path] [-b]")
        sys.exit(2)

    # Default values
    data_date = datetime.datetime.now()
    datestr = data_date.strftime('%Y%m%d')
    inpath = None
    outpath = None
    blflag = False

    for option, argument in opts:
        if option in ("-d", "--date"):
            datestr = argument
        elif option in ("-i", "--inpath"):
            inpath = argument
        elif option in ("-o", "--outpath"):
            outpath = argument
        elif option == "-b":
            blflag = True
        else:
            print(f"Unhandled option: {option}")
            sys.exit(2)

    return datestr, inpath, outpath, blflag

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
        
        segments.append(segment)
    
    print(f"  Split radar data into {len(segments)} segments at time gaps > {max_gap_seconds}s")
    return segments

def make_picasso_rhi_plot(ncfile, figpath, blflag=False):
    """
    Create RHI (Range Height Indicator) quicklook plot for PICASSO data.
    
    Args:
        ncfile: Path to CF-Radial NetCDF file
        figpath: Output directory for figure
        blflag: Boundary layer flag for plot limits
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
    
    for sweep_idx in range(nsweeps):
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
    
    # Reflectivity (DBZ)
    display.plot_rhi("DBZ", ax=axes[0,0], sweep=sweep_idx, vmin=-60, vmax=40,
                     cmap=COLORMAPS['dbz'], colorbar_flag=False,
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[0,0], hmax, xmin, xmax)
    axes[0,0].set_title('Reflectivity')
    plt.colorbar(axes[0,0].collections[-1], ax=axes[0,0], orientation='horizontal', 
                 label='Reflectivity (dBZ)', pad=0.25, aspect=30)
    
    # Velocity (VEL)  
    display.plot_rhi("VEL", ax=axes[1,0], sweep=sweep_idx, vmin=vel_min, vmax=vel_max,
                     cmap=COLORMAPS['vel'], colorbar_flag=False,
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[1,0], hmax, xmin, xmax)
    axes[1,0].set_title('Doppler Velocity')
    plt.colorbar(axes[1,0].collections[-1], ax=axes[1,0], orientation='horizontal',
                 label='Velocity (m/s)', pad=0.25, aspect=30)
    
    # Spectrum width (WIDTH)
    display.plot_rhi("WIDTH", ax=axes[1,1], sweep=sweep_idx,
                     norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                     cmap=COLORMAPS['width'], colorbar_flag=False,
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[1,1], hmax, xmin, xmax)
    axes[1,1].set_title('Spectrum Width')
    plt.colorbar(axes[1,1].collections[-1], ax=axes[1,1], orientation='horizontal',
                 label='Spectrum Width (m/s)', pad=0.25, aspect=30)
    
    # Linear depolarization ratio (LDR)
    display.plot_rhi("LDR", ax=axes[0,1], sweep=sweep_idx, vmin=-35, vmax=5,
                     cmap=COLORMAPS['ldr'], colorbar_flag=False,
                     gatefilter=gatefilter, filter_transitions=True)
    _setup_rhi_axes(axes[0,1], hmax, xmin, xmax)
    axes[0,1].set_title('Linear Depolarization Ratio')
    plt.colorbar(axes[0,1].collections[-1], ax=axes[0,1], orientation='horizontal',
                 label='LDR (dB)', pad=0.25, aspect=30)

def _setup_rhi_axes(ax, hmax, xmin, xmax):
    """Helper function to set up RHI plot axes."""
    ax.set_ylim(0, hmax)
    ax.set_xlim(xmin, xmax)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

def _plot_ppi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max, max_range):
    """Helper function to plot PPI fields on axes for turning scans."""
    
    # Reflectivity (DBZ)
    display.plot_ppi("DBZ", ax=axes[0,0], sweep=sweep_idx, vmin=-60, vmax=40,
                     cmap=COLORMAPS['dbz'], colorbar_flag=False,
                     gatefilter=gatefilter)
    _setup_ppi_axes(axes[0,0], max_range)
    axes[0,0].set_title('Reflectivity')
    plt.colorbar(axes[0,0].collections[-1], ax=axes[0,0], orientation='vertical',
                 label='Reflectivity (dBZ)', pad=0.02, aspect=20)
    
    # Velocity (VEL)  
    display.plot_ppi("VEL", ax=axes[1,0], sweep=sweep_idx, vmin=vel_min, vmax=vel_max,
                     cmap=COLORMAPS['vel'], colorbar_flag=False,
                     gatefilter=gatefilter)
    _setup_ppi_axes(axes[1,0], max_range)
    axes[1,0].set_title('Doppler Velocity')
    plt.colorbar(axes[1,0].collections[-1], ax=axes[1,0], orientation='vertical',
                 label='Velocity (m/s)', pad=0.02, aspect=20)
    
    # Spectrum width (WIDTH)
    display.plot_ppi("WIDTH", ax=axes[1,1], sweep=sweep_idx,
                     norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                     cmap=COLORMAPS['width'], colorbar_flag=False,
                     gatefilter=gatefilter)
    _setup_ppi_axes(axes[1,1], max_range)
    axes[1,1].set_title('Spectrum Width')
    plt.colorbar(axes[1,1].collections[-1], ax=axes[1,1], orientation='vertical',
                 label='Spectrum Width (m/s)', pad=0.02, aspect=20)
    
    # Linear depolarization ratio (LDR)
    display.plot_ppi("LDR", ax=axes[0,1], sweep=sweep_idx, vmin=-35, vmax=5,
                     cmap=COLORMAPS['ldr'], colorbar_flag=False,
                     gatefilter=gatefilter)
    _setup_ppi_axes(axes[0,1], max_range)
    axes[0,1].set_title('Linear Depolarization Ratio')
    plt.colorbar(axes[0,1].collections[-1], ax=axes[0,1], orientation='vertical',
                 label='LDR (dB)', pad=0.02, aspect=20)

def _setup_ppi_axes(ax, max_range):
    """Helper function to set up PPI plot axes."""
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

def _plot_sweep_angles(axes, sweep_times, elevation_data, azimuth_data, mode, 
                       radar_ds, start_ray_idx, end_ray_idx):
    """Plot elevation and azimuth angles vs time for a sweep.
    
    Points marked as antenna_transition are plotted in grey.
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
    
    # Set elevation limits based on phase type
    if mode in ['dwell', 'vertical_pointing']:
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

def make_picasso_man_plot(ncfile, figpath, blflag=False):
    """
    Create MAN (Manual tracking) quicklook plots for PICASSO aircraft tracking data.
    
    This creates special visualizations for MAN scans including:
    - Overview plot with angle timeseries and full-track RHI
    - Individual plots for each phase/sweep (upward, dwell, downward, turning)
    
    Sweeps with >50% of rays flagged as antenna_transition (using CF-Radial
    standard coordinate variable) are skipped to avoid plotting rapid repositioning movements.
    
    Args:
        ncfile: Path to CF-Radial NetCDF file
        figpath: Output directory for figure
        blflag: Boundary layer flag for plot limits
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
    
    # Get product version
    try:
        with nc4.Dataset(ncfile) as ds:
            product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
    except:
        product_version = VERSION
    
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
    
    # Export antenna transition indices to CSV
    _export_antenna_transitions_to_csv(radar_ds, ncfile, man_figpath)
    
    dtime0 = cftime.num2pydate(radar_ds.time['data'][0], radar_ds.time['units'])
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S")
    
    # Create overview plot with angle timeseries
    _make_man_overview_plot(radar_ds, gatefilter, man_figpath, dtime0_str, product_version, blflag)
    
    # Create individual sweep/phase plots
    _make_man_sweep_plots(radar_ds, gatefilter, man_figpath, dtime0_str, product_version, blflag)
    
    print(f"MAN plots saved to: {man_figpath}")

def _make_man_overview_plot(radar_ds, gatefilter, figpath, dtime_str, version, blflag):
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
    sweep_modes = []
    sweep_starts = []
    sweep_ends = []
    
    # Check if phase_sequence metadata exists (from phase-split MAN scans)
    phase_sequence = None
    if 'phase_sequence' in radar_ds.metadata:
        phase_list = radar_ds.metadata['phase_sequence'].split(', ')
        if len(phase_list) == nsweeps:
            phase_sequence = phase_list
    
    for i in range(nsweeps):
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
    fig.suptitle(f'PICASSO MAN Scan - Aircraft Tracking\n{times[0].strftime("%Y-%m-%d %H:%M:%S")} UTC', 
                 fontsize=14, fontweight='bold')
    
    # Save
    figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime_str}_man_overview_l1_{version}.png'
    plt.savefig(os.path.join(figpath, figname), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Overview plot created: {figname}")

def _make_man_sweep_plots(radar_ds, gatefilter, figpath, dtime_str, version, blflag):
    """Create individual plots for each MAN scan phase/sweep."""
    
    display = pyart.graph.RadarDisplay(radar_ds)
    hmax, xmin, xmax = setup_plot_limits(blflag, 'man')
    
    # For PPI plots, use max_range based on boundary layer flag
    max_range = 20 if blflag else 30
    
    vel_field = radar_ds.fields['VEL']
    vel_min = vel_field.get('field_limit_lower', -10)
    vel_max = vel_field.get('field_limit_upper', 10)
    
    nsweeps = radar_ds.nsweeps
    
    # Check if phase_sequence metadata exists (from phase-split MAN scans)
    phase_sequence = None
    if 'phase_sequence' in radar_ds.metadata:
        phase_list = radar_ds.metadata['phase_sequence'].split(', ')
        if len(phase_list) == nsweeps:
            phase_sequence = phase_list
    
    for sweep_idx in range(nsweeps):
        # Get sweep phase info - use phase_sequence if available
        if phase_sequence:
            mode = phase_sequence[sweep_idx]
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
        
        # Check if this sweep consists mostly of antenna transitions
        # using CF-Radial standard antenna_transition coordinate variable
        start_ray_idx = radar_ds.sweep_start_ray_index['data'][sweep_idx]
        end_ray_idx = radar_ds.sweep_end_ray_index['data'][sweep_idx]
        
        skip_sweep = False
        if hasattr(radar_ds, 'antenna_transition') and radar_ds.antenna_transition is not None:
            if 'data' in radar_ds.antenna_transition:
                transition_flags = radar_ds.antenna_transition['data'][start_ray_idx:end_ray_idx+1]
                transition_fraction = np.sum(transition_flags) / len(transition_flags)
                # Skip if more than 50% of rays are in transition
                if transition_fraction > 0.5:
                    print(f"  Skipping sweep {sweep_idx+1}/{nsweeps} (antenna transition: {transition_fraction:.1%} of rays)")
                    skip_sweep = True
        
        if skip_sweep:
            continue
        
        # Get sweep start and end times
        start_ray_idx = radar_ds.sweep_start_ray_index['data'][sweep_idx]
        end_ray_idx = radar_ds.sweep_end_ray_index['data'][sweep_idx]
        
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
        
        if mode == 'turning':
            _plot_ppi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max, max_range)
            subtitle = f'Mean Elevation: {mean_elevation:.1f}° | {time_range_str}'
        elif mode in ['dwell', 'vertical_pointing']:
            # Use VPT-style plotting for near-vertical dwells (time-height)
            # These are essentially vertical pointing and show nothing useful in RHI format
            print(f"    Using VPT format for {mode} phase")
            _plot_vpt_fields_single_sweep(display, axes, sweep_idx, gatefilter, vel_min, vel_max, hmax)
            subtitle = f'Mean Elevation: {mean_elevation:.1f}° | {time_range_str}'
        else:
            _plot_rhi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max, hmax, xmin, xmax)
            subtitle = f'Mean Azimuth: {mean_azimuth:.1f}° | {time_range_str}'
        
        # Add elevation and azimuth panels
        _plot_sweep_angles(axes, sweep_times, elevation_data, azimuth_data, mode,
                          radar_ds, start_ray_idx, end_ray_idx)
        
        # Add phase information to title (use original case for display)
        mode_display = mode.replace('_', ' ').title().replace(' ', '_')
        phase_color = PHASE_COLORS.get(mode, PHASE_COLORS.get('other', '#808080'))
        fig.suptitle(f'MAN Scan - Phase: {mode_display.upper()} (Sweep {sweep_idx+1}/{nsweeps})\n{subtitle}', 
                    fontsize=12, fontweight='bold',
                    color=phase_color)
        
        # Save
        figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime_str}_man_sweep{sweep_idx:02d}_{mode}_l1_{version}.png'
        plt.tight_layout()
        plt.savefig(os.path.join(figpath, figname), dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Phase plot created: sweep {sweep_idx+1}/{nsweeps} ({mode})")


def _plot_vpt_fields_single_sweep(display, axes, sweep_idx, gatefilter, vel_min, vel_max, hmax):
    """Helper function to plot a single sweep in VPT (time-height) format.
    
    Used for dwell and vertical_pointing phases where RHI format shows nothing useful.
    Manually extracts and plots data to avoid PyART dimension issues with single sweeps.
    """
    radar = display._radar
    
    # Get sweep ray indices
    start_ray = radar.sweep_start_ray_index['data'][sweep_idx]
    end_ray = radar.sweep_end_ray_index['data'][sweep_idx]
    
    # Extract time and height data for this sweep
    sweep_times = cftime.num2pydate(radar.time['data'][start_ray:end_ray+1], radar.time['units'])
    sweep_range = radar.range['data'] / 1000.0  # Convert to km
    
    # Create meshgrid for plotting
    time_nums = mdates.date2num(sweep_times)
    X, Y = np.meshgrid(time_nums, sweep_range)
    
    # Plot Reflectivity (DBZ)
    field_data = radar.fields['DBZ']['data'][start_ray:end_ray+1, :].T
    if gatefilter is not None:
        field_mask = gatefilter.gate_excluded[start_ray:end_ray+1, :].T
        field_data = np.ma.masked_where(field_mask, field_data)
    
    pcm = axes[0,0].pcolormesh(X, Y, field_data, vmin=-60, vmax=40, 
                               cmap=COLORMAPS['dbz'], shading='nearest')
    plt.colorbar(pcm, ax=axes[0,0], label='Reflectivity (dBZ)', 
                 orientation='horizontal', pad=0.18, aspect=30)
    axes[0,0].set_ylim(0, hmax)
    axes[0,0].set_title('Reflectivity')
    axes[0,0].set_xlabel('Time (UTC)')
    axes[0,0].set_ylabel('Height (km)')
    axes[0,0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Plot Velocity (VEL)
    field_data = radar.fields['VEL']['data'][start_ray:end_ray+1, :].T
    if gatefilter is not None:
        field_mask = gatefilter.gate_excluded[start_ray:end_ray+1, :].T
        field_data = np.ma.masked_where(field_mask, field_data)
    
    pcm = axes[1,0].pcolormesh(X, Y, field_data, vmin=vel_min, vmax=vel_max,
                               cmap=COLORMAPS['vel'], shading='nearest')
    plt.colorbar(pcm, ax=axes[1,0], label='Velocity (m/s)',
                 orientation='horizontal', pad=0.18, aspect=30)
    axes[1,0].set_ylim(0, hmax)
    axes[1,0].set_title('Doppler Velocity')
    axes[1,0].set_xlabel('Time (UTC)')
    axes[1,0].set_ylabel('Height (km)')
    axes[1,0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Plot Spectrum width (WIDTH)
    field_data = radar.fields['WIDTH']['data'][start_ray:end_ray+1, :].T
    if gatefilter is not None:
        field_mask = gatefilter.gate_excluded[start_ray:end_ray+1, :].T
        field_data = np.ma.masked_where(field_mask, field_data)
    
    pcm = axes[1,1].pcolormesh(X, Y, field_data, 
                               norm=colors.LogNorm(vmin=0.1*np.sqrt(0.1), vmax=np.sqrt(10)),
                               cmap=COLORMAPS['width'], shading='nearest')
    plt.colorbar(pcm, ax=axes[1,1], label='Spectrum Width (m/s)',
                 orientation='horizontal', pad=0.18, aspect=30)
    axes[1,1].set_ylim(0, hmax)
    axes[1,1].set_title('Spectrum Width')
    axes[1,1].set_xlabel('Time (UTC)')
    axes[1,1].set_ylabel('Height (km)')
    axes[1,1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    # Plot Linear depolarization ratio (LDR)
    field_data = radar.fields['LDR']['data'][start_ray:end_ray+1, :].T
    if gatefilter is not None:
        field_mask = gatefilter.gate_excluded[start_ray:end_ray+1, :].T
        field_data = np.ma.masked_where(field_mask, field_data)
    
    pcm = axes[0,1].pcolormesh(X, Y, field_data, vmin=-35, vmax=5,
                               cmap=COLORMAPS['ldr'], shading='nearest')
    plt.colorbar(pcm, ax=axes[0,1], label='LDR (dB)',
                 orientation='horizontal', pad=0.18, aspect=30)
    axes[0,1].set_ylim(0, hmax)
    axes[0,1].set_title('Linear Depolarization Ratio')
    axes[0,1].set_xlabel('Time (UTC)')
    axes[0,1].set_ylabel('Height (km)')
    axes[0,1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))


def make_picasso_vpt_plot(ncfile, figpath, blflag=False):
    """
    Create VPT (Vertical Pointing) quicklook plot for PICASSO data.
    
    Args:
        ncfile: Path to CF-Radial NetCDF file
        figpath: Output directory for figure
        blflag: Boundary layer flag for plot limits
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
    
    # Get product version
    try:
        with nc4.Dataset(ncfile) as ds:
            product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
    except:
        product_version = VERSION
    
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
        
        # First segment - create plots with colorbars
        if seg_idx == 0:
            _plot_vpt_field(display, axes[0,0], 'DBZ', 'Reflectivity (dBZ)', 
                           COLORMAPS['dbz'], -60, 40, hmax, times)
            _plot_vpt_field(display, axes[0,1], 'VEL', 'Velocity (m/s)', 
                           COLORMAPS['vel'], vel_min, vel_max, hmax, times)
            _plot_vpt_field(display, axes[1,0], 'WIDTH', 'Spectrum Width (m/s)', 
                           COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True)
            _plot_vpt_field(display, axes[1,1], 'LDR', 'LDR (dB)', 
                           COLORMAPS['ldr'], -35, 5, hmax, times)
        else:
            # Subsequent segments - overlay without colorbars
            _plot_vpt_field_overlay(display, axes[0,0], 'DBZ', COLORMAPS['dbz'], -60, 40, hmax, times)
            _plot_vpt_field_overlay(display, axes[0,1], 'VEL', COLORMAPS['vel'], vel_min, vel_max, hmax, times)
            _plot_vpt_field_overlay(display, axes[1,0], 'WIDTH', COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True)
            _plot_vpt_field_overlay(display, axes[1,1], 'LDR', COLORMAPS['ldr'], -35, 5, hmax, times)
    
    fig.suptitle(f'PICASSO VPT - {dtime0.strftime("%Y-%m-%d")}', fontsize=14, fontweight='bold')
    
    # Save
    figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{dtime0_str}_vpt_l1_{product_version}.png'
    plt.tight_layout()
    plt.savefig(os.path.join(vpt_figpath, figname), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"VPT plot saved to: {vpt_figpath}")

def _plot_vpt_field(display, ax, field_name, label, cmap, vmin, vmax, hmax, times, log_scale=False):
    """Helper function to plot VPT field."""
    if log_scale:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    
    display.plot_vpt(field_name, ax=ax, time_axis_flag=True, edges=True,
                     vmin=vmin, vmax=vmax, cmap=cmap, 
                     colorbar_label=label, norm=norm, shading='auto')
    
    ax.set_ylim(0, hmax)
    ax.set_xlabel('Time (UTC)')
    ax.set_ylabel('Height (km)')
    ax.set_title(label)
    ax.grid(True, alpha=0.3)

def _plot_vpt_field_overlay(display, ax, field_name, cmap, vmin, vmax, hmax, times, log_scale=False):
    """Helper function to plot VPT field overlay (no colorbar)."""
    if log_scale:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    
    display.plot_vpt(field_name, ax=ax, time_axis_flag=True, edges=True,
                     vmin=vmin, vmax=vmax, cmap=cmap, 
                     colorbar_flag=False, norm=norm, shading='auto')

def make_picasso_vpt_composite(datestr, inpath, figpath, blflag=False):
    """
    Create composite VPT plot combining VPT sweeps and MAN dwell segments for a day.
    
    This function merges all vertical pointing data from:
    1. Dedicated VPT scan files
    2. Dwell segments from MAN (manual tracking) scans
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Output directory path
        blflag: Boundary layer flag for plot limits
    """
    print(f"\nCreating composite VPT plot (VPT + MAN dwells) for {datestr}")
    
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
            else:
                # Multiple segments
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
            print(f"  Indexed VPT file: {os.path.basename(ncfile)}")
            
        except Exception as e:
            print(f"  Error indexing VPT file {ncfile}: {e}")
            continue
    
    # 2. Get MAN dwell segment metadata
    man_files = sorted(glob.glob(os.path.join(date_path, '*_man_*.nc')))
    print(f"Found {len(man_files)} MAN files to search for dwells")
    
    for ncfile in man_files:
        try:
            radar_ds = pyart.io.read_cfradial(ncfile)
            
            # Check if this file has phase_sequence metadata
            if 'phase_sequence' not in radar_ds.metadata:
                del radar_ds
                continue
            
            phase_list = radar_ds.metadata['phase_sequence'].split(', ')
            nsweeps = radar_ds.nsweeps
            
            # Find dwell sweeps
            for sweep_idx in range(nsweeps):
                if sweep_idx >= len(phase_list):
                    continue
                    
                phase = phase_list[sweep_idx]
                
                if phase == 'dwell' or phase == 'vertical_pointing':
                    start_ray_idx = radar_ds.sweep_start_ray_index['data'][sweep_idx]
                    end_ray_idx = radar_ds.sweep_end_ray_index['data'][sweep_idx]
                    times = cftime.num2pydate(radar_ds.time['data'][start_ray_idx:end_ray_idx+1], 
                                             radar_ds.time['units'])
                    
                    vpt_segment_info.append({
                        'file': ncfile,
                        'type': 'man_dwell',
                        'sweep_idx': sweep_idx,
                        'start_time': times[0],
                        'source': f'MAN_{phase}'
                    })
                    
                    print(f"  Found {phase} in {os.path.basename(ncfile)}, sweep {sweep_idx+1}")
            
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
    
    # Get product version
    try:
        first_file = vpt_files[0] if vpt_files else man_files[0]
        with nc4.Dataset(first_file) as ds:
            product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
    except:
        product_version = VERSION
    
    # Set up output directory
    vpt_figpath = os.path.join(figpath, 'vpt')
    os.makedirs(vpt_figpath, exist_ok=True)
    
    # Set height limits
    hmax = 4 if blflag else 12
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    
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
                del radar_ds
                
            elif seg_info['type'] == 'man_dwell':
                radar_ds = pyart.io.read_cfradial(seg_info['file'])
                dwell_sweep = radar_ds.extract_sweeps([seg_info['sweep_idx']])
                radar_segment = split_radar_at_time_gaps(dwell_sweep, max_gap_seconds=600)[0]
                del radar_ds, dwell_sweep
            
            times = cftime.num2pydate(radar_segment.time['data'], radar_segment.time['units'])
            display = pyart.graph.RadarDisplay(radar_segment)
            
            if idx == 0:
                # First segment - create plots with colorbars
                _plot_vpt_field(display, axes[0,0], 'DBZ', 'Reflectivity (dBZ)', 
                               COLORMAPS['dbz'], -60, 40, hmax, times)
                
                vel_field = radar_segment.fields['VEL']
                vel_min = vel_field.get('field_limit_lower', -10)
                vel_max = vel_field.get('field_limit_upper', 10)
                _plot_vpt_field(display, axes[0,1], 'VEL', 'Velocity (m/s)', 
                               COLORMAPS['vel'], vel_min, vel_max, hmax, times)
                
                _plot_vpt_field(display, axes[1,0], 'WIDTH', 'Spectrum Width (m/s)', 
                               COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True)
                
                _plot_vpt_field(display, axes[1,1], 'LDR', 'LDR (dB)', 
                               COLORMAPS['ldr'], -35, 5, hmax, times)
            else:
                # Subsequent segments - overlay without colorbars
                _plot_vpt_field_overlay(display, axes[0,0], 'DBZ', COLORMAPS['dbz'], -60, 40, hmax, times)
                
                vel_field = radar_segment.fields['VEL']
                vel_min = vel_field.get('field_limit_lower', -10)
                vel_max = vel_field.get('field_limit_upper', 10)
                _plot_vpt_field_overlay(display, axes[0,1], 'VEL', COLORMAPS['vel'], vel_min, vel_max, hmax, times)
                
                _plot_vpt_field_overlay(display, axes[1,0], 'WIDTH', COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True)
                
                _plot_vpt_field_overlay(display, axes[1,1], 'LDR', COLORMAPS['ldr'], -35, 5, hmax, times)
            
            # Free memory immediately after plotting
            del radar_segment, display, times
            
        except Exception as e:
            print(f"  Error plotting segment {idx}: {e}")
            continue
    
    # Count sources
    vpt_count = sum(1 for s in vpt_segment_info if s['source'] == 'VPT')
    man_count = sum(1 for s in vpt_segment_info if s['source'] == 'MAN_dwell')
    
    # Set title
    date_obj = datetime.datetime.strptime(datestr, '%Y%m%d')
    fig.suptitle(f'PICASSO Composite VPT - {date_obj.strftime("%Y-%m-%d")}\n'
                 f'{len(vpt_segment_info)} segments (VPT: {vpt_count}, MAN dwells: {man_count})', 
                 fontsize=14, fontweight='bold')
    
    # Save
    figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{datestr}_vpt_composite_l1_{product_version}.png'
    plt.tight_layout()
    plt.savefig(os.path.join(vpt_figpath, figname), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Composite VPT plot saved to: {vpt_figpath}")
    print(f"  Combined {len(vpt_segment_info)} segments: {vpt_count} VPT, {man_count} MAN dwells")

def make_picasso_man_dwell_composite(datestr, inpath, figpath, blflag=False):
    """
    Create composite VPT-style plot from all dwell segments in MAN scans for a day.
    
    This function:
    1. Finds all MAN files for the specified date
    2. Extracts dwell (vertical pointing) segments from each file
    3. Combines them chronologically into a single time-height plot
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Output directory path
        blflag: Boundary layer flag for plot limits
    """
    print(f"\nCreating MAN dwell composite VPT plot for {datestr}")
    
    # Find all MAN files for this date
    date_path = os.path.join(inpath, datestr)
    if not os.path.exists(date_path):
        print(f"Date directory does not exist: {date_path}")
        return
    
    man_files = sorted(glob.glob(os.path.join(date_path, '*_man_*.nc')))
    if not man_files:
        print(f"No MAN files found for {datestr}")
        return
    
    print(f"Found {len(man_files)} MAN files to search for dwells")
    
    # Collect all dwell segments
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
            
            # Find dwell sweeps
            for sweep_idx in range(nsweeps):
                if sweep_idx >= len(phase_list):
                    continue
                    
                phase = phase_list[sweep_idx]
                
                if phase == 'dwell':
                    # Extract this dwell sweep
                    start_idx = radar_ds.sweep_start_ray_index['data'][sweep_idx]
                    end_idx = radar_ds.sweep_end_ray_index['data'][sweep_idx]
                    
                    # Extract subset of radar object for this sweep
                    dwell_sweep = radar_ds.extract_sweeps([sweep_idx])
                    
                    # Split dwell at time gaps to prevent pixel stretching
                    dwell_segments = split_radar_at_time_gaps(dwell_sweep, max_gap_seconds=600)
                    
                    # Add each segment separately
                    for segment in dwell_segments:
                        times = cftime.num2pydate(segment.time['data'], segment.time['units'])
                        
                        dwell_data.append({
                            'radar': segment,
                            'times': times,
                            'start_time': times[0]
                        })
                    
                    print(f"  Found dwell in {os.path.basename(ncfile)}, sweep {sweep_idx+1}, split into {len(dwell_segments)} segment(s)")
                    
        except Exception as e:
            print(f"  Error processing {ncfile}: {e}")
            continue
    
    if not dwell_data:
        print(f"No dwell segments found for {datestr}")
        return
    
    print(f"\nFound {len(dwell_data)} total dwell segments")
    
    # Sort by start time
    dwell_data.sort(key=lambda x: x['start_time'])
    
    # Get product version from first file
    try:
        with nc4.Dataset(man_files[0]) as ds:
            product_version = ds.product_version if 'product_version' in ds.ncattrs() else VERSION
    except:
        product_version = VERSION
    
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
        
        display = pyart.graph.RadarDisplay(radar)
        
        # Plot each field (will overlay all dwells)
        if idx == 0:
            # First dwell - create the plots
            _plot_vpt_field_composite(display, axes[0,0], 'DBZ', 'Reflectivity (dBZ)', 
                           COLORMAPS['dbz'], -60, 40, hmax, times, first=True)
            
            vel_field = radar.fields['VEL']
            vel_min = vel_field.get('field_limit_lower', -10)
            vel_max = vel_field.get('field_limit_upper', 10)
            _plot_vpt_field_composite(display, axes[0,1], 'VEL', 'Velocity (m/s)', 
                           COLORMAPS['vel'], vel_min, vel_max, hmax, times, first=True)
            
            _plot_vpt_field_composite(display, axes[1,0], 'WIDTH', 'Spectrum Width (m/s)', 
                           COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True, first=True)
            
            _plot_vpt_field_composite(display, axes[1,1], 'LDR', 'LDR (dB)', 
                           COLORMAPS['ldr'], -35, 5, hmax, times, first=True)
        else:
            # Subsequent dwells - add to existing plots
            _plot_vpt_field_composite(display, axes[0,0], 'DBZ', 'Reflectivity (dBZ)', 
                           COLORMAPS['dbz'], -60, 40, hmax, times, first=False)
            
            vel_field = radar.fields['VEL']
            vel_min = vel_field.get('field_limit_lower', -10)
            vel_max = vel_field.get('field_limit_upper', 10)
            _plot_vpt_field_composite(display, axes[0,1], 'VEL', 'Velocity (m/s)', 
                           COLORMAPS['vel'], vel_min, vel_max, hmax, times, first=False)
            
            _plot_vpt_field_composite(display, axes[1,0], 'WIDTH', 'Spectrum Width (m/s)', 
                           COLORMAPS['width'], 0.01, 3, hmax, times, log_scale=True, first=False)
            
            _plot_vpt_field_composite(display, axes[1,1], 'LDR', 'LDR (dB)', 
                           COLORMAPS['ldr'], -35, 5, hmax, times, first=False)
    
    # Set title
    date_obj = datetime.datetime.strptime(datestr, '%Y%m%d')
    fig.suptitle(f'PICASSO MAN Dwell Composite VPT - {date_obj.strftime("%Y-%m-%d")}\n'
                 f'{len(dwell_data)} dwell segments', 
                 fontsize=14, fontweight='bold')
    
    # Save
    figname = f'ncas-mobile-ka-band-radar-1_cao_picasso_{datestr}_man_dwell_vpt_l1_{product_version}.png'
    plt.tight_layout()
    plt.savefig(os.path.join(vpt_figpath, figname), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"MAN dwell composite VPT saved to: {vpt_figpath}")
    print(f"  Combined {len(dwell_data)} dwell segments spanning "
          f"{dwell_data[0]['start_time'].strftime('%H:%M')} to "
          f"{dwell_data[-1]['times'][-1].strftime('%H:%M')} UTC")

def _plot_vpt_field_composite(display, ax, field_name, label, cmap, vmin, vmax, hmax, times, log_scale=False, first=True):
    """
    Helper function to plot VPT field for composite plots.
    
    Args:
        display: PyART RadarDisplay object
        ax: Matplotlib axis
        field_name: Name of field to plot
        label: Label for colorbar
        cmap: Colormap name
        vmin, vmax: Color limits
        hmax: Maximum height
        times: Array of datetime objects
        log_scale: Whether to use log scale
        first: If True, set up axes; if False, only add data
    """
    if log_scale:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    
    if first:
        # First dwell - create plot with labels
        display.plot_vpt(field_name, ax=ax, time_axis_flag=True, edges=True,
                        vmin=vmin, vmax=vmax, cmap=cmap, 
                        colorbar_label=label, norm=norm, shading='auto')
        
        ax.set_ylim(0, hmax)
        ax.set_xlabel('Time (UTC)')
        ax.set_ylabel('Height (km)')
        ax.set_title(label)
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    else:
        # Subsequent dwells - just add data without new colorbar
        display.plot_vpt(field_name, ax=ax, time_axis_flag=True, edges=True,
                        vmin=vmin, vmax=vmax, cmap=cmap, 
                        colorbar_flag=False, norm=norm, shading='auto')

def process_day(datestr, inpath, figpath, blflag=False):
    """
    Process all radar files for a given day.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Output directory path
        blflag: Boundary layer flag
    """
    print(f"\n{'='*60}")
    print(f"Processing PICASSO quicklooks for {datestr}")
    print(f"{'='*60}\n")
    
    # Find all files for this date
    date_path = os.path.join(inpath, datestr)
    if not os.path.exists(date_path):
        print(f"Date directory does not exist: {date_path}")
        return
    
    # Process RHI files
    rhi_files = sorted(glob.glob(os.path.join(date_path, '*_rhi_*.nc')))
    print(f"\nFound {len(rhi_files)} RHI files")
    for ncfile in rhi_files:
        try:
            make_picasso_rhi_plot(ncfile, figpath, blflag)
        except Exception as e:
            print(f"Error processing RHI file {ncfile}: {e}")
    
    # Process MAN files
    man_files = sorted(glob.glob(os.path.join(date_path, '*_man_*.nc')))
    print(f"\nFound {len(man_files)} MAN files")
    for ncfile in man_files:
        try:
            make_picasso_man_plot(ncfile, figpath, blflag)
        except Exception as e:
            print(f"Error processing MAN file {ncfile}: {e}")
    
    # Process VPT files
    vpt_files = sorted(glob.glob(os.path.join(date_path, '*_vpt_*.nc')))
    print(f"\nFound {len(vpt_files)} VPT files")
    for ncfile in vpt_files:
        try:
            make_picasso_vpt_plot(ncfile, figpath, blflag)
        except Exception as e:
            print(f"Error processing VPT file {ncfile}: {e}")
    
    # Create composite VPT plot from MAN dwell segments
    try:
        make_picasso_man_dwell_composite(datestr, inpath, figpath, blflag)
    except Exception as e:
        print(f"Error creating MAN dwell composite: {e}")
    
    # Create merged composite VPT plot (VPT files + MAN dwells)
    try:
        make_picasso_vpt_composite(datestr, inpath, figpath, blflag)
    except Exception as e:
        print(f"Error creating composite VPT plot: {e}")
    
    # Process PPI files (if any)
    ppi_files = sorted(glob.glob(os.path.join(date_path, '*_ppi_*.nc')))
    if ppi_files:
        print(f"\nFound {len(ppi_files)} PPI files")
        # PPI plots would be similar to RHI but plan view
        # (Implementation can be added if needed)
    
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
    datestr, inpath, outpath, blflag = parse_command_line()
    
    # Setup paths
    inpath, figpath = setup_paths(datestr, inpath, outpath)
    
    print(f"PICASSO Quicklook Generator v{VERSION}")
    print(f"Date: {datestr}")
    print(f"Input path: {inpath}")
    print(f"Output path: {figpath}")
    print(f"Boundary layer mode: {blflag}")
    
    # Process the day
    process_day(datestr, inpath, figpath, blflag)

if __name__ == "__main__":
    main()
