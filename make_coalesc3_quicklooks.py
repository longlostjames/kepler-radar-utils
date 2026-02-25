#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
make_coalesc3_quicklooks.py

Generate quicklook plots for COALESC3 campaign radar data.

This script creates visualization plots for NCAS Mobile Ka-band Radar data
from the COALESC3 (Combined Observations of the Atmospheric Boundary Layer to
study the Evolution of StratoCumulus) campaign at Weybourne Atmospheric Observatory.

It generates VPT (Vertical Pointing), RHI (Range Height Indicator), and 
VAD (Velocity Azimuth Display) plots for boundary layer cloud observations.

Usage:
    python make_coalesc3_quicklooks.py -d YYYYMMDD -i input_path -o output_path [-b]

Arguments:
    -d, --date:     Date string in YYYYMMDD format
    -i, --inpath:   Input directory containing CF-Radial files
    -o, --outpath:  Output directory for quicklook images
    -b:             Boundary layer flag (limits plots to 4km height)

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
Last modified: 24-02-2026
Version: 0.1
"""

import getopt
import sys
import os
import datetime
import glob
from pathlib import Path

import netCDF4 as nc4
import pyart
import numpy as np
import numpy.ma as ma
import cftime

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib import colors
import cmocean

# Configuration
VERSION = 0.1
TRACKING_TAG = 'AMF_07092016101810'
CAMPAIGN = 'coalesc3'

# Default paths
NCAS_OBS_PROC_PATH = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1'

# Plot configuration
COLORMAPS = {
    'dbz': 'HomeyerRainbow',
    'vel': 'balance',
    'ldr': 'SpectralExtended',
    'width': 'SpectralExtended'
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
        print("Usage: python make_coalesc3_quicklooks.py -d YYYYMMDD -i input_path -o output_path [-b]")
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
        inpath = os.path.join(NCAS_OBS_PROC_PATH, CAMPAIGN, 'L1_v1.0.1')
    
    # Output path for figures
    if outpath is None:
        figpath = os.path.join(inpath, 'quicklooks')
    else:
        figpath = outpath

    return inpath, figpath

def make_coalesc3_vpt_plot_day(datestr, inpath, figpath, blflag=False):
    """
    Generate daily vertical profiling time (VPT) plots.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Figure output directory path
        blflag: Boundary layer flag (True limits to 4km height)
    """
    if blflag:
        hmax = 4
    else:
        hmax = 12

    velmin = -5.0
    velmax = 5.0

    current_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    prev_date = current_date - datetime.timedelta(days=1)
    prevstr = prev_date.strftime('%Y%m%d')

    inpath_date = os.path.join(inpath, datestr)
    
    try:
        os.chdir(inpath_date)
        vpt_file = [os.path.join(inpath_date, f) for f in glob.glob(f'*{datestr}*vpt*.nc')][0]
    except (IndexError, FileNotFoundError):
        print(f"No VPT file found for {datestr}")
        return

    DS = nc4.Dataset(vpt_file)
    product_version = DS.product_version
    print(f'Product version = {product_version}')
    DS.close()

    RadarDS_VPT = pyart.io.read_cfradial(vpt_file)
    nsweeps = RadarDS_VPT.nsweeps

    # Try to load previous day's data for continuity
    nsweeps_prev = 0
    try:
        inpath_prev = os.path.join(inpath, prevstr)
        os.chdir(inpath_prev)
        vpt_file_prev = [os.path.join(inpath_prev, f) for f in glob.glob(f'*{prevstr}*vpt*.nc')][0]
        RadarDS_VPT_prev = pyart.io.read_cfradial(vpt_file_prev)
        nsweeps_prev = RadarDS_VPT_prev.nsweeps
    except:
        pass
    
    fig, ax = plt.subplots(4, 1, figsize=(12, 18), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=1/72, h_pad=1/72, hspace=0.2, wspace=0.2)

    # Process first sweep
    RadarSweepDS = RadarDS_VPT.extract_sweeps([0])
    print("Sweep = 0")
    gatefilter = pyart.correct.GateFilter(RadarSweepDS)
    gatefilter.exclude_below('SNR', -20)
    display = pyart.graph.RadarDisplay(RadarSweepDS)

    dtime = cftime.num2pydate(RadarDS_VPT.time['data'], RadarDS_VPT.time['units'])
    dt_min = dtime[0].replace(hour=0, minute=0, second=0)
    dt_max = dt_min + datetime.timedelta(days=1)
    time_str = dtime[0].strftime("%Y-%m-%d")
    
    # Plot DBZ
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "DBZ")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str}\n{field_name}"
    display.plot_vpt("DBZ", ax=ax[0], time_axis_flag=True, title=title, edges=False, 
                     gatefilter=gatefilter, vmin=-40, vmax=40, norm=None, 
                     filter_transitions=True, cmap=COLORMAPS['dbz'], 
                     colorbar_orient='horizontal')
    ax[0].set_ylim(0, hmax)
    ax[0].grid(True)

    # Plot VEL
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "VEL")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str}\n{field_name}"
    display.plot_vpt("VEL", ax=ax[1], time_axis_flag=True, title=title, edges=False, 
                     gatefilter=gatefilter, vmin=velmin, vmax=velmax, norm=None, 
                     cmap=COLORMAPS['vel'], colorbar_orient='horizontal')
    ax[1].set_ylim(0, hmax)
    ax[1].grid(True)

    # Plot WIDTH
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "WIDTH")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str}\n{field_name}"
    display.plot_vpt("WIDTH", ax=ax[2], time_axis_flag=True, title=title, edges=False, 
                     gatefilter=gatefilter, 
                     norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1), vmax=np.sqrt(1e1)), 
                     cmap=COLORMAPS['width'], colorbar_orient='horizontal')
    ax[2].set_ylim(0, hmax)
    ax[2].grid(True)

    # Plot LDR
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "LDR")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str}\n{field_name}"
    display.plot_vpt("LDR", ax=ax[3], time_axis_flag=True, title=title, edges=False, 
                     gatefilter=gatefilter, vmin=-35, vmax=5, norm=None, 
                     cmap=COLORMAPS['ldr'], colorbar_orient='horizontal')
    ax[3].set_ylim(0, hmax)
    ax[3].grid(True)

    # Set time limits
    ax[0].set_xlim(dt_min, dt_max)
    ax[1].set_xlim(dt_min, dt_max)
    ax[2].set_xlim(dt_min, dt_max)
    ax[3].set_xlim(dt_min, dt_max)

    ax[0].set_xlabel('Time (UTC)')
    ax[1].set_xlabel('Time (UTC)')
    ax[2].set_xlabel('Time (UTC)')
    ax[3].set_xlabel('Time (UTC)')

    # Process remaining sweeps
    for s in range(1, nsweeps):
        RadarSweepDS = RadarDS_VPT.extract_sweeps([s])
        print(f"Sweep = {s}")
        gatefilter = pyart.correct.GateFilter(RadarSweepDS)
        gatefilter.exclude_below('SNR', -20)
        display = pyart.graph.RadarDisplay(RadarSweepDS)
        display.plot_vpt("DBZ", ax=ax[0], time_axis_flag=True, edges=False, 
                         gatefilter=gatefilter, vmin=-40, vmax=40, norm=None, 
                         filter_transitions=True, cmap=COLORMAPS['dbz'], 
                         colorbar_flag=False, title_flag=False)
        display.plot_vpt("VEL", ax=ax[1], time_axis_flag=True, edges=False, 
                         gatefilter=gatefilter, vmin=velmin, vmax=velmax, norm=None, 
                         cmap=COLORMAPS['vel'], colorbar_flag=False, title_flag=False)
        display.plot_vpt("WIDTH", ax=ax[2], time_axis_flag=True, edges=False, 
                         gatefilter=gatefilter, 
                         norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1), vmax=np.sqrt(1e1)), 
                         cmap=COLORMAPS['width'], colorbar_flag=False, title_flag=False)
        display.plot_vpt("LDR", ax=ax[3], time_axis_flag=True, edges=False, 
                         gatefilter=gatefilter, vmin=-35, vmax=5, norm=None, 
                         cmap=COLORMAPS['ldr'], colorbar_flag=False, title_flag=False)

    # Add previous day's last sweep if available
    if nsweeps_prev > 0:
        for s in range(nsweeps_prev-1, nsweeps_prev):
            RadarSweepDS = RadarDS_VPT_prev.extract_sweeps([s])
            print(f"Previous day sweep = {s}")
            gatefilter = pyart.correct.GateFilter(RadarSweepDS)
            gatefilter.exclude_below('SNR', -20)
            display = pyart.graph.RadarDisplay(RadarSweepDS)
            display.plot_vpt("DBZ", ax=ax[0], time_axis_flag=True, edges=False, 
                             gatefilter=gatefilter, vmin=-40, vmax=40, norm=None, 
                             filter_transitions=True, cmap=COLORMAPS['dbz'], 
                             colorbar_flag=False, title_flag=False)
            display.plot_vpt("VEL", ax=ax[1], time_axis_flag=True, edges=False, 
                             gatefilter=gatefilter, vmin=velmin, vmax=velmax, norm=None, 
                             cmap=COLORMAPS['vel'], colorbar_flag=False, title_flag=False)
            display.plot_vpt("WIDTH", ax=ax[2], time_axis_flag=True, edges=False, 
                             gatefilter=gatefilter, 
                             norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1), vmax=np.sqrt(1e1)), 
                             cmap=COLORMAPS['width'], colorbar_flag=False, title_flag=False)
            display.plot_vpt("LDR", ax=ax[3], time_axis_flag=True, edges=False, 
                             gatefilter=gatefilter, vmin=-35, vmax=5, norm=None, 
                             cmap=COLORMAPS['ldr'], colorbar_flag=False, title_flag=False)
     
    ax[0].grid(True)
    ax[1].grid(True)
    ax[2].grid(True)
    ax[3].grid(True)

    ax[0].set_xlabel('Time (UTC)')
    ax[1].set_xlabel('Time (UTC)')
    ax[2].set_xlabel('Time (UTC)')
    ax[3].set_xlabel('Time (UTC)')

    figname = f'ncas-mobile-ka-band-radar-1_wao_{CAMPAIGN}_{datestr}_vpt_l1_{product_version}.png'
    
    if blflag:
        figname = figname.replace('.png', '_bl.png')

    figpath_vpt = os.path.join(figpath, 'vpt')
    if not os.path.isdir(figpath_vpt): 
        os.makedirs(figpath_vpt)

    plt.savefig(os.path.join(figpath_vpt, figname), dpi=300)
    print(f"Saved VPT plot: {figname}")
    plt.close()

def make_coalesc3_vad_plot_day(datestr, inpath, figpath, zlevels, blflag=False):
    """
    Generate daily VAD (Velocity Azimuth Display) wind plots.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Figure output directory path
        zlevels: Height levels for VAD retrieval (numpy array)
        blflag: Boundary layer flag (True limits to 4km height)
    """
    if blflag:
        hmax = 4
    else:
        hmax = 12

    inpath_date = os.path.join(inpath, datestr)

    try:
        os.chdir(inpath_date)
        vad_files = sorted([os.path.join(inpath_date, f) for f in glob.glob(f'*{datestr}*vad*.nc')])
    except (FileNotFoundError, OSError):
        print(f"No VAD files found for {datestr}")
        return
    
    if not vad_files:
        print(f"No VAD files found for {datestr}")
        return
    
    print(f"Found {len(vad_files)} VAD files")

    # Read first file to get product version
    DS = nc4.Dataset(vad_files[0])
    product_version = DS.product_version
    print(f'Product version = {product_version}')
    DS.close()

    # Process all VAD files
    vad_ray_index_start = []
    vad_ray_index_end = []
    u_allsweeps = []
    v_allsweeps = []
    dt_vad_start_all = []
    dt_vad_end_all = []
    elevation_angles = []  # Track all elevation angles
    
    for vad_file in vad_files:
        try:
            radar = pyart.io.read_cfradial(vad_file)
            
            # Track elevation angle from this file
            if radar.nsweeps > 0:
                elv_angle = radar.fixed_angle['data'][0]
                if elv_angle not in elevation_angles:
                    elevation_angles.append(elv_angle)
            
            # Process each sweep in the file
            for s in range(radar.nsweeps):
                sweep_start = radar.sweep_start_ray_index['data'][s]
                sweep_end = radar.sweep_end_ray_index['data'][s]
                
                dt_start = cftime.num2pydate(radar.time['data'][sweep_start], radar.time['units'])
                dt_end = cftime.num2pydate(radar.time['data'][sweep_end], radar.time['units'])
                
                dt_vad_start_all.append(dt_start)
                dt_vad_end_all.append(dt_end)
                
                # Extract the sweep and perform VAD
                radar_1sweep = radar.extract_sweeps([s])
                gatefilter = pyart.correct.GateFilter(radar_1sweep)
                gatefilter.exclude_below('SNR', 0)
                
                try:
                    vad = pyart.retrieve.vad_browning(radar_1sweep, "VEL", z_want=zlevels, gatefilter=gatefilter)
                    u_allsweeps.append(vad.u_wind)
                    v_allsweeps.append(vad.v_wind)
                except Exception as e:
                    print(f"VAD retrieval failed for {os.path.basename(vad_file)} sweep {s}: {e}")
                    # Add NaN arrays to maintain indexing
                    u_allsweeps.append(np.full_like(zlevels, np.nan, dtype=float))
                    v_allsweeps.append(np.full_like(zlevels, np.nan, dtype=float))
                    
        except Exception as e:
            print(f"Error processing VAD file {vad_file}: {e}")
            continue
    
    if len(u_allsweeps) == 0:
        print(f"No valid VAD data for {datestr}")
        return
    
    dt_vad_start = np.array(dt_vad_start_all)
    dt_vad_end = np.array(dt_vad_end_all)

    u_vel = np.array(u_allsweeps)
    v_vel = np.array(v_allsweeps)
    
    print(f"Processed {len(u_allsweeps)} sweeps")
    print(f"u_vel shape: {u_vel.shape}, zlevels shape: {zlevels.shape}")
    
    orientation = np.rad2deg(np.arctan2(-u_vel, -v_vel)) % 360
    speed = np.sqrt(u_vel**2 + v_vel**2)

    speed = ma.masked_where(speed > 100., speed)
    orientation = ma.masked_where(speed > 100., orientation)

    vad_duration = (dt_vad_end - dt_vad_start)
    dt_vad_mid = dt_vad_start + 0.5 * vad_duration 
    vad_duration[:] = datetime.timedelta(minutes=12)

    myFmt = mdates.DateFormatter('%H:%M')

    fig, ax = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=2/72, h_pad=2/72, hspace=0.2, wspace=0.2)

    dtime = cftime.num2pydate(radar.time['data'], radar.time['units'])

    dt_min = dtime[0].replace(hour=0, minute=0, second=0)
    dt_max = dt_min + datetime.timedelta(days=1)

    # Plot wind direction
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].grid(True)
    h1 = ax[0].pcolormesh([dt_vad_start[0] - 0.5*vad_duration[0], dt_vad_end[0] + 0.5*vad_duration[0]],
                          zlevels/1000., orientation[0:1, :-1].transpose(), 
                          cmap='twilight_shifted', vmin=0, vmax=360)
    ax[0].set_ylim(0, 12)
    ax[0].set_ylabel('Distance above radar [km]')
    cb0 = plt.colorbar(h1, ax=ax[0], orientation='horizontal', shrink=0.8)
    cb0.ax.set_xlabel("Wind from direction (deg)")
    ax[0].set_facecolor('white')
    
    for s in range(1, len(dt_vad_start)):
        ax[0].pcolormesh([dt_vad_start[s] - 0.5*vad_duration[s], dt_vad_end[s] + 0.5*vad_duration[s]],
                         zlevels/1000., orientation[s:s+1, :-1].transpose(), 
                         cmap='twilight_shifted', vmin=0, vmax=360)

    ax[0].set_xlabel("Time (UTC)")

    # Plot wind speed
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].grid(True)
    h2 = ax[1].pcolormesh([dt_vad_start[0] - 0.5*vad_duration[0], dt_vad_end[0] + 0.5*vad_duration[0]],
                          zlevels/1000., speed[0:1, :-1].transpose(), 
                          cmap='viridis', vmin=0, vmax=50)
    ax[1].grid(True)
    ax[1].set_ylim(0, 12)
    ax[1].set_ylabel('Distance above radar (km)')
    cb1 = plt.colorbar(h2, ax=ax[1], orientation='horizontal', shrink=0.8)
    cb1.ax.set_xlabel("Wind speed (m/s)")
    ax[1].set_facecolor('gainsboro')

    for s in range(1, len(dt_vad_start)):
        ax[1].pcolormesh([dt_vad_start[s] - 0.5*vad_duration[s], dt_vad_end[s] + 0.5*vad_duration[s]], 
                         zlevels/1000., speed[s:s+1, :-1].transpose(), 
                         cmap='viridis', vmin=0, vmax=50)

    nlevels = zlevels.shape[0]
    nsweeps = len(dt_vad_mid)  # Use actual number of processed sweeps

    x = np.tile(dt_vad_mid, [nlevels, 1]).transpose()
    y = np.tile(zlevels/1000., [nsweeps, 1])
    
    ax[1].barbs(x[:, 5::10], y[:, 5::10], u_vel[:, 5::10], v_vel[:, 5::10], 
                np.sqrt(u_vel[:, 5::10]**2 + v_vel[:, 5::10]**2), 
                sizes=dict(emptybarb=0.), length=6, cmap='Grays', clim=[0, 50])
    ax[1].set_facecolor('gainsboro')
    ax[1].set_xlabel("Time (UTC)")

    dtime = cftime.num2pydate(radar.time['data'], radar.time['units'])[0]
    time_str = dtime.strftime("%Y-%m-%d")
    
    # Create elevation string for title
    if len(elevation_angles) == 1:
        elev_str = f"at elevation {elevation_angles[0]:.1f}°"
    else:
        elev_list = ", ".join([f"{e:.1f}°" for e in sorted(elevation_angles)])
        elev_str = f"at elevations {elev_list}"
    
    title0 = f"{pyart.graph.common.generate_radar_name(radar)} {time_str}\n"
    title0 += f"Wind direction from VAD {elev_str}"
    title1 = f"{pyart.graph.common.generate_radar_name(radar)} {time_str}\n"
    title1 += f"Wind speed from VAD {elev_str}"

    ax[0].set_xlim(dt_min, dt_max)
    ax[1].set_xlim(dt_min, dt_max)
    ax[0].grid(True)
    ax[1].grid(True)

    ax[0].set_title(title0)
    ax[1].set_title(title1)

    figpath_vad = os.path.join(figpath, 'vad')
    if not os.path.isdir(figpath_vad): 
        os.makedirs(figpath_vad)

    figname = f'ncas-mobile-ka-band-radar-1_wao_{CAMPAIGN}_{datestr}_vad_l1_{product_version}.png'
    plt.savefig(os.path.join(figpath_vad, figname), dpi=300)
    print(f"Saved VAD plot: {figname}")
    plt.close()

def make_coalesc3_rhi_plots_day(datestr, inpath, figpath, blflag=False):
    """
    Generate RHI (Range Height Indicator) plots for all RHI scans in a day.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        figpath: Figure output directory path
        blflag: Boundary layer flag (True limits to 4km height)
    """
    if blflag:
        hmax = 4
        xmin = -20
        xmax = 20
    else:
        hmax = 12
        xmin = -25
        xmax = 25

    inpath_date = os.path.join(inpath, datestr)
    
    try:
        os.chdir(inpath_date)
        rhi_files = sorted(glob.glob(f'*{datestr}*rhi*.nc'))
    except (FileNotFoundError, OSError):
        print(f"No RHI files found for {datestr}")
        return
    
    if not rhi_files:
        print(f"No RHI files found for {datestr}")
        return
    
    print(f"Found {len(rhi_files)} RHI files")
    
    # Create output directory
    figpath_rhi = os.path.join(figpath, 'rhi', datestr)
    if not os.path.isdir(figpath_rhi):
        os.makedirs(figpath_rhi)
    
    # Process each RHI file
    for rhi_file in rhi_files:
        try:
            rhi_path = os.path.join(inpath_date, rhi_file)
            
            # Read the file
            DS = nc4.Dataset(rhi_path)
            product_version = DS.product_version
            DS.close()
            
            RadarDS = pyart.io.read_cfradial(rhi_path)
            
            # Get time and azimuth
            dtime0 = cftime.num2pydate(RadarDS.time['data'][0], RadarDS.time['units'])
            dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S")
            rhi_az = RadarDS.get_azimuth(0)[0]
            
            # Create gatefilter
            gatefilter = pyart.correct.GateFilter(RadarDS)
            gatefilter.exclude_below('SNR', -20)
            
            # Create display
            display = pyart.graph.RadarDisplay(RadarDS)
            
            # Create figure
            fig, ax = plt.subplots(2, 2, figsize=(15, 15))
            
            # Plot DBZ
            display.plot_rhi("DBZ", ax=ax[0, 0], sweep=0, vmin=-40, vmax=40, 
                           norm=None, gatefilter=gatefilter,
                           cmap=COLORMAPS['dbz'], colorbar_orient='horizontal',
                           reverse_xaxis=False, filter_transitions=True)
            ax[0, 0].set_ylim(0, hmax)
            ax[0, 0].set_xlim(xmin, xmax)
            ax[0, 0].grid(True)
            ax[0, 0].set_aspect('equal')
            
            # Plot VEL
            display.plot_rhi("VEL", ax=ax[1, 0], sweep=0, vmin=-10, vmax=10,
                           gatefilter=gatefilter, norm=None, cmap=COLORMAPS['vel'],
                           colorbar_orient='horizontal', reverse_xaxis=False)
            ax[1, 0].set_ylim(0, hmax)
            ax[1, 0].set_xlim(xmin, xmax)
            ax[1, 0].grid(True)
            ax[1, 0].set_aspect('equal', 'box')
            
            # Plot WIDTH
            display.plot_rhi("WIDTH", ax=ax[1, 1], sweep=0,
                           norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1), vmax=np.sqrt(1e1)),
                           gatefilter=gatefilter, cmap=COLORMAPS['width'],
                           colorbar_orient='horizontal', reverse_xaxis=False)
            ax[1, 1].set_ylim(0, hmax)
            ax[1, 1].set_xlim(xmin, xmax)
            ax[1, 1].grid(True)
            ax[1, 1].set_aspect('equal', 'box')
            
            # Plot LDR
            display.plot_rhi("LDR", ax=ax[0, 1], sweep=0, vmin=-35, vmax=5,
                           norm=None, gatefilter=gatefilter, cmap=COLORMAPS['ldr'],
                           colorbar_orient='horizontal', reverse_xaxis=False)
            ax[0, 1].set_ylim(0, hmax)
            ax[0, 1].set_xlim(xmin, xmax)
            ax[0, 1].grid(True)
            ax[0, 1].set_aspect('equal', 'box')
            
            # Save figure
            figname = f'ncas-mobile-ka-band-radar-1_wao_{CAMPAIGN}_{dtime0_str}_rhi_az{rhi_az:0.2f}_l1_{product_version}.png'
            plt.savefig(os.path.join(figpath_rhi, figname), dpi=300)
            plt.close()
            
            print(f"Saved RHI plot: {figname}")
            
        except Exception as e:
            print(f"Error processing RHI file {rhi_file}: {e}")
            continue

def main():
    """Main execution function."""
    datestr, inpath, outpath, blflag = parse_command_line()
    inpath, figpath = setup_paths(datestr, inpath, outpath)
    
    print(f"COALESC3 Quicklooks Generation")
    print(f"Date: {datestr}")
    print(f"Input path: {inpath}")
    print(f"Output path: {figpath}")
    print(f"Boundary layer mode: {blflag}")
    print("-" * 60)
    
    # Generate VPT quicklooks
    print("Generating VPT plots...")
    try:
        make_coalesc3_vpt_plot_day(datestr, inpath, figpath, blflag=blflag)
    except Exception as e:
        print(f"Error generating VPT plots: {e}")
    
    # Generate RHI quicklooks
    print("Generating RHI plots...")
    try:
        make_coalesc3_rhi_plots_day(datestr, inpath, figpath, blflag=blflag)
    except Exception as e:
        print(f"Error generating RHI plots: {e}")
        import traceback
        traceback.print_exc()
    
    # Generate VAD quicklooks
    print("Generating VAD plots...")
    try:
        zlevels = np.arange(100, 15000, 100)  # height above radar in meters
        make_coalesc3_vad_plot_day(datestr, inpath, figpath, zlevels, blflag=blflag)
    except Exception as e:
        print(f"Error generating VAD plots: {e}")
        import traceback
        traceback.print_exc()
    
    print("-" * 60)
    print("Quicklook generation complete!")

if __name__ == "__main__":
    main()
