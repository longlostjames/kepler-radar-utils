#!/usr/bin/env python
# coding: utf-8

"""
campaign_processing.py

Campaign-specific processing functions for Kepler MIRA-35 radar data

This module provides processing functions for different campaigns:
- WOEST (WESCON Observing the Evolving Structures of Turbulence)
- CCREST (Characterising CiRrus and icE cloud acrosS the specTrum - Microwave)
- COBALT (Contrail OBservations And Lifecycle Tracking)
- KASBEX (Ka- and S-Band EXperiment)

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
Last modified: 10-08-2025
Version: 1.0.0
"""

from typing import List, Dict, Tuple, Optional, Any
import datetime
import os
import glob
import yaml
import pyart
import numpy as np
from kepler_utils import (
    read_mira35_mmclx, 
    multi_mmclx2cfrad,
    cfradial_add_ncas_metadata,
    update_history_attribute,
    find_mmclxfiles,
    find_mmclx_rhi_files,
    find_mmclx_ppi_files,
    split_monotonic_sequence
)

def process_kepler_woest_day_step1(
    datestr: str, 
    inpath: str, 
    outpath: str, 
    yaml_project_file: str, 
    yaml_instrument_file: str, 
    azimuth_offset: float = -6.85, 
    revised_northangle: float = 302.15, 
    gzip_flag: bool = False, 
    data_version: str = "1.0.0"
) -> None:
    """
    Process WOEST campaign data for a single day - Step 1.
    
    This function processes RHI scans from the WOEST campaign, converting
    mmclx files to CF-Radial format and organizing by azimuth sectors.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing mmclx files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        azimuth_offset: Azimuth offset correction (default: -6.85)
        revised_northangle: North angle correction (default: 302.15)
        gzip_flag: Whether input files are gzip compressed
        data_version: Data version string
    """
    print(f"Processing WOEST day: {datestr}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Define WOEST azimuth sectors
    woest_sectors = [
        (315, 345),  # Northwest sector
        (345, 15),   # North sector  
        (15, 45),    # Northeast sector
        (45, 75),    # East sector
        (75, 105),   # Southeast sector
        (105, 135),  # South sector
        (135, 165),  # Southwest sector
        (165, 195),  # West sector
    ]
    
    # Time range for the day
    start_time = f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:8]} 00:00:00"
    end_time = f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:8]} 23:59:59"
    
    # Process each azimuth sector
    for sector_idx, (azim_min, azim_max) in enumerate(woest_sectors):
        print(f"Processing sector {sector_idx + 1}: {azim_min}°-{azim_max}°")
        
        # Find RHI files in this azimuth range
        rhi_files = find_mmclx_rhi_files(
            start_time, end_time, azim_min, azim_max, inpath,
            gzip_flag=gzip_flag, azimuth_offset=azimuth_offset,
            revised_northangle=revised_northangle
        )
        
        if not rhi_files:
            print(f"No RHI files found for sector {sector_idx + 1}")
            continue
        
        # Process files in this sector
        for rhi_file in rhi_files:
            try:
                _process_single_woest_file(
                    rhi_file, outdir, yaml_project_file, yaml_instrument_file,
                    revised_northangle, gzip_flag, data_version
                )
            except Exception as e:
                print(f"Error processing {rhi_file}: {e}")

def _process_single_woest_file(
    infile: str,
    outdir: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    revised_northangle: float,
    gzip_flag: bool,
    data_version: str
) -> None:
    """
    Process a single WOEST file.
    
    Args:
        infile: Input mmclx file path
        outdir: Output directory
        yaml_project_file: Project YAML file
        yaml_instrument_file: Instrument YAML file
        revised_northangle: North angle correction
        gzip_flag: Whether file is gzip compressed
        data_version: Data version string
    """
    # Read radar data
    radar_ds = read_mira35_mmclx(infile, gzip_flag=gzip_flag, revised_northangle=revised_northangle)
    
    # Generate output filename
    base_filename = os.path.basename(infile).replace('.mmclx', '').replace('.gz', '')
    outfile = os.path.join(outdir, f"{base_filename}_l1_v{data_version}.nc")
    
    # Write CF-Radial file
    pyart.io.write_cfradial(outfile, radar_ds, format='NETCDF4', time_reference=True)
    
    # Add NCAS metadata
    tracking_tag = 'AMOF_20220301000000'  # WOEST tracking tag
    cfradial_add_ncas_metadata(outfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version)
    
    # Update history
    update_history_attribute(outfile, f"WOEST Step 1: Convert mmclx to CF-Radial, revised_northangle={revised_northangle}")
    
    print(f"Created: {outfile}")

def process_kepler_ccrest_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    azimuth_offset: float = 0.0,
    revised_northangle: float = 55.7,
    gzip_flag: bool = False,
    data_version: str = "1.0.0"
) -> None:
    """
    Process CCREST campaign data for a single day - Step 1.
    
    This function processes both RHI and PPI scans from the CCREST campaign.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing mmclx files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        azimuth_offset: Azimuth offset correction (default: 0.0)
        revised_northangle: North angle correction (default: 55.7)
        gzip_flag: Whether input files are gzip compressed
        data_version: Data version string
    """
    print(f"Processing CCREST day: {datestr}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Time range for the day
    start_time = f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:8]} 00:00:00"
    end_time = f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:8]} 23:59:59"
    
    tracking_tag = 'AMOF_20230401000000'  # CCREST tracking tag
    
    # Process RHI files
    print("Processing RHI files...")
    rhi_files = find_mmclxfiles(start_time, end_time, "rhi", inpath, gzip_flag=gzip_flag)
    
    for rhi_file in rhi_files:
        try:
            _process_single_ccrest_file(
                rhi_file, outdir, yaml_project_file, yaml_instrument_file,
                revised_northangle, gzip_flag, data_version, tracking_tag
            )
        except Exception as e:
            print(f"Error processing RHI {rhi_file}: {e}")
    
    # Process PPI files
    print("Processing PPI files...")
    ppi_files = find_mmclxfiles(start_time, end_time, "ppi", inpath, gzip_flag=gzip_flag)
    
    for ppi_file in ppi_files:
        try:
            _process_single_ccrest_file(
                ppi_file, outdir, yaml_project_file, yaml_instrument_file,
                revised_northangle, gzip_flag, data_version, tracking_tag
            )
        except Exception as e:
            print(f"Error processing PPI {ppi_file}: {e}")

def _process_single_ccrest_file(
    infile: str,
    outdir: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    revised_northangle: float,
    gzip_flag: bool,
    data_version: str,
    tracking_tag: str
) -> None:
    """
    Process a single CCREST file.
    """
    # Read radar data
    radar_ds = read_mira35_mmclx(infile, gzip_flag=gzip_flag, revised_northangle=revised_northangle)
    
    # Generate output filename
    base_filename = os.path.basename(infile).replace('.mmclx', '').replace('.gz', '')
    outfile = os.path.join(outdir, f"{base_filename}_l1_v{data_version}.nc")
    
    # Write CF-Radial file
    pyart.io.write_cfradial(outfile, radar_ds, format='NETCDF4', time_reference=True)
    
    # Add NCAS metadata
    cfradial_add_ncas_metadata(outfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version)
    
    # Update history
    update_history_attribute(outfile, f"CCREST Step 1: Convert mmclx to CF-Radial, revised_northangle={revised_northangle}")
    
    print(f"Created: {outfile}")

def process_kepler_cobalt_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    azimuth_offset: float = 0.0,
    revised_northangle: float = 55.9,
    gzip_flag: bool = True,
    data_version: str = "1.0.3",
    single_sweep: bool = False,
    tracking_tag: str = 'AMOF_20231120125118',
    no_vpt: bool = False  # Add no_vpt argument
) -> None:
    """
    Process COBALT campaign data for a single day - Step 1.
    """
    print(f"Processing COBALT day: {datestr}")
    print(f"Single sweep mode: {single_sweep}")
    print(f"Data version: {data_version}")
    print(f"Tracking tag: {tracking_tag}")
    print(f"No VPT processing: {no_vpt}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Define the start and end times for the loop
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    end_date = start_date + datetime.timedelta(days=1)
    
    print(f"start_date = {start_date} end_date = {end_date}")
    
    # Skip VPT processing if no_vpt is True
    if no_vpt:
        print("Skipping VPT file processing (--no-vpt option enabled)")
    else:
        print("Processing VPT files...")
        try:
            vpt_files = find_mmclxfiles(
                start_date.strftime('%Y-%m-%d %H:%M:%S'),
                end_date.strftime('%Y-%m-%d %H:%M:%S'),
                'vert', 
                inpath,
                gzip_flag=gzip_flag
            )
            print(vpt_files)
            print(f"Found {len(vpt_files)} VPT files")
            
            if len(vpt_files) > 0:
                print("VPT files will always create multi-sweep output (ignoring single_sweep flag)")
                RadarDS_VPT = multi_mmclx2cfrad(
                    vpt_files, outdir, scan_name='VPT', gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset, 
                    tracking_tag=tracking_tag,
                    campaign='cobalt', 
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,  # ALWAYS FALSE FOR VPT
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(vpt_files)} VPT files into multi-sweep dataset")
        except Exception as e:
            print(f"VPT problem: {e}")
            pass
    
    # Process RHI files (use single_sweep setting)
    print("Processing RHI files...")
    rhi_files = find_mmclx_rhi_files(
        start_date.strftime('%Y-%m-%d %H:%M:%S'), 
        end_date.strftime('%Y-%m-%d %H:%M:%S'), 
        0, 370, 
        inpath,
        gzip_flag=gzip_flag, 
        azimuth_offset=azimuth_offset,
        revised_northangle=revised_northangle
    )
    
    print(rhi_files)
    print(f"Found {len(rhi_files)} RHI files")
    
    if len(rhi_files) > 0:
        if single_sweep:
            # Process each RHI file separately
            for f in rhi_files:
                try:
                    RadarDS_RHI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='RHI', gzip_flag=gzip_flag,
                        azimuth_offset=azimuth_offset, 
                        tracking_tag=tracking_tag,
                        campaign='cobalt', 
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single RHI file: {f}")
                except Exception as e:
                    print(f"Error processing RHI file {f}: {e}")
        else:
            # Process all RHI files together (original behavior)
            try:
                RadarDS_RHI = multi_mmclx2cfrad(
                    rhi_files, outdir, scan_name='RHI', gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset, 
                    tracking_tag=tracking_tag,
                    campaign='cobalt', 
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(rhi_files)} RHI files into combined dataset")
            except Exception as e:
                print(f"Error processing RHI files: {e}")
    
    # Process VAD files for whole day - ALWAYS MULTI-SWEEP  
    print("Processing VAD files...")
    try:
        from kepler_utils import find_mmclx_vad_files
        vad_files = find_mmclx_vad_files(
            start_date.strftime('%Y-%m-%d %H:%M:%S'),
            end_date.strftime('%Y-%m-%d %H:%M:%S'),
            80, 90, 
            inpath,
            gzip_flag=gzip_flag
        )
        print(f"Found {len(vad_files)} VAD files")
        
        if len(vad_files) > 0:
            print("VAD files will always create multi-sweep output (ignoring single_sweep flag)")
            RadarDS_VAD = multi_mmclx2cfrad(
                vad_files, outdir, scan_name='VAD', gzip_flag=gzip_flag,
                azimuth_offset=azimuth_offset, 
                tracking_tag=tracking_tag,
                campaign='cobalt', 
                revised_northangle=revised_northangle,
                data_version=data_version,
                single_sweep=False,  # ALWAYS FALSE FOR VAD
                yaml_project_file=yaml_project_file,
                yaml_instrument_file=yaml_instrument_file
            )
            print(f"Processed {len(vad_files)} VAD files into multi-sweep dataset")
    except Exception as e:
        print(f"VAD problem: {e}")
        pass
    
    print(f"Completed COBALT processing for {datestr}")

def process_kepler_kasbex_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    azimuth_offset: float = 0.0,
    revised_northangle: float = 55.7,
    gzip_flag: bool = False,
    data_version: str = "1.0.0",
    single_sweep: bool = False,
    tracking_tag: str = 'AMOF_20250508133639',
    min_elevation: float = 0.0,
    max_elevation: float = 89.5,
    no_vpt: bool = False,  # Add no_vpt argument
    max_age: Optional[float] = None  # Add max_age argument
) -> None:
    """
    Process KASBEX campaign data for a single day - Step 1.
    Now handles VPT files with different range gate dimensions.
    """
    print(f"Processing KASBEX day: {datestr}")
    print(f"Using revised_northangle: {revised_northangle}° (from YAML configuration)")
    print(f"Single sweep mode: {single_sweep}")
    print(f"Data version: {data_version}")
    print(f"Tracking tag: {tracking_tag}")
    print(f"Input path: {inpath}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Define the start and end times for the loop
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    end_date = start_date + datetime.timedelta(days=1)
    
    print(f"start_date = {start_date} end_date = {end_date}")
    
    # Check if inpath already includes the date (i.e., when using --latest)
    path_parts = inpath.rstrip('/').split('/')
    if datestr in path_parts or 'latest' in path_parts:
        # inpath already includes date/latest structure
        search_path = inpath
        print(f"Using direct path for file search: {search_path}")
        is_latest_mode = 'latest' in path_parts
    else:
        # Traditional structure - need to add date
        search_path = os.path.join(inpath, datestr)
        print(f"Using traditional date-based path: {search_path}")
        is_latest_mode = False
    
    # Process PPI files - prefer gzipped in latest mode
    print(f"Processing PPI files (elevation range: {min_elevation}° to {max_elevation}°)...")
    
    if is_latest_mode:
        print("Latest mode: Preferring gzipped PPI files")
        ppi_files = find_files_prefer_gzip(
        inpath, '*ppi*.mmclx*', force_gzip=gzip_flag, max_age_hours=max_age, use_max_age=True
        )
        #ppi_files = find_files_prefer_gzip(search_path, '*ppi*.mmclx*', force_gzip=True)
        # All files should be gzipped in this mode
        ppi_gzip_flag = True
    else:
        # Traditional mode - use glob and global gzip_flag
        import glob
        ppi_pattern = os.path.join(search_path, '*ppi*.mmclx*')
        ppi_files = glob.glob(ppi_pattern)
        ppi_files.sort()
        ppi_gzip_flag = gzip_flag
    
    print(f"Found {len(ppi_files)} PPI files")
    if ppi_files:
        print(f"Sample PPI files: {[os.path.basename(f) for f in ppi_files[:3]]}")
        if is_latest_mode:
            print(f"Using gzip mode for PPI files: {ppi_gzip_flag}")
    
    if len(ppi_files) > 0:
        if single_sweep:
            # Process each PPI file separately
            for f in ppi_files:
                try:
                    RadarDS_PPI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='PPI', gzip_flag=ppi_gzip_flag,
                        azimuth_offset=azimuth_offset, 
                        tracking_tag=tracking_tag,
                        campaign='kasbex', 
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single PPI file: {os.path.basename(f)}")
                except Exception as e:
                    print(f"Error processing PPI file {f}: {e}")
        else:
            # Process all PPI files together
            try:
                RadarDS_PPI = multi_mmclx2cfrad(
                    ppi_files, outdir, scan_name='PPI', gzip_flag=ppi_gzip_flag,
                    azimuth_offset=azimuth_offset, 
                    tracking_tag=tracking_tag,
                    campaign='kasbex', 
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(ppi_files)} PPI files into combined dataset")
            except Exception as e:
                print(f"Error processing PPI files: {e}")
    
    # Process RHI files - prefer gzipped in latest mode
    print("Processing RHI files...")
    
    if is_latest_mode:
        print("Latest mode: Preferring gzipped RHI files")
        rhi_files = find_files_prefer_gzip(
            inpath, '*rhi*.mmclx*', force_gzip=gzip_flag, max_age_hours=max_age, use_max_age=True
        )
        print(f"Found {len(rhi_files)} RHI files")
        #rhi_files = find_files_prefer_gzip(search_path, '*rhi*.mmclx*', force_gzip=True)
        # All files should be gzipped in this mode
        rhi_gzip_flag = True
    else:
        # Traditional mode - use glob and global gzip_flag
        import glob
        rhi_pattern = os.path.join(search_path, '*rhi*.mmclx*')
        rhi_files = glob.glob(rhi_pattern)
        rhi_files.sort()
        rhi_gzip_flag = gzip_flag
    
    print(f"Found {len(rhi_files)} RHI files")
    if rhi_files:
        print(f"Sample RHI files: {[os.path.basename(f) for f in rhi_files[:3]]}")
        if is_latest_mode:
            print(f"Using gzip mode for RHI files: {rhi_gzip_flag}")
    
    if len(rhi_files) > 0:
        if single_sweep:
            # Process each RHI file separately
            for f in rhi_files:
                try:
                    RadarDS_RHI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='RHI', gzip_flag=rhi_gzip_flag,
                        azimuth_offset=azimuth_offset, 
                        tracking_tag=tracking_tag,
                        campaign='kasbex', 
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single RHI file: {os.path.basename(f)}")
                except Exception as e:
                    print(f"Error processing RHI file {f}: {e}")
        else:
            # Process all RHI files together
            try:
                RadarDS_RHI = multi_mmclx2cfrad(
                    rhi_files, outdir, scan_name='RHI', gzip_flag=rhi_gzip_flag,
                    azimuth_offset=azimuth_offset, 
                    tracking_tag=tracking_tag,
                    campaign='kasbex', 
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(rhi_files)} RHI files into combined dataset")
            except Exception as e:
                print(f"Error processing RHI files: {e}")
    
    # Process VPT files with range gate grouping - allow both gzipped and uncompressed
    if no_vpt:
        print("Skipping VPT file processing (--no-vpt option enabled)")
    else:
        print("Processing VPT files...")
        try:
            if is_latest_mode:
                print("Latest mode: Allowing both gzipped and uncompressed VPT files")
                # Look for both compressed and uncompressed VPT files, but don't force gzip
                vpt_files = find_files_prefer_gzip(search_path, '*vert*.mmclx*', force_gzip=False)
            
                # Separate by compression status for processing
                vpt_files_gz = [f for f in vpt_files if f.endswith('.gz')]
                vpt_files_normal = [f for f in vpt_files if not f.endswith('.gz')]
            
                print(f"Found {len(vpt_files_gz)} gzipped VPT files")
                print(f"Found {len(vpt_files_normal)} normal VPT files")
                print(f"Total VPT files: {len(vpt_files)}")
            
                if vpt_files_gz:
                    print(f"Sample gzipped VPT files: {[os.path.basename(f) for f in vpt_files_gz[:3]]}")
                if vpt_files_normal:
                    print(f"Sample normal VPT files: {[os.path.basename(f) for f in vpt_files_normal[:3]]}")
            
                # Process gzipped VPT files with range gate grouping
                if len(vpt_files_gz) > 0:
                    print("Processing gzipped VPT files (multi-sweep output)")
                    vpt_files_gz.sort()
                
                    # Group by range gates
                    gz_file_groups = group_files_by_range_gates(vpt_files_gz, gzip_flag=True)
                
                    for n_gates, file_group in gz_file_groups.items():
                        if n_gates == 'error':
                            print(f"Skipping {len(file_group)} problematic gzipped VPT files")
                            continue
                        
                        print(f"Processing {len(file_group)} gzipped VPT files with {n_gates} range gates")
                        try:
                            RadarDS_VPT_gz = multi_mmclx2cfrad(
                                file_group, outdir, scan_name='VPT', gzip_flag=True,
                                azimuth_offset=azimuth_offset, 
                                tracking_tag=tracking_tag,
                                campaign='kasbex', 
                                revised_northangle=revised_northangle,
                                data_version=data_version,
                                single_sweep=False,  # Always FALSE for VPT
                                yaml_project_file=yaml_project_file,
                                yaml_instrument_file=yaml_instrument_file
                            )
                            print(f"Successfully processed {len(file_group)} gzipped VPT files ({n_gates} gates)")
                        except Exception as e:
                            print(f"Error processing gzipped VPT files with {n_gates} gates: {e}")
            
                # Process normal VPT files with range gate grouping
                if len(vpt_files_normal) > 0:
                    print("Processing uncompressed VPT files (multi-sweep output)")
                    vpt_files_normal.sort()
                
                    # Group by range gates
                    normal_file_groups = group_files_by_range_gates(vpt_files_normal, gzip_flag=False)
                
                    for n_gates, file_group in normal_file_groups.items():
                        if n_gates == 'error':
                            print(f"Skipping {len(file_group)} problematic normal VPT files")
                            continue
                        
                        print(f"Processing {len(file_group)} normal VPT files with {n_gates} range gates")
                        try:
                            RadarDS_VPT_normal = multi_mmclx2cfrad(
                                file_group, outdir, scan_name='VPT', gzip_flag=False,
                                azimuth_offset=azimuth_offset, 
                                tracking_tag=tracking_tag,
                                campaign='kasbex', 
                                revised_northangle=revised_northangle,
                                data_version=data_version,
                                single_sweep=False,  # Always FALSE for VPT
                                yaml_project_file=yaml_project_file,
                                yaml_instrument_file=yaml_instrument_file
                            )
                            print(f"Successfully processed {len(file_group)} normal VPT files ({n_gates} gates)")
                        except Exception as e:
                            print(f"Error processing normal VPT files with {n_gates} gates: {e}")
        
            else:
                # Traditional mode - use global gzip_flag for all VPT files with range gate grouping
                import glob
                vpt_pattern = os.path.join(search_path, '*vert*.mmclx*')
                vpt_files = glob.glob(vpt_pattern)
                vpt_files.sort()
            
                print(f"Found {len(vpt_files)} VPT files")
                if vpt_files:
                    print(f"Sample VPT files: {[os.path.basename(f) for f in vpt_files[:3]]}")
            
                if len(vpt_files) > 0:
                    print("Processing VPT files (multi-sweep output) with range gate grouping")
                
                    # Group by range gates
                    file_groups = group_files_by_range_gates(vpt_files, gzip_flag=gzip_flag)
                
                    for n_gates, file_group in file_groups.items():
                        if n_gates == 'error':
                            print(f"Skipping {len(file_group)} problematic VPT files")
                            continue
                        
                        print(f"Processing {len(file_group)} VPT files with {n_gates} range gates")
                        try:
                            RadarDS_VPT = multi_mmclx2cfrad(
                                file_group, outdir, scan_name='VPT', gzip_flag=gzip_flag,
                                azimuth_offset=azimuth_offset, 
                                tracking_tag=tracking_tag,
                                campaign='kasbex', 
                                revised_northangle=revised_northangle,
                                data_version=data_version,
                                single_sweep=False,  # Always FALSE for VPT
                                yaml_project_file=yaml_project_file,
                                yaml_instrument_file=yaml_instrument_file
                            )
                            print(f"Successfully processed {len(file_group)} VPT files ({n_gates} gates)")
                        except Exception as e:
                            print(f"Error processing VPT files with {n_gates} gates: {e}")
            
        except Exception as e:
            print(f"VPT problem: {e}")
            import traceback
            traceback.print_exc()
            pass
    
    print(f"Completed KASBEX processing for {datestr} with north_angle={revised_northangle}°")

def find_files_prefer_gzip(search_path, pattern, force_gzip=False, max_age_hours=None, use_max_age=False):
    """
    Find files matching pattern, preferring gzipped versions over uncompressed ones,
    and optionally filter by maximum age using timestamps from filenames.

    Args:
        search_path (str): Directory to search in
        pattern (str): Glob pattern for files (e.g., '*ppi*.mmclx*')
        force_gzip (bool): If True, only return gzipped files
        max_age_hours (float): Maximum file age in hours (None for no limit)
        use_max_age (bool): Whether to apply max-age filtering (only if --latest is used)

    Returns:
        list: List of files with gzipped versions preferred and filtered by age
    """
    import glob
    import os
    import datetime
    import re

    # Find all matching files
    all_files = glob.glob(os.path.join(search_path, pattern))
    all_files.sort()

    # Apply max-age filtering only if use_max_age is True
    if use_max_age and max_age_hours is not None:
        current_time = datetime.datetime.utcnow()
        max_age_seconds = max_age_hours * 3600
        filtered_files = []
        time_pattern = r'(\d{8})-(\d{6})'  # Match YYYYMMDD-HHMMSS in filenames

        for file_path in all_files:
            try:
                filename = os.path.basename(file_path)
                match = re.search(time_pattern, filename)
                if match:
                    date_str = match.group(1)  # YYYYMMDD
                    time_str = match.group(2)  # HHMMSS
                    file_time = datetime.datetime.strptime(f"{date_str}-{time_str}", "%Y%m%d-%H%M%S")
                    file_age_seconds = (current_time - file_time).total_seconds()

                    if file_age_seconds <= max_age_seconds:
                        filtered_files.append(file_path)
                    else:
                        print(f"Skipping old file: {filename} (age: {file_age_seconds / 3600:.1f} hours)")
                else:
                    print(f"Skipping file with no timestamp: {filename}")
            except Exception as e:
                print(f"Error checking file age for {file_path}: {e}")
                # Include file if age check fails
                filtered_files.append(file_path)
        all_files = filtered_files

    if force_gzip:
        # Only return gzipped files
        gzip_files = [f for f in all_files if f.endswith('.gz')]
        print(f"Force gzip mode: found {len(gzip_files)} gzipped files out of {len(all_files)} total")
        return gzip_files

    # Create dictionary to group files by base name
    file_groups = {}

    for file_path in all_files:
        if file_path.endswith('.gz'):
            # Remove .gz to get base name
            base_name = file_path[:-3]
            if base_name not in file_groups:
                file_groups[base_name] = {'gz': None, 'normal': None}
            file_groups[base_name]['gz'] = file_path
        else:
            # This is the base name
            base_name = file_path
            if base_name not in file_groups:
                file_groups[base_name] = {'gz': None, 'normal': None}
            file_groups[base_name]['normal'] = file_path

    # Build final file list, preferring gzipped versions
    final_files = []
    for base_name, versions in file_groups.items():
        if versions['gz'] is not None:
            # Prefer gzipped version
            final_files.append(versions['gz'])
            if versions['normal'] is not None:
                print(f"Found both versions, preferring gzipped: {os.path.basename(versions['gz'])}")
        elif versions['normal'] is not None:
            # Use uncompressed version
            final_files.append(versions['normal'])

    final_files.sort()
    return final_files

def load_yaml_config(yaml_project_file, yaml_instrument_file):
    """
    Load configuration from YAML files, including north_angle.
    
    Args:
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        
    Returns:
        Dictionary of configuration parameters
    """
    config = {
        'north_angle': 55.62,  # Default fallback
        'data_version': '1.0.0',
        'location': 'cao'
    }
    
    # Load project YAML
    try:
        with open(yaml_project_file, 'r') as f:
            project_data = yaml.safe_load(f)
        
        print(f"Loading configuration from: {yaml_project_file}")
        
        # Navigate through your YAML structure: ncas_instruments -> ncas-mobile-ka-band-radar-1
        if 'ncas_instruments' in project_data:
            for instrument in project_data['ncas_instruments']:
                if 'ncas-mobile-ka-band-radar-1' in instrument:
                    radar_config = instrument['ncas-mobile-ka-band-radar-1']
                    
                    # Extract north_angle
                    if 'north_angle' in radar_config:
                        config['north_angle'] = float(radar_config['north_angle'])
                        print(f"Loaded north_angle from YAML: {config['north_angle']}°")
                    
                    # Extract other parameters if available
                    if 'processing_software' in radar_config and 'version' in radar_config['processing_software']:
                        config['data_version'] = radar_config['processing_software']['version']
                    
                    if 'platform' in radar_config and 'location' in radar_config['platform']:
                        config['location'] = radar_config['platform']['location'].lower()
                    
                    break
        
    except Exception as e:
        print(f"Warning: Could not load project YAML config: {e}")
        print(f"Using default north_angle: {config['north_angle']}°")
    
    return config

def get_campaign_info(campaign: str) -> Dict[str, Any]:
    """
    Get campaign-specific configuration parameters.
    
    Args:
        campaign: Campaign name
        
    Returns:
        Dictionary of campaign-specific parameters
    """
    campaign_configs = {
        'cobalt': {
            'revised_northangle': 55.62,
            'tracking_tag': 'AMOF_20231120125118'
        },
        'woest': {
            'revised_northangle': 302.15,  # Different north angle for WOEST
            'tracking_tag': 'AMOF_20220301000000'
        },
        'ccrest': {
            'revised_northangle': 55.7,
            'tracking_tag': 'AMOF_20230401000000'
        },
        'kasbex': {
            'revised_northangle': 55.62,
            'tracking_tag': 'AMOF_20250508133639'  # KASBEX-specific tracking tag
        }
    }
    
    return campaign_configs.get(campaign, {
        'revised_northangle': 55.62,  # Default
        'tracking_tag': f'AMOF_{campaign.upper()}'
    })

def group_files_by_range_gates(files, gzip_flag=False):
    """
    Group mmclx files by their range gate dimensions to avoid concatenation errors.
    
    Args:
        files (list): List of mmclx file paths
        gzip_flag (bool): Whether files are gzipped
        
    Returns:
        dict: Dictionary with range_gates as keys and lists of files as values
    """
    from kepler_utils import read_mira35_mmclx
    import os
    
    file_groups = {}
    
    for file_path in files:
        try:
            # Read just the first sweep to get dimensions
            radar_ds = read_mira35_mmclx(file_path, gzip_flag=gzip_flag, revised_northangle=55.62)
            n_gates = radar_ds.ngates
            
            if n_gates not in file_groups:
                file_groups[n_gates] = []
            file_groups[n_gates].append(file_path)
            
            print(f"File {os.path.basename(file_path)}: {n_gates} gates")
            
        except Exception as e:
            print(f"Warning: Could not read {file_path} for grouping: {e}")
            # Put problematic files in a separate group
            if 'error' not in file_groups:
                file_groups['error'] = []
            file_groups['error'].append(file_path)
    
    return file_groups

# Example usage functions
def process_campaign_day(
    campaign: str,
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    gzip_flag: bool = True,
    data_version: str = "1.0.0",
    single_sweep: bool = False,
    revised_northangle: float = None,
    no_vpt: bool = False,  # Add no_vpt argument
    max_age: Optional[float] = None,
    **kwargs
) -> None:
    """
    Process a single day of campaign data.
    
    Args:
        campaign: Campaign name ('woest', 'ccrest', 'cobalt', 'kasbex')
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        outpath: Output directory path
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        gzip_flag: Whether input files are gzip compressed
        data_version: Data version string
        single_sweep: If True, create separate files for each sweep
        revised_northangle: North angle override (if None, load from YAML)
        **kwargs: Additional campaign-specific arguments
    """
    
    # Load configuration from YAML - THIS IS THE KEY MISSING PIECE
    yaml_config = load_yaml_config(yaml_project_file, yaml_instrument_file)
    
    # Use revised_northangle from command line if provided, otherwise use YAML value
    if revised_northangle is None:
        revised_northangle = yaml_config['north_angle']
        print(f"Using north_angle from YAML: {revised_northangle}°")
    else:
        print(f"Using north_angle from command line override: {revised_northangle}°")
    
    # Get campaign-specific processor and parameters
    campaign_lower = campaign.lower()
    
    # Define available processors
    processors = {
        'woest': process_kepler_woest_day_step1,
        'ccrest': process_kepler_ccrest_day_step1,
        'cobalt': process_kepler_cobalt_day_step1,
        'kasbex': process_kepler_kasbex_day_step1
    }
    
    if campaign_lower not in processors:
        raise ValueError(f"Unknown campaign: {campaign}. Available: {list(processors.keys())}")
    
    # Get campaign-specific configuration
    campaign_info = get_campaign_info(campaign_lower)
    
    # Set up processing arguments
    processing_args = {
        'datestr': datestr,
        'inpath': inpath,
        'outpath': outpath,
        'yaml_project_file': yaml_project_file,
        'yaml_instrument_file': yaml_instrument_file,
        'gzip_flag': gzip_flag,
        'data_version': data_version,
        'single_sweep': single_sweep,
        'revised_northangle': revised_northangle,  # Use the north angle from YAML or override
        'tracking_tag': campaign_info.get('tracking_tag', f'AMOF_{campaign.upper()}'),
        'no_vpt': no_vpt,  # Pass no_vpt argument
        **kwargs
    }
    
    print(f"Processing {campaign.upper()} campaign data for {datestr}")
    print(f"Configuration:")
    print(f"  North angle: {revised_northangle}°")
    print(f"  Data version: {data_version}")
    print(f"  Single sweep: {single_sweep}")
    print(f"  Tracking tag: {processing_args['tracking_tag']}")
    
    # Call the appropriate processor
    processor = processors[campaign_lower]
    processor(**processing_args)
    
    print(f"Completed {campaign.upper()} processing for {datestr}")