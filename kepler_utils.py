#!/usr/bin/env python
# coding: utf-8

"""
kepler_utils_new.py

Module for processing mmclx radar files from Kepler (MIRA-35) radar

This module provides functions to:
- Read MIRA-35 radar data in mmclx format
- Convert to CF-Radial format
- Add NCAS metadata
- Process data for various campaigns (WOEST, CCREST, COBALT, KASBEX)

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
Last modified: 10-08-2025
Version: 1.1.0
"""

from typing import List, Dict, Tuple, Optional, Union, Any
import datetime
import cftime
import netCDF4 as nc4
import numpy as np
import getpass
import socket
import pyart
from pyart.config import FileMetadata, get_fillvalue, get_metadata
from pyart.core.radar import Radar
from pyart.io.common import _test_arguments, make_time_unit_str
import yaml
import os
import glob
import fnmatch
import re
import gzip
from io import StringIO
from pathlib import Path

module_version = "1.1.0"

def read_mira35_mmclx_hsrhi(mmclxfiles: List[str], **kwargs) -> Radar:
    """
    Read a set of single-sweep netCDF mmclx files from MIRA-35 radar recorded as part of HSRHI scan strategy.
    
    Args:
        mmclxfiles: List of paths to mmclx netCDF files to read data from
        **kwargs: Additional keyword arguments
        
    Returns:
        Radar object containing the merged data
        
    Note:
        This function is currently incomplete and needs implementation.
    """
    mmclxfiles.sort()
    nsweep = len(mmclxfiles)
    # TODO: Implement HSRHI reading functionality
    pass

def read_mira35_mmclx_vpt_multi(mmclxfiles: List[str], **kwargs) -> Radar:
    """
    Read a set of netCDF mmclx files from MIRA-35 radar recorded as separate vertical pointing sweeps.
    
    Args:
        mmclxfiles: List of paths to mmclx netCDF files to read data from
        **kwargs: Additional keyword arguments
        
    Returns:
        Radar object containing the merged VPT data
        
    Note:
        This function is currently incomplete and needs implementation.
    """
    mmclxfiles.sort()
    nsweep = len(mmclxfiles)
    # TODO: Implement VPT multi-sweep reading functionality
    pass

def read_mira35_mmclx(
    filename: str,
    gzip_flag: bool = False,
    revised_northangle: float = 55.9,
    logp_dir: str = None,
    **kwargs
) -> Radar:
    """
    Read a netCDF mmclx file from MIRA-35 radar and convert to PyART Radar object.
    
    This function reads a single mmclx file (either compressed or uncompressed) and converts
    it to a PyART Radar object with appropriate metadata and field mappings.
    
    Args:
        filename: Path to the mmclx netCDF file
        gzip_flag: If True, treat file as gzip-compressed
        revised_northangle: North angle correction in degrees (default: 55.7)
        logp_dir: Path to directory of axis-log files (*.axis.gz). When provided,
            used to reconstruct azimuth/elevation for scans where the mmclx angle
            variables contain fill values (e.g. axis controller offline).
        **kwargs: Additional keyword arguments passed to PyART
        
    Returns:
        PyART Radar object containing the radar data
        
    Raises:
        IOError: If file cannot be opened
        ValueError: If required variables are missing from file
        
    Example:
        >>> radar = read_mira35_mmclx('data.mmclx', gzip_flag=False, revised_northangle=56.0)
        >>> print(f"Number of gates: {radar.ngates}")
    """
    # Validate arguments
    _test_arguments(kwargs)
    
    # Initialize metadata handler
    filemetadata = FileMetadata('mmclx')
    
    # Open file (compressed or uncompressed)
    print(f"Reading file: {filename}, gzip_flag={gzip_flag}")
    
    if gzip_flag:
        with gzip.open(filename) as gz:
            ncobj = nc4.Dataset('dummy', mode='r', memory=gz.read())
    else:
        ncobj = nc4.Dataset(filename)
    
    try:
        # Get dimensions
        nrays = len(ncobj.dimensions["time"])
        ngates = len(ncobj.dimensions["range"])
        nsweeps = 1  # Single sweep files only
        
        print(f"File dimensions - nrays: {nrays}, ngates: {ngates}, nsweeps: {nsweeps}")
        
        # Extract location information from global attributes
        latitude, longitude, altitude = _extract_location_info(ncobj, filemetadata)
        
        # Extract metadata
        metadata = _extract_metadata(ncobj, filemetadata)
        
        # Extract time information
        time = _extract_time_info(ncobj, filemetadata)
        
        # Extract range information
        _range = _extract_range_info(ncobj, filemetadata)
        
        # Determine scan type and fixed angle
        scan_name, sweep_mode, fixed_angle = _determine_scan_type(filename, ncobj, revised_northangle, filemetadata)
        
        # Extract angle information
        angle_info = _extract_angle_info(
            ncobj, scan_name, time, revised_northangle, filemetadata, logp_dir=logp_dir
        )
        azimuth = angle_info['azimuth']
        elevation = angle_info['elevation']
        scan_rate = angle_info.get('scan_rate')
        antenna_transition = angle_info.get('antenna_transition')
        target_scan_rate = angle_info.get('target_scan_rate')
        elevation_scan_rate = angle_info.get('elevation_scan_rate')
        azimuth_scan_rate = angle_info.get('azimuth_scan_rate')

        # Append provenance notes to global metadata comment when angles or scan
        # rates were reconstructed from logp axis files or corrected from fill values.
        existing_comment = metadata.get('comment', '') or ''
        extra_comments = []
        if angle_info.get('logp_reconstructed'):
            extra_comments.append(
                "Azimuth, elevation, and scan_rate were reconstructed from the MIRA "
                "axis-controller log (logp *.axis.gz files) because the mmclx angle "
                "variables (azi, elv, elvv, aziv) contained fill values indicating "
                "the axis controller was not communicating with the radar software."
            )
        if angle_info.get('vpt_fill_corrected'):
            extra_comments.append(
                "VPT elevation fill values (< -900\u00b0) were replaced with 90.0\u00b0 and "
                "azimuth fill values were replaced with 0.0\u00b0 (undefined for vertical "
                "pointing). The axis controller was not communicating with the radar "
                "software; logp axis files confirm the antenna was at vertical pointing."
            )
        if extra_comments:
            joined = ' '.join(extra_comments)
            metadata['comment'] = (existing_comment.strip() + ' ' + joined).strip() if existing_comment.strip() else joined
        
        # Create sweep information
        sweep_info = _create_sweep_info(nrays, filemetadata)
        
        # Extract radar fields
        fields = _extract_radar_fields(ncobj, scan_name, time, filemetadata)
        
        # Extract instrument parameters
        instrument_parameters = _extract_instrument_parameters(ncobj, nrays, filemetadata)
        
    finally:
        ncobj.close()
        if gzip_flag:
            gz.close()
    
    # Create and return Radar object
    radar = Radar(
        time, _range, fields, metadata, scan_name,
        latitude, longitude, altitude,
        sweep_info['sweep_number'], sweep_mode, fixed_angle,
        sweep_info['sweep_start_ray_index'], sweep_info['sweep_end_ray_index'],
        azimuth, elevation,
        target_scan_rate=target_scan_rate,
        scan_rate=scan_rate,
        antenna_transition=antenna_transition,
        instrument_parameters=instrument_parameters,
        radar_calibration={}  # Empty for now
    )
    
    # Add elevation and azimuth scan rates for MAN scans (as coordinate variables)
    if elevation_scan_rate is not None:
        radar.elevation_scan_rate = elevation_scan_rate
    if azimuth_scan_rate is not None:
        radar.azimuth_scan_rate = azimuth_scan_rate
    
    print("Successfully created Radar object")
    return radar

def _extract_location_info(ncobj: nc4.Dataset, filemetadata: FileMetadata) -> Tuple[Dict, Dict, Dict]:
    """
    Extract latitude, longitude, and altitude from netCDF object.
    
    Args:
        ncobj: Open netCDF4 Dataset
        filemetadata: PyART FileMetadata object
        
    Returns:
        Tuple of (latitude_dict, longitude_dict, altitude_dict)
    """
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude') 
    altitude = filemetadata('altitude')
    
    # Parse latitude - use float64 for geolocation precision
    lat_str = ncobj.getncattr('Latitude')
    z = StringIO(lat_str)
    lat_data = np.genfromtxt(z, dtype=None, names=['lat', 'direction'])
    
    if lat_data['direction'] == b'S' and lat_data['lat'] > 0:
        latitude['data'] = np.array([-lat_data['lat']], dtype='f8')
    else:
        latitude['data'] = np.array([lat_data['lat']], dtype='f8')
    
    # Parse longitude - use float64 for geolocation precision
    lon_str = ncobj.getncattr('Longitude')
    z = StringIO(lon_str)
    lon_data = np.genfromtxt(z, dtype=None, names=['lon', 'direction'])
    
    if lon_data['direction'] == b'W' and lon_data['lon'] > 0:
        longitude['data'] = np.array([-lon_data['lon']], dtype='f8')
    else:
        longitude['data'] = np.array([lon_data['lon']], dtype='f8')
    
    # Parse altitude - use float64 for geolocation precision
    alt_str = ncobj.getncattr('Altitude')
    z = StringIO(alt_str)
    alt_data = np.genfromtxt(z, dtype=None, names=['alt', 'units'])
    altitude['data'] = np.array([alt_data['alt']], dtype='f8')
    
    return latitude, longitude, altitude

def _extract_metadata(ncobj: nc4.Dataset, filemetadata: FileMetadata) -> Dict:
    """
    Extract global metadata from netCDF object.
    
    Args:
        ncobj: Open netCDF4 Dataset
        filemetadata: PyART FileMetadata object
        
    Returns:
        Dictionary containing metadata
    """
    metadata = filemetadata('metadata')
    
    # Extract standard attributes
    for attr in ['institution', 'title', 'used_algorithms']:
        if attr in ncobj.ncattrs():
            metadata[attr] = ncobj.getncattr(attr)
    
    metadata['instrument_name'] = 'ncas-mobile-ka-band-radar-1'
    
    return metadata

def _extract_time_info(ncobj: nc4.Dataset, filemetadata: FileMetadata) -> Dict:
    """
    Extract and process time information from netCDF object.
    
    Args:
        ncobj: Open netCDF4 Dataset
        filemetadata: PyART FileMetadata object
        
    Returns:
        Dictionary containing time information
    """
    time = filemetadata('time')
    
    # Convert time with microsecond precision
    dtime = cftime.num2pydate(ncobj.variables['time'][:], 'seconds since 1970-01-01 00:00:00')
    
    # Add microsecond information
    for idx, dt in np.ndenumerate(dtime):
        dtime[idx] = dt.replace(microsecond=ncobj.variables['microsec'][idx])
    
    # Set time reference to start of day
    base_time = dtime[0].replace(hour=0, minute=0, second=0, microsecond=0)
    
    time['units'] = make_time_unit_str(base_time)
    time['data'] = cftime.date2num(dtime, time['units'])
    
    print(f"Time data range: {time['data'][0]:.3f} to {time['data'][-1]:.3f} {time['units']}")
    
    return time

def _extract_range_info(ncobj: nc4.Dataset, filemetadata: FileMetadata) -> Dict:
    """
    Extract range gate information from netCDF object.
    
    Args:
        ncobj: Open netCDF4 Dataset  
        filemetadata: PyART FileMetadata object
        
    Returns:
        Dictionary containing range information
    """
    _range = filemetadata('range')
    _range['data'] = ncobj.variables['range'][:]
    _range['units'] = 'metres'
    _range.pop('standard_name', None)
    _range['proposed_standard_name'] = "projection_range_coordinate"
    _range['long_name'] = "distance to centre of each range gate"
    
    return _range

def _determine_scan_type(filename: str, ncobj: nc4.Dataset, revised_northangle: float, filemetadata: FileMetadata) -> Tuple[str, Dict, Dict]:
    """
    Determine scan type and calculate fixed angle from filename and data.
    
    Args:
        filename: Path to the data file
        ncobj: Open netCDF4 Dataset
        revised_northangle: North angle correction in degrees
        filemetadata: PyART FileMetadata object
        
    Returns:
        Tuple of (scan_name, sweep_mode_dict, fixed_angle_dict)
    """
    sweep_modes = {
        'ppi': 'ppi',
        'rhi': 'rhi', 
        'vert': 'vertical_pointing',
        'man': 'manual_rhi'
    }
    
    sweep_mode = filemetadata("sweep_mode")
    fixed_angle = filemetadata("fixed_angle")
    
    # Determine scan type from filename
    scan_name = None
    filename_lower = filename.lower()
    
    for key, value in sweep_modes.items():
        if key in filename_lower:
            scan_name = value
            sweep_mode["data"] = np.array([value])
            break
    
    if scan_name is None:
        # No keyword matched in the filename (e.g. plain YYYYMMDD_HHMMSS.mmclx[.gz] files).
        # Fall back to elevation: if mean elevation ≥ 80° treat as VPT.
        if 'elv' in ncobj.variables:
            mean_elv = float(ncobj.variables['elv'][:].mean())
            if mean_elv >= 80.0:
                print(
                    f"No scan-type keyword in filename; mean elevation={mean_elv:.1f}° "
                    f"→ treating as vertical_pointing"
                )
                scan_name = 'vertical_pointing'
                sweep_mode["data"] = np.array(['vertical_pointing'])
                fixed_angle["data"] = np.array([90.0], dtype='f')
                print(f"Detected scan type: {scan_name}, fixed angle: 90.0")
                return scan_name, sweep_mode, fixed_angle
            else:
                print(
                    f"No scan-type keyword in filename; mean elevation={mean_elv:.1f}° "
                    f"→ scan type unknown (manual classification required)"
                )
        sweep_mode["data"] = np.array([None])
        fixed_angle["data"] = np.array([None])
        return scan_name, sweep_mode, fixed_angle
    
    # Calculate fixed angle based on scan type
    nrays = len(ncobj.dimensions['time'])
    i4fixed_angle = min(4, nrays - 1)  # Use 5th ray or last ray if fewer than 5
    
    ncvars = ncobj.variables
    
    if scan_name in ['rhi', 'manual_rhi']:
        # For RHI, fixed angle is the azimuth
        fixed_angle_value = np.round((ncvars['azi'][i4fixed_angle] + revised_northangle) % 360, 2)
        fixed_angle["data"] = np.array([fixed_angle_value], dtype='f')
    elif scan_name == 'vertical_pointing':
        # For VPT, fixed angle is always 90° regardless of the stored elevation value
        # (elevation may contain fill values such as -1000 when the mmclx was missing data)
        fixed_angle["data"] = np.array([90.0], dtype='f')
    elif scan_name == 'ppi':
        # For PPI, fixed angle is the elevation
        fixed_angle_value = np.round(ncvars['elv'][i4fixed_angle], 2)
        fixed_angle["data"] = np.array([fixed_angle_value], dtype='f')
    
    print(f"Detected scan type: {scan_name}, fixed angle: {fixed_angle['data'][0]}")
    
    return scan_name, sweep_mode, fixed_angle


def _load_logp_angles(
    t_mmclx: np.ndarray,
    logp_dir: str
) -> tuple:
    """
    Load axis-controller angle data from logp files and interpolate onto mmclx timestamps.

    Logp files are named YYYYMMDDHHMI.axis.gz and record antenna azimuth and elevation
    at ~2 Hz. This function finds all logp files covering the mmclx time window, loads
    them, and linearly interpolates azimuth and elevation onto each mmclx ray timestamp.

    Args:
        t_mmclx: 1-D array of mmclx ray timestamps (Unix seconds, float/int)
        logp_dir: Path to the directory containing *.axis.gz files

    Returns:
        Tuple (azi_interp, elv_interp) of float32 arrays the same length as t_mmclx,
        or (None, None) if no usable logp data is found.
    """
    import datetime as _dt

    t0 = float(t_mmclx[0])
    t1 = float(t_mmclx[-1])

    # Collect candidate logp files: check the date of t0 and also the previous and
    # next days to handle midnight crossings.
    candidate_dates = set()
    for ts in (t0 - 86400, t0, t1):
        d = _dt.datetime.utcfromtimestamp(ts)
        candidate_dates.add(d.strftime('%Y%m%d'))

    logp_ts_all = []
    logp_azi_all = []
    logp_elv_all = []

    for datestr in sorted(candidate_dates):
        pattern = os.path.join(logp_dir, f'{datestr}*.axis.gz')
        for fpath in sorted(glob.glob(pattern)):
            with gzip.open(fpath, 'rt') as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    try:
                        ts = float(parts[0])
                    except ValueError:
                        continue
                    # Only keep rows within a generous window around the mmclx sweep
                    if ts < t0 - 60 or ts > t1 + 60:
                        continue
                    if parts[1] == 'nan' or parts[2] == 'nan':
                        continue
                    try:
                        logp_ts_all.append(ts)
                        logp_azi_all.append(float(parts[1]))
                        logp_elv_all.append(float(parts[2]))
                    except ValueError:
                        continue

    if len(logp_ts_all) < 2:
        return None, None

    logp_ts = np.array(logp_ts_all, dtype=np.float64)
    logp_azi = np.array(logp_azi_all, dtype=np.float64)
    logp_elv = np.array(logp_elv_all, dtype=np.float64)

    # Sort in case files were not chronological
    sort_idx = np.argsort(logp_ts)
    logp_ts = logp_ts[sort_idx]
    logp_azi = logp_azi[sort_idx]
    logp_elv = logp_elv[sort_idx]

    t_query = np.array(t_mmclx, dtype=np.float64)
    azi_interp = np.interp(t_query, logp_ts, logp_azi).astype(np.float32)
    elv_interp = np.interp(t_query, logp_ts, logp_elv).astype(np.float32)

    return azi_interp, elv_interp


def _extract_angle_info(
    ncobj: nc4.Dataset,
    scan_name: str,
    time: Dict,
    revised_northangle: float,
    filemetadata: FileMetadata,
    logp_dir: str = None
) -> Dict:
    """
    Extract azimuth, elevation and scan rate information.
    
    Args:
        ncobj: Open netCDF4 Dataset
        scan_name: Type of scan (ppi, rhi, manual_rhi, etc.)
        time: Time information dictionary
        revised_northangle: North angle correction in degrees
        filemetadata: PyART FileMetadata object
        
    Returns:
        Dictionary containing azimuth, elevation, scan_rate, antenna_transition, 
        target_scan_rate, and for MAN scans: elevation_scan_rate, azimuth_scan_rate
    """
    ncvars = ncobj.variables
    
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    
    # Basic angle assignments - explicitly cast to float32 for NCAS compliance
    raw_azi = ncvars['azi'][:]
    azimuth['data'] = ((raw_azi + revised_northangle) % 360).astype(np.float32)
    elevation['data'] = ncvars['elv'][:].astype(np.float32)

    # Detect MIRA axis-controller fill values (known values: -1000, -1001; use
    # threshold < -900 to be robust against any similar fill).
    angles_are_fill = (elevation['data'] < -900.0).all()

    # Calculate ray durations for scan rate correction
    ray_duration = np.diff(time['data'])
    
    result = {
        'azimuth': azimuth,
        'elevation': elevation
    }
    
    if scan_name in ['ppi', 'rhi', 'manual_rhi']:
        # If the mmclx angle variables contain fill values and a logp directory is
        # available, reconstruct azimuth and elevation by interpolating the axis-log.
        logp_reconstructed = False
        if angles_are_fill and logp_dir is not None and os.path.isdir(logp_dir):
            print(f"  [logp] Angle fill detected for {scan_name.upper()} — attempting logp reconstruction")
            # Use the raw mmclx timestamps (Unix seconds) not time['data'] which is
            # in "seconds since" PyART format and would not match the logp Unix times.
            azi_logp, elv_logp = _load_logp_angles(ncvars['time'][:], logp_dir)
            if azi_logp is not None:
                # Apply north-angle correction to logp azimuth (same reference frame as mmclx)
                azimuth['data'] = ((azi_logp + revised_northangle) % 360).astype(np.float32)
                elevation['data'] = elv_logp
                logp_reconstructed = True
                _logp_comment = (
                    "Values reconstructed by interpolating the MIRA axis-controller log "
                    "(logp *.axis.gz) onto mmclx ray timestamps. "
                    "The mmclx elvv/aziv scan-rate variables contained fill values "
                    "(axis controller offline) and were not used."
                )
                azimuth['comment'] = _logp_comment
                elevation['comment'] = _logp_comment
                print(f"  [logp] Reconstructed: az range [{azimuth['data'].min():.2f}, {azimuth['data'].max():.2f}], "
                      f"el range [{elevation['data'].min():.2f}, {elevation['data'].max():.2f}]")
            else:
                print(f"  [logp] No usable logp data found — angles remain as fill values")

        # For scanning modes, calculate scan rates and antenna transitions
        antenna_transition = filemetadata("antenna_transition") 
        target_scan_rate = filemetadata("target_scan_rate")
        
        if scan_name == 'manual_rhi':
            # MAN scans: track both elevation and azimuth rates
            elevation_scan_rate = filemetadata("elevation_scan_rate")
            azimuth_scan_rate = filemetadata("azimuth_scan_rate")
            
            if logp_reconstructed:
                # elvv/aziv are fill values — derive rates from logp-reconstructed angles
                elevation_scan_rate['data'] = np.gradient(
                    elevation['data'].astype(np.float64), time['data']).astype(np.float32)
                # Azimuth is nominally fixed for MAN RHI; any variation is noise
                azimuth_scan_rate['data'] = np.zeros(len(elevation['data']), dtype=np.float32)
            else:
                elevation_scan_rate['data'] = ncvars['elvv'][:]
                azimuth_scan_rate['data'] = ncvars['aziv'][:]
            elevation_scan_rate['units'] = 'degrees_per_second'
            elevation_scan_rate['long_name'] = 'antenna elevation angle scan rate'
            elevation_scan_rate['standard_name'] = 'platform_elevation_scan_rate'
            
            azimuth_scan_rate['units'] = 'degrees_per_second'
            azimuth_scan_rate['long_name'] = 'antenna azimuth angle scan rate'
            azimuth_scan_rate['standard_name'] = 'platform_azimuth_scan_rate'
            
            # For antenna_transition detection, use combined rate magnitude
            combined_rate = np.sqrt(elevation_scan_rate['data']**2 + azimuth_scan_rate['data']**2)
            antenna_transition['data'] = np.where(combined_rate < 0.01, 1, 0).astype('int8')
            antenna_transition['long_name'] = "antenna is in transition between sweeps"
            antenna_transition['comment'] = "1 if antenna is in transition, 0 otherwise"
            
            # For target scan rate, use mean elevation rate from scanning periods
            scanning_indices = np.where(antenna_transition['data'] == 0)[0]
            if len(scanning_indices) > 0:
                target_rate = np.mean(np.abs(elevation_scan_rate['data'][scanning_indices]))
                target_scan_rate['data'] = np.round(target_rate, 2)
            else:
                target_scan_rate['data'] = np.array([1.0], dtype="f4")
            
            target_scan_rate['long_name'] = 'target elevation scan rate for sweep'
            
            # For backward compatibility, set scan_rate to elevation rate
            scan_rate = filemetadata("scan_rate")
            scan_rate['data'] = elevation_scan_rate['data']
            scan_rate['units'] = 'degrees_per_second'
            scan_rate['long_name'] = 'antenna angle scan rate (elevation for MAN scans)'
            if logp_reconstructed:
                _sr_comment = (
                    "Derived from np.gradient of logp-reconstructed elevation angles "
                    "over mmclx ray timestamps; mmclx elvv/aziv contained fill values."
                )
                scan_rate['comment'] = _sr_comment
                elevation_scan_rate['comment'] = _sr_comment
                azimuth_scan_rate['comment'] = (
                    "Set to zero: azimuth is nominally fixed for MAN RHI scans; "
                    "mmclx aziv contained fill values."
                )
            
            result.update({
                'scan_rate': scan_rate,
                'antenna_transition': antenna_transition,
                'target_scan_rate': target_scan_rate,
                'elevation_scan_rate': elevation_scan_rate,
                'azimuth_scan_rate': azimuth_scan_rate
            })
            
        elif scan_name == 'ppi':
            # PPI: azimuth scanning only
            scan_rate = filemetadata("scan_rate")
            if logp_reconstructed:
                # aziv is a fill value — derive azimuth scan rate from logp-reconstructed angles
                scan_rate['data'] = np.gradient(
                    azimuth['data'].astype(np.float64), time['data']).astype(np.float32)
                scan_rate['comment'] = (
                    "Derived from np.gradient of logp-reconstructed azimuth angles "
                    "over mmclx ray timestamps; mmclx aziv contained fill values."
                )
            else:
                scan_rate['data'] = ncvars['aziv'][:]
            scan_rate['units'] = 'degrees_per_second'
            scan_rate['long_name'] = 'antenna angle scan rate'
            
            # Identify antenna transitions using multiple criteria:
            # 1. Very low scan rates (< 0.01 deg/s)
            low_rate_flags = np.abs(scan_rate['data']) < 0.01
            
            # 2. Non-monotonic azimuth (detect direction reversals)
            azimuth_data = azimuth['data']
            azimuth_diff = np.diff(azimuth_data)
            # Handle 360° wrap-around
            azimuth_diff = np.where(azimuth_diff > 180, azimuth_diff - 360, azimuth_diff)
            azimuth_diff = np.where(azimuth_diff < -180, azimuth_diff + 360, azimuth_diff)
            
            # Determine predominant scan direction
            median_diff = np.median(azimuth_diff)
            if median_diff > 0:
                # Scanning clockwise - flag counter-clockwise rays as transitions
                non_monotonic = np.concatenate([[False], azimuth_diff < -1.0])  # -1° tolerance
            else:
                # Scanning counter-clockwise - flag clockwise rays as transitions  
                non_monotonic = np.concatenate([[False], azimuth_diff > 1.0])   # +1° tolerance
            
            # Combine criteria
            antenna_transition['data'] = np.where(low_rate_flags | non_monotonic, 1, 0).astype('int8')
            antenna_transition['long_name'] = "antenna is in transition between sweeps"
            antenna_transition['comment'] = "1 if antenna is in transition (low rate or non-monotonic), 0 otherwise"
            
            # Calculate target scan rate from non-transitioning periods
            scanning_indices = np.where(antenna_transition['data'] == 0)[0]
            if len(scanning_indices) > 0:
                target_rate = np.mean(scan_rate['data'][scanning_indices])
                target_scan_rate['data'] = np.round(target_rate, 2)
            else:
                target_scan_rate['data'] = np.array([4.0], dtype="f4")
                
            target_scan_rate['long_name'] = 'target scan rate for sweep'
            
            result.update({
                'scan_rate': scan_rate,
                'antenna_transition': antenna_transition,
                'target_scan_rate': target_scan_rate
            })
            
        else:  # rhi
            # RHI: elevation scanning only
            scan_rate = filemetadata("scan_rate")
            if logp_reconstructed:
                # elvv is a fill value — derive elevation scan rate from logp-reconstructed angles
                scan_rate['data'] = np.gradient(
                    elevation['data'].astype(np.float64), time['data']).astype(np.float32)
                scan_rate['comment'] = (
                    "Derived from np.gradient of logp-reconstructed elevation angles "
                    "over mmclx ray timestamps; mmclx elvv contained fill values."
                )
            else:
                scan_rate['data'] = ncvars['elvv'][:]
            scan_rate['units'] = 'degrees_per_second'
            scan_rate['long_name'] = 'antenna angle scan rate'
            
            # Identify antenna transitions (very low scan rates)
            antenna_transition['data'] = np.where(np.abs(scan_rate['data']) < 0.01, 1, 0).astype('int8')
            antenna_transition['long_name'] = "antenna is in transition between sweeps"
            antenna_transition['comment'] = "1 if antenna is in transition, 0 otherwise"
            
            # Calculate target scan rate from non-transitioning periods
            scanning_indices = np.where(antenna_transition['data'] == 0)[0]
            if len(scanning_indices) > 0:
                target_rate = np.mean(scan_rate['data'][scanning_indices])
                target_scan_rate['data'] = np.round(target_rate, 2)
            else:
                target_scan_rate['data'] = np.array([4.0], dtype="f4")
                
            target_scan_rate['long_name'] = 'target scan rate for sweep'
            
            result.update({
                'scan_rate': scan_rate,
                'antenna_transition': antenna_transition,
                'target_scan_rate': target_scan_rate
            })
        
        # Propagate the reconstruction flag so callers can annotate global metadata.
        result['logp_reconstructed'] = logp_reconstructed

        # Apply timing corrections to angles (pass antenna_transition for RHI and MAN)
        # Skip if logp reconstruction was used — the scan rates (elvv/aziv) are also
        # fill values when the axis controller is offline, so timing corrections would
        # corrupt the reconstructed angles.
        if not logp_reconstructed:
            _apply_timing_corrections(azimuth, elevation, ncvars, ray_duration,
                                      revised_northangle, antenna_transition['data'])

    elif scan_name in ['vpt', 'vertical_pointing']:
        # VPT: elevation should be ~90°. MIRA uses fill values around -1000 when
        # the axis controller data is not available (observed values: -1000 when
        # the controller is completely offline, -1001 when it is logging but not
        # communicating with the radar software). Both are caught by the < -900
        # threshold. Logp axis files confirm these scans were genuinely at 90°.
        fill_mask = elevation['data'] < -900.0
        if fill_mask.any():
            elevation['data'] = np.where(fill_mask, np.float32(90.0), elevation['data'])
            _vpt_comment = (
                "Fill values (< -900°) replaced with 90.0°: the MIRA axis controller "
                "was not communicating with the radar software but logp axis files "
                "confirm the antenna was at vertical pointing."
            )
            elevation['comment'] = _vpt_comment
            # Azimuth is undefined for VPT; the fill value produces a spurious
            # wrapped angle — set to 0.0 (north) for affected rays.
            azi_fill_mask = raw_azi < -900.0
            if azi_fill_mask.any():
                azimuth['data'] = np.where(azi_fill_mask, np.float32(0.0), azimuth['data'])
                azimuth['comment'] = (
                    "Fill values (< -900°) replaced with 0.0° (north): azimuth is "
                    "undefined for vertical pointing; the MIRA axis controller was "
                    "not communicating with the radar software."
                )
            result['vpt_fill_corrected'] = True
        # Mark any ray below 85° as antenna_transition=1
        # (radar was not yet pointing vertically).
        antenna_transition = filemetadata("antenna_transition")
        antenna_transition['data'] = np.where(
            elevation['data'] < 85.0, 1, 0
        ).astype('int8')
        antenna_transition['long_name'] = "antenna is in transition between sweeps"
        antenna_transition['comment'] = (
            "1 if elevation < 85° (antenna not yet at vertical pointing position), 0 otherwise"
        )
        result['antenna_transition'] = antenna_transition

    # Normalise azimuth/elevation metadata for all scan types:
    # PyART's FileMetadata('mmclx') seeds azimuth with standard_name='beam_azimuth_angle'
    # and elevation with standard_name='beam_elevation_angle'.  These are not the NCAS
    # names and must not appear in the output file.  Overwrite with proposed_standard_name
    # here so the PyART write step never emits a stray standard_name.
    azimuth.pop('standard_name', None)
    azimuth['proposed_standard_name'] = 'ray_azimuth_angle'
    azimuth['units'] = 'degrees'
    azimuth['long_name'] = 'azimuth_angle_from_true_north'

    elevation.pop('standard_name', None)
    elevation['proposed_standard_name'] = 'ray_elevation_angle'
    elevation['units'] = 'degrees'
    elevation['long_name'] = 'elevation_angle_from_horizontal_plane'

    return result

def _apply_timing_corrections(
    azimuth: Dict, 
    elevation: Dict, 
    ncvars: Dict, 
    ray_duration: np.ndarray,
    revised_northangle: float,
    antenna_transition: np.ndarray = None
) -> None:
    """
    Apply timing corrections to azimuth and elevation angles.
    
    This corrects for the fact that angles are recorded at the start of each ray,
    but we want angles at the center of each ray.
    
    Args:
        azimuth: Azimuth data dictionary (modified in place)
        elevation: Elevation data dictionary (modified in place) 
        ncvars: NetCDF variables
        ray_duration: Array of ray durations
        revised_northangle: North angle correction
        antenna_transition: Optional array of antenna transition flags (1=transition, 0=scanning)
    """
    # Calculate target ray duration using only non-transition rays if available
    if antenna_transition is not None:
        scanning_rays = np.where(antenna_transition == 0)[0]
        if len(scanning_rays) > 1:
            ray_duration_for_calc = ray_duration[scanning_rays[:-1]]  # Exclude last (no duration after it)
            target_ray_duration = np.round(np.mean(ray_duration_for_calc), 3)
            ok_duration = np.where(ray_duration_for_calc / target_ray_duration < 1.5)[0]
            if len(ok_duration) > 0:
                target_ray_duration = np.round(np.mean(ray_duration_for_calc[ok_duration]), 3)
        else:
            # Fallback if no scanning rays found
            target_ray_duration = np.round(np.mean(ray_duration), 3)
            ok_duration = np.where(ray_duration / target_ray_duration < 1.5)[0]
            if len(ok_duration) > 0:
                target_ray_duration = np.round(np.mean(ray_duration[ok_duration]), 3)
    else:
        # Original logic for scans without antenna_transition array
        target_ray_duration = np.round(np.mean(ray_duration), 3)
        ok_duration = np.where(ray_duration / target_ray_duration < 1.5)[0]
        
        if len(ok_duration) > 0:
            target_ray_duration = np.round(np.mean(ray_duration[ok_duration]), 3)
    
    # Prepend target duration for first ray
    ray_duration_extended = np.insert(ray_duration, 0, target_ray_duration)
    
    # Find rays with anomalously long durations
    long_duration = np.where(ray_duration_extended / target_ray_duration > 1.5)[0]
    
    print(f'Long-duration ray indices: {long_duration}')
    
    # Apply corrections for azimuth
    azimuth['data'] -= 0.5 * ray_duration_extended * ncvars['aziv'][:]
    azimuth['units'] = "degrees"
    azimuth['proposed_standard_name'] = "ray_azimuth_angle"
    azimuth['long_name'] = "azimuth_angle_from_true_north"
    
    # Special handling for long-duration rays and antenna transitions
    # For antenna_transition rays, don't apply timing correction (just use raw angles)
    special_rays = long_duration if antenna_transition is None else np.union1d(long_duration, np.where(antenna_transition == 1)[0])
    
    for idx in special_rays:
        if idx < len(azimuth['data']):
            azimuth['data'][idx] = (ncvars['azi'][idx] + revised_northangle) % 360
            # Only apply timing correction for long-duration rays (not for antenna_transition)
            if antenna_transition is None or antenna_transition[idx] == 0:
                azimuth['data'][idx] -= 0.5 * target_ray_duration * ncvars['aziv'][idx]
    
    # Apply corrections for elevation
    elevation['data'] -= 0.5 * ray_duration_extended * ncvars['elvv'][:]
    elevation['units'] = "degrees"
    elevation['proposed_standard_name'] = "ray_elevation_angle"
    elevation['long_name'] = "elevation_angle_from_horizontal_plane"
    
    # Special handling for long-duration rays and antenna transitions
    for idx in special_rays:
        if idx < len(elevation['data']):
            elevation['data'][idx] = ncvars['elv'][idx]
            # Only apply timing correction for long-duration rays (not for antenna_transition)
            if antenna_transition is None or antenna_transition[idx] == 0:
                elevation['data'][idx] -= 0.5 * target_ray_duration * ncvars['elvv'][idx]

def _create_sweep_info(nrays: int, filemetadata: FileMetadata) -> Dict:
    """
    Create sweep indexing information.
    
    Args:
        nrays: Number of rays in the sweep
        filemetadata: PyART FileMetadata object
        
    Returns:
        Dictionary containing sweep indexing information
    """
    sweep_info = {}
    
    sweep_info['sweep_start_ray_index'] = filemetadata("sweep_start_ray_index")
    sweep_info['sweep_end_ray_index'] = filemetadata("sweep_end_ray_index")
    sweep_info['sweep_number'] = filemetadata("sweep_number")
    
    sweep_info['sweep_start_ray_index']["data"] = np.array([0], dtype="int32")
    sweep_info['sweep_end_ray_index']["data"] = np.array([nrays - 1], dtype="int32") 
    sweep_info['sweep_number']["data"] = np.array([0], dtype="int32")
    
    return sweep_info

def _extract_radar_fields(ncobj: nc4.Dataset, scan_name: str, time: Dict, filemetadata: FileMetadata) -> Dict:
    """
    Extract radar fields (DBZ, VEL, WIDTH, etc.) from netCDF object.
    
    Args:
        ncobj: Open netCDF4 Dataset
        scan_name: Type of scan
        time: Time information for identifying bad rays
        filemetadata: PyART FileMetadata object
        
    Returns:
        Dictionary of radar fields
    """
    ncvars = ncobj.variables
    fields = {}
    
    # Field mappings from mmclx to CF-Radial names
    fields_keymap = {
        "Zg": "DBZ",      # Equivalent reflectivity factor
        "VELg": "VEL",    # Radial velocity  
        "RMSg": "WIDTH",  # Spectrum width
        "LDRg": "LDR",    # Linear depolarization ratio
        "SNRg": "SNR",    # Signal-to-noise ratio
        "RHO": "RHOHX",   # Correlation coefficient
        "DPS": "PHIHX",   # Differential phase
    }
    
    # Identify problematic rays for scanning modes
    long_duration_rays = []
    if scan_name in ['ppi', 'rhi']:
        ray_duration = np.diff(time['data'])
        target_duration = np.mean(ray_duration)
        long_duration_rays = np.where(ray_duration / target_duration > 1.5)[0]
    
    # Process each available field
    for mmclx_name, cfradial_name in fields_keymap.items():
        if mmclx_name in ncvars:
            print(f"Processing field: {mmclx_name} -> {cfradial_name}")
            
            field_dict = filemetadata(cfradial_name)
            field_dict['_FillValue'] = get_fillvalue()
            
            # Handle different field types
            if mmclx_name in ["Zg", "LDRg", "SNRg"]:
                # Convert linear to log scale.  Use np.float32(10) to avoid
                # the Python float (float64) silently upcasting the result.
                field_dict['data'] = (np.float32(10) * np.log10(ncvars[mmclx_name][:])).astype(np.float32)
                field_dict['units'] = 'dBZ' if mmclx_name == "Zg" else 'dB'
            else:
                # Keep in original units
                field_dict['data'] = ncvars[mmclx_name][:]

                if mmclx_name in ["VELg", "RMSg"]:
                    field_dict['units'] = 'm s-1'
                elif mmclx_name == "RHO":
                    field_dict['units'] = ''
                elif mmclx_name == "DPS":
                    field_dict['units'] = 'degrees'
            
            # Handle missing/invalid values
            invalid_mask = np.isnan(field_dict['data'])
            field_dict['data'][invalid_mask] = field_dict['_FillValue']
            
            # Mask problematic rays in scanning modes  
            if scan_name in ['ppi', 'rhi'] and len(long_duration_rays) > 0:
                field_dict['data'][long_duration_rays, :] = field_dict['_FillValue']
                field_dict['data'][long_duration_rays - 1, :] = field_dict['_FillValue']
            
            # Add field-specific metadata
            _add_field_metadata(field_dict, cfradial_name, ncvars)
            
            fields[cfradial_name] = field_dict
            
        else:
            print(f"Field {mmclx_name} not found in file")
    
    return fields

def _add_field_metadata(field_dict: Dict, cfradial_name: str, ncvars: Dict) -> None:
    """
    Add field-specific metadata and attributes.
    
    Args:
        field_dict: Field dictionary to modify
        cfradial_name: CF-Radial field name
        ncvars: NetCDF variables for extracting additional info
    """
    metadata_map = {
        "DBZ": {
            "long_name": "radar equivalent reflectivity factor",
            "standard_name": "equivalent_reflectivity_factor",
            "proposed_standard_name": "radar_equivalent_reflectivity_factor"
        },
        "VEL": {
            "long_name": "radial velocity of scatterers away from instrument", 
            "standard_name": "radial_velocity_of_scatterers_away_from_instrument",
            "field_folds": "true"
        },
        "WIDTH": {
            "long_name": "radar doppler spectrum width",
            "proposed_standard_name": "radar_doppler_spectrum_width"
        },
        "LDR": {
            "long_name": "radar linear depolarization ratio",
            "proposed_standard_name": "radar_linear_depolarization_ratio"  
        },
        "SNR": {
            "long_name": "radar signal to noise ratio",
            "proposed_standard_name": "radar_signal_to_noise_ratio"
        },
        "RHOHX": {
            "long_name": "co- to cross-polar correlation ratio for horizontal transmitted polarization",
            "proposed_standard_name": "co_to_cross_polar_correlation_ratio_h"
        },
        "PHIHX": {
            "long_name": "cross-polar differential phase",
            "proposed_standard_name": "cross_polar_differential_phase"
        }
    }
    
    if cfradial_name in metadata_map:
        field_dict.update(metadata_map[cfradial_name])
    
    # Add velocity folding information
    if cfradial_name == "VEL" and 'prf' in ncvars and 'lambda' in ncvars:
        vfold = ncvars['prf'][:] * ncvars['lambda'][:] / 4.0
        field_dict['field_limit_lower'] = -vfold
        field_dict['field_limit_upper'] = vfold

def _extract_instrument_parameters(ncobj: nc4.Dataset, nrays: int, filemetadata: FileMetadata) -> Dict:
    """
    Extract instrument parameters from netCDF object.
    
    Args:
        ncobj: Open netCDF4 Dataset
        nrays: Number of rays
        filemetadata: PyART FileMetadata object
        
    Returns:
        Dictionary of instrument parameters
    """
    ncvars = ncobj.variables
    instrument_parameters = {}
    
    # Pulse repetition time
    if "prf" in ncvars:
        prt_dict = filemetadata("prt")
        prt_value = float(1.0 / ncvars["prf"][:])  # Extract scalar, ensure float64 intermediate
        prt_dict["data"] = np.full((nrays,), prt_value, dtype=np.float32)
        instrument_parameters["prt"] = prt_dict
    
    # Parse hardware parameters from 'hrd' attribute
    hrd_variables = {}
    if 'hrd' in ncobj.ncattrs():
        hrd_attribute_value = ncobj.getncattr('hrd')
        hrd_variables = _parse_hrd_attribute(hrd_attribute_value)
    
    # Pulse width
    if "PULSE_WIDTH" in hrd_variables:
        pw_dict = filemetadata("pulse_width")
        pulse_width = float(hrd_variables["PULSE_WIDTH"])  # Ensure scalar
        pw_dict["data"] = np.full((nrays,), pulse_width, dtype=np.float32)
        instrument_parameters["pulse_width"] = pw_dict
    
    # Transmit power
    if "tpow" in ncvars:
        tp_dict = filemetadata("radar_measured_transmit_power_h")
        txpower = 10.0 * np.log10(ncvars["tpow"][:]) + 30  # Convert to dBm
        tp_dict["data"] = txpower.astype(np.float32)
        instrument_parameters["radar_measured_transmit_power_h"] = tp_dict
    
    return instrument_parameters

def _parse_hrd_attribute(hrd_attribute_value: str) -> Dict:
    """
    Parse the 'hrd' attribute string into a dictionary of hardware parameters.
    
    Args:
        hrd_attribute_value: String containing hardware parameters
        
    Returns:
        Dictionary of parsed hardware parameters
    """
    hrd_variables = {}
    
    if isinstance(hrd_attribute_value, str):
        for line in hrd_attribute_value.split('\n'):
            if not line.startswith("DESCR") and ':' in line:
                try:
                    key, value = line.split(':', 1)
                    hrd_variables[key.strip()] = convert_to_numerical(value.strip())
                except ValueError:
                    continue  # Skip malformed lines
    
    return hrd_variables

def convert_to_numerical(value: str) -> Union[float, int, str]:
    """
    Convert a string value to numerical type if possible.
    
    Args:
        value: String value to convert
        
    Returns:
        Converted value (float, int, or original string)
    """
    try:
        return float(value)
    except ValueError:
        try:
            return int(value)
        except ValueError:
            return value

# ===================
# CONVERSION ROUTINES  
# ===================

def convert_kepler_mmclx2l1(
    infile: str,
    outpath: str, 
    yaml_project_file: str,
    yaml_instrument_file: str,
    tracking_tag: str,
    data_version: str
) -> None:
    """
    Convert mmclx data from NCAS Mobile Ka-band Radar to Level 1 cfradial format.
    
    This function converts a single mmclx file to CF-Radial format compliant with
    the NCAS Radar Data Standard v1.0.0, adding appropriate metadata from YAML files.
    
    Args:
        infile: Full path to input mmclx file
        outpath: Directory where output CF-Radial file will be written  
        yaml_project_file: Path to YAML file with project metadata
        yaml_instrument_file: Path to YAML file with instrument metadata
        tracking_tag: AMOF tracking tag for the project
        data_version: Data version string (format: "n.m.p")
        
    Raises:
        FileNotFoundError: If input file or YAML files don't exist
        ValueError: If tracking tag not found in project file
        
    Example:
        >>> convert_kepler_mmclx2l1(
        ...     'data.mmclx', 
        ...     '/output/path',
        ...     'project.yml',
        ...     'instruments.yml', 
        ...     'AMOF_20220101000000',
        ...     '1.0.0'
        ... )
    """
    instrument_tagname = "ncas-mobile-ka-band-radar-1"
    
    # Load instrument metadata
    with open(yaml_instrument_file, "r") as stream:
        instruments = yaml.safe_load(stream)
    
    instrument = None
    for elem in instruments:
        if instrument_tagname in elem:
            instrument = elem[instrument_tagname]
            break
    
    if instrument is None:
        raise ValueError(f"Instrument {instrument_tagname} not found in {yaml_instrument_file}")
    
    # Load project metadata
    with open(yaml_project_file, "r") as stream:
        projects = yaml.safe_load(stream)
    
    project = None
    for p in projects:
        if tracking_tag in p:
            project = p[tracking_tag]
            break
    
    if project is None:
        raise ValueError(f"Tracking tag {tracking_tag} not found in {yaml_project_file}")
    
    # Read radar data
    radar_dataset = read_mira35_mmclx(infile)
    scan_name = radar_dataset.scan_name
    
    # Generate output filename
    radar_name = instrument["instrument_name"].lower()
    
    # Find location from project metadata
    _instr_key = 'chilbolton_instruments' if 'chilbolton_instruments' in project else 'ncas_instruments'
    for n in project[_instr_key]:
        if radar_name in n:
            project_instrument = n[radar_name]
            break
    
    location = project_instrument['platform']['location'].lower()
    
    file_timestamp = datetime.datetime.strptime(
        radar_dataset.metadata["time_coverage_start"], 
        '%Y-%m-%dT%H:%M:%SZ'
    )
    dtstr = file_timestamp.strftime('%Y%m%d-%H%M%S')
    
    outfile = os.path.join(
        outpath, 
        f'{radar_name}_{location}_{dtstr}_{scan_name.replace("_", "-", 1)}_l1_v{data_version}.nc'
    )
    
    # Write CF-Radial file
    pyart.io.write_cfradial(outfile, radar_dataset, format='NETCDF4', time_reference=True)
    
    # Add elevation and azimuth scan rates for MAN scans
    add_scan_rate_coordinates(outfile, radar_dataset)
    
    # Add NCAS metadata (revised_northangle=None will cause function to extract from YAML)
    cfradial_add_ncas_metadata(outfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version, revised_northangle=None)
    
    # Update history
    update_history_attribute(outfile, "Convert mmclx to CF-Radial with NCAS metadata")
    
    print(f"Successfully created {outfile}")

def update_history_attribute(filename: str, history_entry: str) -> None:
    """
    Update the history attribute of a NetCDF file.
    
    Args:
        filename: Path to NetCDF file
        history_entry: History entry to add
    """
    try:
        import netCDF4 as nc4
        import datetime
        
        with nc4.Dataset(filename, 'a') as ds:
            # Get current history if it exists
            current_history = ""
            if hasattr(ds, 'history'):
                current_history = ds.getncattr('history')
            
            # Create new history entry with timestamp and user info
            timestamp = datetime.datetime.utcnow().strftime('%a %b %d %H:%M:%S %Y')
            username = os.environ.get('USER', 'unknown')
            hostname = os.environ.get('HOSTNAME', 'unknown')
            
            # Replace degree symbols with "deg" in the history entry
            safe_history_entry = history_entry.replace('°', 'deg')
            
            new_entry = f"{timestamp} - user:{username} machine:{hostname} {safe_history_entry}"
            
            # Combine with existing history
            if current_history:
                updated_history = f"{new_entry}\n{current_history}"
            else:
                updated_history = new_entry
            
            # Set the updated history
            ds.setncattr('history', updated_history)
            
        print(f"Updated history attribute in {filename}")
        
    except Exception as e:
        print(f"Warning: Could not update history attribute: {e}")

def add_scan_rate_coordinates(filename: str, radar_obj: Radar) -> None:
    """
    Add elevation_scan_rate and azimuth_scan_rate coordinate variables to NetCDF file.
    
    This is needed for MAN (manual_rhi) scans where both elevation and azimuth
    change independently during aircraft tracking.
    
    Args:
        filename: Path to CF-Radial NetCDF file
        radar_obj: PyART Radar object containing the scan rate data
    """
    try:
        import netCDF4 as nc4
        
        # Check if radar object has the scan rate attributes
        has_elevation_rate = (hasattr(radar_obj, 'elevation_scan_rate') and 
                             radar_obj.elevation_scan_rate is not None and
                             'data' in radar_obj.elevation_scan_rate)
        has_azimuth_rate = (hasattr(radar_obj, 'azimuth_scan_rate') and 
                           radar_obj.azimuth_scan_rate is not None and
                           'data' in radar_obj.azimuth_scan_rate)
        
        if not (has_elevation_rate or has_azimuth_rate):
            return  # Nothing to add
        
        with nc4.Dataset(filename, 'a') as ds:
            # Add elevation_scan_rate if present
            if has_elevation_rate:
                if 'elevation_scan_rate' not in ds.variables:
                    var = ds.createVariable('elevation_scan_rate', 'f4', ('time',))
                else:
                    var = ds.variables['elevation_scan_rate']
                
                var[:] = radar_obj.elevation_scan_rate['data']
                var.units = radar_obj.elevation_scan_rate.get('units', 'degrees_per_second')
                var.long_name = radar_obj.elevation_scan_rate.get('long_name', 
                                                                   'antenna elevation angle scan rate')
                if 'standard_name' in radar_obj.elevation_scan_rate:
                    var.standard_name = radar_obj.elevation_scan_rate['standard_name']
                
                print(f"  Added elevation_scan_rate to {filename}")
            
            # Add azimuth_scan_rate if present
            if has_azimuth_rate:
                if 'azimuth_scan_rate' not in ds.variables:
                    var = ds.createVariable('azimuth_scan_rate', 'f4', ('time',))
                else:
                    var = ds.variables['azimuth_scan_rate']
                
                var[:] = radar_obj.azimuth_scan_rate['data']
                var.units = radar_obj.azimuth_scan_rate.get('units', 'degrees_per_second')
                var.long_name = radar_obj.azimuth_scan_rate.get('long_name', 
                                                                 'antenna azimuth angle scan rate')
                if 'standard_name' in radar_obj.azimuth_scan_rate:
                    var.standard_name = radar_obj.azimuth_scan_rate['standard_name']
                
                print(f"  Added azimuth_scan_rate to {filename}")
                
    except Exception as e:
        print(f"Warning: Could not add scan rate coordinates: {e}")

def cfradial_get_bbox(cfradfile: str) -> str:
    """
    Calculate bounding box from CF-Radial file.
    
    Args:
        cfradfile: Path to CF-Radial file
        
    Returns:
        String describing the bounding box
    """
    print(cfradfile)
    radar = pyart.io.read_cfradial(cfradfile)
    latmin = np.min(radar.gate_latitude['data'])
    lonmin = np.min(radar.gate_longitude['data'])
    latmax = np.max(radar.gate_latitude['data'])
    lonmax = np.max(radar.gate_longitude['data'])
    print(latmin, latmax, lonmin, lonmax)
    boundingbox = f"Bounding box: {latmin:.2f}N {lonmin:.2f}E, {latmax:.2f}N {lonmax:.2f}E"
    return boundingbox

def cfradial_add_instrument_parameters(
    mmclxfile: str, 
    cfradfile: str, 
    yaml_project_file: str, 
    yaml_instrument_file: str, 
    tracking_tag: str, 
    data_version: str
) -> None:
    """
    Add instrument parameters to CF-Radial file from original mmclx file.
    
    Args:
        mmclxfile: Path to original mmclx file
        cfradfile: Path to CF-Radial file to modify
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        tracking_tag: AMOF tracking tag
        data_version: Data version string
    """
    with nc4.Dataset(cfradfile, 'r+') as ds:
        print(cfradfile)
        
        if mmclxfile.endswith('.mmclx.gz'):
            with gzip.open(mmclxfile) as gz:
                dsin = nc4.Dataset('dummy', mode='r', memory=gz.read())
        else:
            dsin = nc4.Dataset(mmclxfile, 'r')
        
        try:
            print('Creating frequency dimension')
            if 'frequency' not in ds.dimensions:
                frequency = ds.createDimension("frequency", 1)
            
            if 'lambda' in dsin.variables:
                varin = dsin['lambda']
                lightspeed = 299792458
                tx_freq = lightspeed / varin[:]
                print(tx_freq)
                
                if 'frequency' not in ds.variables:
                    varout = ds.createVariable('frequency', varin.datatype, ("frequency",))
                    varout.standard_name = 'radiation_frequency'
                    varout.long_name = 'frequency of transmitted radiation'
                    varout.units = 's-1'
                    varout[:] = tx_freq
                    print('Creating meta_group attribute')
                    varout.meta_group = "instrument_parameters"
            
            if 'radar_measured_transmit_power_h' in ds.variables:
                varout = ds['radar_measured_transmit_power_h']
                varout.long_name = "radar_measured_transmit_power_h"
                varout.units = 'dBm'
                varout.meta_group = "radar_parameters"
                
        finally:
            dsin.close()

def cfradial_add_geometry_correction(cfradfile: str, revised_northangle: float) -> None:
    """
    Add geometry correction information to CF-Radial file.
    
    Args:
        cfradfile: Path to CF-Radial file to modify
        revised_northangle: North angle correction applied
    """
    with nc4.Dataset(cfradfile, 'r+') as ds:
        print(cfradfile)
        print(f'Revised northangle = {revised_northangle}')
        
        print('creating azimuth_correction')
        if 'azimuth_correction' not in ds.variables:
            varout = ds.createVariable('azimuth_correction', 'float32')
            print('created')
            varout.long_name = "azimuth correction applied"
            varout.units = 'degrees'
            varout.meta_group = "geometry_correction"
            varout.comment = "Azimuth correction applied. North angle relative to instrument home azimuth."
            varout[:] = revised_northangle

def amend_unitless(cfradfile: str) -> None:
    """
    Replace 'unitless' and 'count' units with empty strings.
    
    Args:
        cfradfile: Path to CF-Radial file to modify
    """
    with nc4.Dataset(cfradfile, 'r+') as ds:
        # Loop through each variable in the NetCDF file
        for var_name in ds.variables:
            var = ds.variables[var_name]
            if hasattr(var, 'units') and var.units == 'unitless':
                var.units = ""
            if hasattr(var, 'units') and var.units == 'count':
                var.units = ""

def lowercase_long_names(cfradfile: str) -> None:
    """
    Convert long_name attributes to lowercase (except UTC).
    
    Args:
        cfradfile: Path to CF-Radial file to modify
    """
    with nc4.Dataset(cfradfile, 'r+') as ds:
        # Loop through each variable in the NetCDF file
        for var_name in ds.variables:
            var = ds.variables[var_name]
            if hasattr(var, 'long_name'):
                var.long_name = var.long_name.lower()
                if "utc" in var.long_name:
                    var.long_name = var.long_name.replace("utc", "UTC")

def time_long_name(cfradfile: str) -> None:
    """
    Update time variable long_name if time_reference exists.
    
    Args:
        cfradfile: Path to CF-Radial file to modify
    """
    with nc4.Dataset(cfradfile, 'r+') as ds:
        time_var = ds.variables['time']
        if 'time_reference' in ds.variables:
            time_var.long_name = "time_since_time_reference"

def find_mmclxfiles(
    start_time: str, 
    end_time: str, 
    sweep_type: str, 
    inpath: str, 
    gzip_flag: bool = False
) -> List[str]:
    """
    Find mmclx files within a time range and sweep type.
    
    Args:
        start_time: Start time as 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'
        end_time: End time as 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'
        sweep_type: Type of sweep ('rhi', 'ppi', 'vert', 'vad')
        inpath: Directory path to search
        gzip_flag: Whether files are gzip compressed
        
    Returns:
        List of matching file paths
    """
    # Handle both time formats
    def parse_time_string(time_str):
        """Parse time string in either format"""
        try:
            # Try ISO format first
            return datetime.datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ')
        except ValueError:
            try:
                # Try space-separated format
                return datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
            except ValueError:
                raise ValueError(f"Time string '{time_str}' does not match expected formats: "
                               "'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'")
    
    # Convert the input times to datetime objects
    start_datetime = parse_time_string(start_time)
    end_datetime = parse_time_string(end_time)
    
    print(f"Searching for {sweep_type.upper()} files from {start_time} to {end_time}")
    
    # Get the date string from start time
    date_str = start_datetime.strftime('%Y%m%d')
    
    # Check if we need to add date subdirectory or use path directly
    # Try date subdirectory first (standard case for most campaigns)
    # If that doesn't exist, use inpath directly (WOEST case with scan-type subdirs)
    inpath_obj = Path(inpath)
    date_subdir = inpath_obj / date_str
    
    if date_subdir.exists():
        # Standard case: files organized in YYYYMMDD subdirectories
        search_path = date_subdir
        print(f"Looking in date directory: {search_path}")
    elif inpath_obj.exists():
        # WOEST case: path already points to scan-type directory
        search_path = inpath_obj
        print(f"Using path directly (no date subdir): {search_path}")
    else:
        # Neither exists
        print(f"Path does not exist: {inpath_obj} (with or without date {date_str})")
        return []
    
    # Define search patterns based on sweep type
    if gzip_flag:
        patterns = {
            'rhi': "*rhi*.mmclx.gz",
            'ppi': "*ppi*.mmclx.gz", 
            'vert': "*vert*.mmclx.gz",
            'vad': "*ppi*.mmclx.gz",  # VAD files are typically PPI scans
            'man': "*man*.mmclx.gz"   # MAN (manual tracking) scans
        }
    else:
        patterns = {
            'rhi': "*rhi*.mmclx",
            'ppi': "*ppi*.mmclx",
            'vert': "*vert*.mmclx", 
            'vad': "*ppi*.mmclx",  # VAD files are typically PPI scans
            'man': "*man*.mmclx"   # MAN (manual tracking) scans
        }
    
    pattern = patterns.get(sweep_type.lower(), f"*{sweep_type}*.mmclx{'.gz' if gzip_flag else ''}")
    
    # Find candidate files
    candidate_files = list(search_path.glob(pattern))
    print(f"Found {len(candidate_files)} candidate {sweep_type.upper()} files")
    
    # Filter files by time range
    matching_files = []
    
    for file_path in candidate_files:
        try:
            if gzip_flag:
                with gzip.open(file_path, 'rb') as gz:
                    with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                        if 'time' in nc.dimensions and len(nc.dimensions['time']) > 0:
                            file_time = _safe_parse_time(nc, 0)
                            if file_time is None:
                                print(f"Could not parse time from {file_path.name}, skipping")
                                continue
                                
                            if start_datetime <= file_time <= end_datetime:
                                matching_files.append(str(file_path))
                                print(f"Added {sweep_type.upper()} file: {file_path.name} (time: {file_time})")
            else:
                with nc4.Dataset(file_path, 'r') as nc:
                    if 'time' in nc.dimensions and len(nc.dimensions['time']) > 0:
                        file_time = _safe_parse_time(nc, 0)
                        if file_time is None:
                            print(f"Could not parse time from {file_path.name}, skipping")
                            continue
                            
                        if start_datetime <= file_time <= end_datetime:
                            matching_files.append(str(file_path))
                            print(f"Added {sweep_type.upper()} file: {file_path.name} (time: {file_time})")
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            continue
    
    print(f"Found {len(matching_files)} matching {sweep_type.upper()} files")
    return sorted(matching_files)

def find_mmclx_rhi_files(
    start_time: str,
    end_time: str,
    azim_min: float,
    azim_max: float,
    inpath: str,
    gzip_flag: bool = False,
    azimuth_offset: float = -6.85,
    revised_northangle: float = 302.15
) -> List[str]:
    """
    Find RHI mmclx files within time and azimuth ranges.
    
    Args:
        start_time: Start time as 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'
        end_time: End time as 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'
        azim_min: Minimum azimuth angle
        azim_max: Maximum azimuth angle
        inpath: Directory path to search
        gzip_flag: Whether files are gzip compressed
        azimuth_offset: Azimuth offset for searching
        revised_northangle: North angle correction
        
    Returns:
        List of matching RHI file paths
    """
    # Handle both time formats
    def parse_time_string(time_str):
        """Parse time string in either format"""
        try:
            # Try ISO format first
            return datetime.datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ')
        except ValueError:
            try:
                # Try space-separated format
                return datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
            except ValueError:
                raise ValueError(f"Time string '{time_str}' does not match expected formats: "
                               "'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'")
    
    # Convert the input times to datetime objects
    start_datetime = parse_time_string(start_time)
    end_datetime = parse_time_string(end_time)
    
    print(f"Searching for RHI files from {start_time} to {end_time}")
    print(f"Azimuth range: {azim_min}° to {azim_max}°")
    
    # Get the date string from start time to look in the right subdirectory
    date_str = start_datetime.strftime('%Y%m%d')
    
    # Try date subdirectory first (standard campaigns); fall back to inpath directly (WOEST)
    inpath_obj = Path(inpath)
    date_subdir = inpath_obj / date_str
    if date_subdir.exists():
        search_path = date_subdir
        print(f"Looking in date directory: {search_path}")
    elif inpath_obj.exists():
        search_path = inpath_obj
        print(f"Using existing path directly (no date subdir): {search_path}")
    else:
        print(f"Path does not exist: {inpath_obj} (with or without date {date_str})")
        return []
    
    if not search_path.exists():
        print(f"Search path does not exist: {search_path}")
        return []
    
    # Define a list to store the found files
    matching_files = []
    
    az_search_offset = -8.0  # Before July
    
    # Look for RHI files (use rglob for recursive search in subdirectories like WOEST's hsrhi/)
    if gzip_flag:
        rhi_pattern = "*rhi*.mmclx.gz"
    else:
        rhi_pattern = "*rhi*.mmclx"
    
    candidate_files = list(search_path.rglob(rhi_pattern))
    print(f"Found {len(candidate_files)} candidate RHI files")
    
    # Check each file's timestamp and azimuth
    for file_path in candidate_files:
        try:
            if gzip_flag:
                with gzip.open(file_path, 'rb') as gz:
                    with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                        nrays = len(nc.dimensions.get('time', []))
                        if nrays > 0:
                            # Use a safe index - don't assume index 3 exists
                            azim_index = min(3, nrays - 1)  # Use index 3 or last available
                            
                            # Use safe time parsing with microseconds
                            file_time = _safe_parse_time(nc, 0)
                            if file_time is None:
                                print(f"Could not parse time from {file_path.name}, skipping")
                                continue
                                
                            azim = (nc['azi'][azim_index] + revised_northangle) % 360
                            
                            if start_datetime <= file_time <= end_datetime:
                                if azim_min <= convert_angle(azim, az_search_offset) < azim_max:
                                    print(f'{file_path.name}: {file_time} {azim_min} {convert_angle(azim, az_search_offset)} {azim_max}')
                                    matching_files.append(str(file_path))
            else:
                with nc4.Dataset(file_path, 'r') as nc:
                    nrays = len(nc.dimensions.get('time', []))
                    if nrays > 0:
                        # Use a safe index - don't assume index 3 exists
                        azim_index = min(3, nrays - 1)  # Use index 3 or last available
                        
                        # Use safe time parsing with microseconds
                        file_time = _safe_parse_time(nc, 0)
                        if file_time is None:
                            print(f"Could not parse time from {file_path.name}, skipping")
                            continue
                            
                        azim = (nc['azi'][azim_index] + revised_northangle) % 360
                        
                        if start_datetime <= file_time <= end_datetime:
                            if azim_min <= convert_angle(azim, az_search_offset) < azim_max:
                                print(f'{file_path.name}: {file_time} {azim_min} {convert_angle(azim, az_search_offset)} {azim_max}')
                                matching_files.append(str(file_path))
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            continue
                    
    print(f"Found {len(matching_files)} matching RHI files")
    return sorted(matching_files)

def convert_angle(angle: float, offset: float) -> float:
    """
    Convert angle with offset, handling 360-degree wraparound.
    
    Args:
        angle: Input angle in degrees
        offset: Offset to apply
        
    Returns:
        Converted angle
    """
    print(angle, offset)
    if angle >= 360 + round(offset):
        angle -= 360
    return angle

def find_mmclx_ppi_files(
    start_time: str,
    end_time: str,
    elev_min: float,
    elev_max: float,
    inpath: str,
    gzip_flag: bool = False,
    azimuth_offset: float = 0.0,
    revised_northangle: float = 55.7
) -> List[str]:
    """
    Find PPI mmclx files within time and elevation ranges.
    
    Args:
        start_time: Start time as 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'
        end_time: End time as 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'
        elev_min: Minimum elevation angle
        elev_max: Maximum elevation angle
        inpath: Directory path to search
        gzip_flag: Whether files are gzip compressed
        azimuth_offset: Azimuth offset for searching
        revised_northangle: North angle correction
        
    Returns:
        List of matching PPI file paths
    """
    # Handle both time formats
    def parse_time_string(time_str):
        """Parse time string in either format"""
        try:
            # Try ISO format first
            return datetime.datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ')
        except ValueError:
            try:
                # Try space-separated format
                return datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
            except ValueError:
                raise ValueError(f"Time string '{time_str}' does not match expected formats")
    
    # Convert the input times to datetime objects
    start_datetime = parse_time_string(start_time)
    end_datetime = parse_time_string(end_time)
    
    print(f"Searching for PPI files from {start_time} to {end_time}")
    print(f"Elevation range: {elev_min}° to {elev_max}°")
    
    # Get the date string from start time
    date_str = start_datetime.strftime('%Y%m%d')
    
    # Try date subdirectory first (standard campaigns); fall back to inpath directly (WOEST)
    inpath_obj = Path(inpath)
    date_subdir = inpath_obj / date_str
    if date_subdir.exists():
        search_path = date_subdir
        print(f"Looking in date directory: {search_path}")
    elif inpath_obj.exists():
        search_path = inpath_obj
        print(f"Using existing path directly (no date subdir): {search_path}")
    else:
        print(f"Path does not exist: {inpath_obj} (with or without date {date_str})")
        return []
    
    if not search_path.exists():
        print(f"Date directory does not exist: {search_path}")
        return []
    
    # Define a list to store the found files
    matching_files = []
    
    # Look for PPI files
    if gzip_flag:
        ppi_pattern = "*ppi*.mmclx.gz"
    else:
        ppi_pattern = "*ppi*.mmclx"
    
    candidate_files = list(search_path.glob(ppi_pattern))
    print(f"Found {len(candidate_files)} candidate PPI files")
    
    # Check each file's timestamp and elevation
    for file_path in candidate_files:
        try:
            if gzip_flag:
                with gzip.open(file_path, 'rb') as gz:
                    with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                        nrays = len(nc.dimensions.get('time', []))
                        if nrays > 0:
                            # Use safe time parsing
                            file_time = _safe_parse_time(nc, 0)
                            if file_time is None:
                                print(f"Could not parse time from {file_path.name}, skipping")
                                continue
                                
                            # Check elevation angle
                            if 'elv' in nc.variables:
                                elev = nc['elv'][0]  # First elevation angle
                                
                                if (start_datetime <= file_time <= end_datetime and 
                                    elev_min <= elev <= elev_max):
                                    matching_files.append(str(file_path))
                                    print(f"Added PPI file: {file_path.name} (elev: {elev:.1f}°)")
            else:
                with nc4.Dataset(file_path, 'r') as nc:
                    nrays = len(nc.dimensions.get('time', []))
                    if nrays > 0:
                        # Use safe time parsing
                        file_time = _safe_parse_time(nc, 0)
                        if file_time is None:
                            print(f"Could not parse time from {file_path.name}, skipping")
                            continue
                            
                        # Check elevation angle
                        if 'elv' in nc.variables:
                            elev = nc['elv'][0]  # First elevation angle
                            
                            if (start_datetime <= file_time <= end_datetime and 
                                elev_min <= elev <= elev_max):
                                matching_files.append(str(file_path))
                                print(f"Added PPI file: {file_path.name} (elev: {elev:.1f}°)")
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            continue
    
    print(f"Found {len(matching_files)} matching PPI files")
    return sorted(matching_files)

def cfradial_add_ncas_metadata(
    cfradial_file: str,
    yaml_project_file: str, 
    yaml_instrument_file: str,
    tracking_tag: str,
    data_version: str,
    revised_northangle: float = None
) -> None:
    """
    Add NCAS metadata to a CF-Radial file.
    
    Args:
        cfradial_file: Path to CF-Radial file to modify
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file  
        tracking_tag: AMOF tracking tag
        data_version: Data version string
        revised_northangle: North angle correction (if None, will try to extract from YAML)
        
    Raises:
        FileNotFoundError: If YAML files don't exist
        ValueError: If YAML files are invalid or missing required content
    """
    print(f"Adding NCAS metadata to {cfradial_file}")
    print(f"Using project file: {yaml_project_file}")
    print(f"Using instrument file: {yaml_instrument_file}")
    print(f"Tracking tag: {tracking_tag}")
    print(f"Data version: {data_version}")
    
    # Validate that YAML files exist BEFORE attempting any processing
    if not os.path.exists(yaml_project_file):
        raise FileNotFoundError(f"Project YAML file not found: {yaml_project_file}")
    if not os.path.exists(yaml_instrument_file):
        raise FileNotFoundError(f"Instrument YAML file not found: {yaml_instrument_file}")
    
    # Add NCAS metadata manually by reading YAML project and instrument files
    # The ncas_amof_netcdf_template library expects a different workflow (single metadata CSV/YAML),
    # so we use our custom implementation that properly handles project + instrument YAML files
    try:
        _add_ncas_metadata_manually(cfradial_file, yaml_project_file, yaml_instrument_file, tracking_tag, data_version, revised_northangle)
        print(f"Successfully added NCAS metadata")
    except Exception as e:
        raise RuntimeError(f"Error adding NCAS metadata to {cfradial_file}: {e}")


def _add_ncas_metadata_manually(
    cfradial_file: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    tracking_tag: str,
    data_version: str,
    revised_northangle: float = None
) -> None:
    """
    Manually add NCAS metadata to a CF-Radial file by reading YAML files and updating attributes.
    
    Raises:
        FileNotFoundError: If YAML files don't exist
        ValueError: If YAML files are invalid or missing required content
        RuntimeError: If NetCDF file cannot be modified
    """
    import netCDF4 as nc4
    import yaml
    from datetime import datetime
    
    # Load YAML metadata files - FAIL if they don't exist or are invalid
    try:
        with open(yaml_project_file, 'r') as f:
            project_metadata = yaml.safe_load(f)
        print(f"Loaded project metadata from {yaml_project_file}")
        
        if project_metadata is None:
            raise ValueError(f"Project YAML file is empty or invalid: {yaml_project_file}")
            
    except FileNotFoundError:
        raise FileNotFoundError(f"Project YAML file not found: {yaml_project_file}")
    except yaml.YAMLError as e:
        raise ValueError(f"Invalid YAML syntax in project file {yaml_project_file}: {e}")
    except Exception as e:
        raise RuntimeError(f"Error loading project YAML file {yaml_project_file}: {e}")
    
    try:
        with open(yaml_instrument_file, 'r') as f:
            instrument_metadata = yaml.safe_load(f)
        print(f"Loaded instrument metadata from {yaml_instrument_file}")
        
        if instrument_metadata is None:
            raise ValueError(f"Instrument YAML file is empty or invalid: {yaml_instrument_file}")
            
    except FileNotFoundError:
        raise FileNotFoundError(f"Instrument YAML file not found: {yaml_instrument_file}")
    except yaml.YAMLError as e:
        raise ValueError(f"Invalid YAML syntax in instrument file {yaml_instrument_file}: {e}")
    except Exception as e:
        raise RuntimeError(f"Error loading instrument YAML file {yaml_instrument_file}: {e}")
    
    # Debug: Print the structure of loaded YAML files
    print(f"Project metadata type: {type(project_metadata)}")
    print(f"Instrument metadata type: {type(instrument_metadata)}")
    
    if isinstance(project_metadata, list):
        print(f"Project metadata is a list with {len(project_metadata)} items")
        if project_metadata:
            print(f"First item keys: {list(project_metadata[0].keys()) if isinstance(project_metadata[0], dict) else 'Not a dict'}")
    elif isinstance(project_metadata, dict):
        print(f"Project metadata keys: {list(project_metadata.keys())}")
    
    if isinstance(instrument_metadata, list):
        print(f"Instrument metadata is a list with {len(instrument_metadata)} items")
        if instrument_metadata:
            print(f"First item keys: {list(instrument_metadata[0].keys()) if isinstance(instrument_metadata[0], dict) else 'Not a dict'}")
    elif isinstance(instrument_metadata, dict):
        print(f"Instrument metadata keys: {list(instrument_metadata.keys())}")
    
    # Handle different YAML file structures flexibly
    project_info = None
    
    # Handle project metadata (could be list or dict)
    if isinstance(project_metadata, list):
        # Search through list for project info
        for item in project_metadata:
            if isinstance(item, dict):
                # Look for tracking tag
                if tracking_tag in item:
                    tag_data = item[tracking_tag]
                    if isinstance(tag_data, dict):
                        if 'project' in tag_data:
                            project_info = tag_data['project']
                        else:
                            project_info = tag_data
                    break
                # Look for direct project info
                elif 'project' in item:
                    project_info = item['project']
                    break
                # Look for project-like content
                elif 'title' in item or 'principal_investigator' in item:
                    project_info = item
                    break
    elif isinstance(project_metadata, dict):
        # Original dict handling logic
        if 'project' in project_metadata:
            project_info = project_metadata['project']
        elif tracking_tag in project_metadata:
            tag_data = project_metadata[tracking_tag]
            if isinstance(tag_data, dict) and 'project' in tag_data:
                project_info = tag_data['project']
            else:
                project_info = tag_data
        else:
            for key, value in project_metadata.items():
                if isinstance(value, dict) and ('title' in value or 'principal_investigator' in value):
                    project_info = value
                    print(f"Found project info under key: {key}")
                    break
    
    # Handle instrument metadata (could be list or dict)
    instrument_info = None
    
    if isinstance(instrument_metadata, list):
        # Search through list for instrument info
        for item in instrument_metadata:
            if isinstance(item, dict):
                # Look for direct instrument match
                if item.get('name') == 'ncas-mobile-ka-band-radar-1':
                    instrument_info = item
                    break
                # Look for nested instrument structure
                elif 'ncas-mobile-ka-band-radar-1' in item:
                    instrument_info = item['ncas-mobile-ka-band-radar-1']
                    break
                # Look for instruments list
                elif 'instruments' in item:
                    instruments_list = item['instruments']
                    if isinstance(instruments_list, list):
                        for instrument in instruments_list:
                            if isinstance(instrument, dict) and instrument.get('name') == 'ncas-mobile-ka-band-radar-1':
                                instrument_info = instrument
                                break
                    if instrument_info:
                        break
    elif isinstance(instrument_metadata, dict):
        # Original dict handling logic
        if 'instruments' in instrument_metadata:
            instruments_list = instrument_metadata['instruments']
            if isinstance(instruments_list, list):
                for instrument in instruments_list:
                    if isinstance(instrument, dict) and instrument.get('name') == 'ncas-mobile-ka-band-radar-1':
                        instrument_info = instrument
                        break
        else:
            for key, value in instrument_metadata.items():
                if isinstance(value, dict):
                    if value.get('name') == 'ncas-mobile-ka-band-radar-1':
                        instrument_info = value
                        break
                    elif 'ncas-mobile-ka-band-radar-1' in str(key).lower():
                        instrument_info = value
                        break
    
    # Also look for instrument info in the project YAML file under ncas_instruments/chilbolton_instruments
    project_instrument_info = None
    _KEPLER_RADAR_NAMES = ('ncas-mobile-ka-band-radar-1', 'reading-mobile-ka-band-radar-1')

    # Check if project_info has ncas_instruments or chilbolton_instruments
    _instr_key = next((k for k in ('chilbolton_instruments', 'ncas_instruments') if isinstance(project_info, dict) and k in project_info), None)
    if project_info and isinstance(project_info, dict) and _instr_key:
        instruments_list = project_info[_instr_key]
        if isinstance(instruments_list, list):
            for instrument in instruments_list:
                if isinstance(instrument, dict):
                    # Handle dict with instrument name as key
                    for inst_name, inst_data in instrument.items():
                        # Check for exact match or partial match (case insensitive)
                        if inst_name in _KEPLER_RADAR_NAMES or 'ka-band-radar-1' in inst_name.lower():
                            project_instrument_info = inst_data
                            print(f"Found instrument info under key: {inst_name}")
                            break
                if project_instrument_info:
                    break
    
    print(f"Found project info: {project_info is not None}")
    print(f"Found instrument info: {instrument_info is not None}")
    print(f"Found project instrument info: {project_instrument_info is not None}")
    
    if project_info:
        print(f"Project info keys: {list(project_info.keys())}")
    if instrument_info:
        print(f"Instrument info keys: {list(instrument_info.keys())}")
    if project_instrument_info:
        print(f"Project instrument info keys: {list(project_instrument_info.keys())}")
    
    # Extract north_angle from project metadata for azimuth_correction
    # Use the parameter if provided, otherwise try to extract from YAML
    if revised_northangle is None:
        revised_northangle = 55.9  # Default value
        if project_instrument_info and 'north_angle' in project_instrument_info:
            revised_northangle = float(project_instrument_info['north_angle'])
            print(f"Using north_angle from project YAML: {revised_northangle}°")
        else:
            print(f"Using default north_angle: {revised_northangle}°")
    else:
        print(f"Using north_angle from parameter: {revised_northangle}°")
    
    # Open NetCDF file and update attributes
    try:
        with nc4.Dataset(cfradial_file, 'a') as ds:
            # Remove PyART's default "version" attribute
            if hasattr(ds, 'version'):
                ds.delncattr('version')
                print("Removed PyART 'version' attribute")
            
            # Add time coverage attributes FIRST
            if 'time' in ds.variables:
                time_var = ds.variables['time']
                if hasattr(time_var, 'units') and len(time_var) > 0:
                    try:
                        start_time = cftime.num2pydate(time_var[0], time_var.units)
                        end_time = cftime.num2pydate(time_var[-1], time_var.units)
                        
                        ds.setncattr('time_coverage_start', start_time.strftime('%Y-%m-%dT%H:%M:%SZ'))
                        ds.setncattr('time_coverage_end', end_time.strftime('%Y-%m-%dT%H:%M:%SZ'))
                        
                        print(f"Added time coverage: {start_time.strftime('%Y-%m-%dT%H:%M:%SZ')} to {end_time.strftime('%Y-%m-%dT%H:%M:%SZ')}")
                        
                    except Exception as e:
                        print(f"Warning: Could not add time coverage: {e}")
            
            # Update Conventions
            ds.setncattr('Conventions', 'NCAS-Radar-1.0 CfRadial-1.4 instrument_parameters radar_parameters geometry_correction')

            # Set CF featureType for vertical-pointing files
            _is_vpt = False
            if 'sweep_mode' in ds.variables:
                _sm = ds.variables['sweep_mode'][:]
                _sm_str = b''.join(_sm.flatten()).decode('ascii', errors='ignore').strip('\x00 ')
                if 'vertical_pointing' in _sm_str:
                    _is_vpt = True
            if _is_vpt:
                ds.setncattr('featureType', 'timeSeriesProfile')
                print("Set featureType = 'timeSeriesProfile' (VPT file)")

            # Update core attributes - set defaults first
            ds.setncattr('title', 'Moment data from NCAS Mobile Ka-band Radar (Kepler)')
            ds.setncattr('institution', 'National Centre for Atmospheric Science (NCAS) and Science and Technology Facilities Council (STFC) as part of UK Research and Innovation (UKRI)')
            ds.setncattr('source', 'NCAS Ka-Band Mobile Cloud Radar unit 1')
            # Save any comment already written by the reader (e.g. logp reconstruction note)
            # so we can append it after the YAML project comment is applied.
            _reader_comment = ''
            try:
                _existing = ds.getncattr('comment')
                if _existing and _existing.strip():
                    _reader_comment = _existing.strip()
            except AttributeError:
                pass
            ds.setncattr('references', 'www.metek.de')
            
            # Override source, title, references from YAML if available
            # Check project_instrument_info first (from project YAML)
            if project_instrument_info:
                if 'source' in project_instrument_info:
                    ds.setncattr('source', project_instrument_info['source'])
                    print(f"Set source from project YAML: {project_instrument_info['source']}")
                if 'title' in project_instrument_info:
                    ds.setncattr('title', project_instrument_info['title'])
                    print(f"Set title from project YAML: {project_instrument_info['title']}")
                if 'comment' in project_instrument_info:
                    _yaml_comment = project_instrument_info['comment'].rstrip()
                    # Extract any extra content (e.g. logp reconstruction note) from the
                    # existing file comment that goes beyond the YAML comment.
                    # This is idempotent: if cfradial_add_ncas_metadata has already been
                    # called once (e.g. inside multi_mmclx2cfrad), _reader_comment will
                    # start with the YAML comment followed by the extra note; strip the
                    # YAML prefix to avoid duplicating it.
                    _extra = ''
                    if _reader_comment:
                        _stripped_reader = _reader_comment.strip()
                        _stripped_yaml = _yaml_comment.strip()
                        if _stripped_reader.startswith(_stripped_yaml):
                            _extra = _stripped_reader[len(_stripped_yaml):].strip()
                        elif _stripped_reader != _stripped_yaml:
                            _extra = _stripped_reader
                    if _extra:
                        ds.setncattr('comment', (_yaml_comment + '\n' + _extra).encode('utf-8'))
                    else:
                        ds.setncattr('comment', _yaml_comment.encode('utf-8'))
                    print(f"Set comment from project YAML")
                else:
                    # project_instrument_info present but no 'comment' key — preserve
                    # any reader comment (e.g. logp note) or fall back to a blank.
                    _fb = _reader_comment if _reader_comment else ' '
                    ds.setncattr('comment', _fb.encode('utf-8'))
            else:
                # No project_instrument_info — preserve any reader comment or fall back.
                _fb = _reader_comment if _reader_comment else ' '
                ds.setncattr('comment', _fb.encode('utf-8'))
            
            # Check instrument_info (from instrument YAML) as fallback
            if instrument_info:
                if 'source' in instrument_info and not (project_instrument_info and 'source' in project_instrument_info):
                    ds.setncattr('source', instrument_info['source'])
                    print(f"Set source from instrument YAML: {instrument_info['source']}")
                if 'references' in instrument_info:
                    ds.setncattr('references', instrument_info['references'])
                    print(f"Set references from instrument YAML: {instrument_info['references']}")
            ds.setncattr('product_version', f'v{data_version}')
            ds.setncattr('processing_level', '1')
            ds.setncattr('licence', 'This dataset is released for use under the Open Government Licence, OGL-UK-3.0 (see https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).')
            
            # Add creator information
            ds.setncattr('creator_name', 'Chris Walden')
            ds.setncattr('creator_email', 'chris.walden@ncas.ac.uk')
            ds.setncattr('creator_url', 'https://orcid.org/0000-0002-5718-466X')
            
            # Add instrument metadata - always set defaults, then override from YAML if available
            ds.setncattr('instrument_manufacturer', 'Meteorologische Messtechnik (Metek) GmbH')
            ds.setncattr('instrument_model', 'MIRA-35S')
            ds.setncattr('instrument_serial_number', 'BX3C')
            ds.setncattr('instrument_pid', 'https://hdl.handle.net/21.12132/3.b7a8298b7a54405f')
            
            if instrument_info:
                print("Overriding instrument metadata from YAML")
                # Check for both prefixed (instrument_*) and non-prefixed key names
                if 'instrument_manufacturer' in instrument_info:
                    ds.setncattr('instrument_manufacturer', instrument_info['instrument_manufacturer'])
                elif 'manufacturer' in instrument_info:
                    ds.setncattr('instrument_manufacturer', instrument_info['manufacturer'])
                
                if 'instrument_model' in instrument_info:
                    ds.setncattr('instrument_model', instrument_info['instrument_model'])
                elif 'model' in instrument_info:
                    ds.setncattr('instrument_model', instrument_info['model'])
                
                if 'instrument_serial_number' in instrument_info:
                    ds.setncattr('instrument_serial_number', instrument_info['instrument_serial_number'])
                elif 'serial_number' in instrument_info:
                    ds.setncattr('instrument_serial_number', instrument_info['serial_number'])
                
                if 'instrument_pid' in instrument_info:
                    ds.setncattr('instrument_pid', instrument_info['instrument_pid'])
                elif 'pid' in instrument_info:
                    ds.setncattr('instrument_pid', instrument_info['pid'])
            else:
                print("Using default instrument metadata")
            
            # Add project metadata - always set default, then override from YAML if available
            ds.setncattr('project', 'Contrail Observations and Lifecycle Tracking (COBALT)')
            
            if project_info:
                print("Overriding project metadata from YAML")
                # Check for project_name first (standard key), then title (fallback)
                if 'project_name' in project_info:
                    ds.setncattr('project', project_info['project_name'])
                elif 'title' in project_info:
                    ds.setncattr('project', project_info['title'])
                
                if 'principal_investigator' in project_info:
                    pi = project_info['principal_investigator']
                    if isinstance(pi, dict):
                        if 'name' in pi:
                            ds.setncattr('project_principal_investigator', pi['name'])
                        if 'email' in pi:
                            ds.setncattr('project_principal_investigator_email', pi['email'])
                        if 'pid' in pi:
                            ds.setncattr('project_principal_investigator_url', pi['pid'])
                
                # Try to get acknowledgement from project_instrument_info first (from project YAML)
                acknowledgement_set = False
                if project_instrument_info:
                    if 'acknowledgement' in project_instrument_info:
                        ds.setncattr('acknowledgement', project_instrument_info['acknowledgement'])
                        acknowledgement_set = True
                        print("Set acknowledgement from project YAML instrument info")
                    elif 'acknowledgment' in project_instrument_info:  # Alternative spelling
                        ds.setncattr('acknowledgement', project_instrument_info['acknowledgment'])
                        acknowledgement_set = True
                        print("Set acknowledgement from project YAML instrument info (alternative spelling)")
                
                # Fall back to project_info if not found in instrument info
                if not acknowledgement_set:
                    if 'acknowledgement' in project_info:
                        ds.setncattr('acknowledgement', project_info['acknowledgement'])
                        print("Set acknowledgement from project info")
                    elif 'acknowledgment' in project_info:  # Alternative spelling
                        ds.setncattr('acknowledgement', project_info['acknowledgment'])
                        print("Set acknowledgement from project info (alternative spelling)")
            else:
                print("Using default project metadata for COBALT")
                # Default project metadata for COBALT (project already set above)
                ds.setncattr('project_principal_investigator', 'Edward Gryspeerdt')
                ds.setncattr('project_principal_investigator_email', 'e.gryspeerdt@imperial.ac.uk')
                ds.setncattr('project_principal_investigator_url', 'https://orcid.org/0000-0002-3815-4756')
                ds.setncattr('acknowledgement', 'This dataset was developed as part of the activity "Contrail Observations and Lifecycle Tracking (COBALT)", funded by\nNatural Environment Research Council (NERC) Grant NE/Z503794/1.\nIt uses instrumentation provided by the Atmospheric Measurement and Observation Facility (AMOF), part of NERC National Capability.\nUsers should acknowledge the National Centre for Atmospheric Science (NCAS) as the data provider.')
            
            # Platform information - try to get from project_instrument_info first
            platform_set = False
            if project_instrument_info and 'platform' in project_instrument_info:
                platform_info = project_instrument_info['platform']
                print("Setting platform info from project YAML instrument section")
                if 'location' in platform_info:
                    ds.setncattr('platform', platform_info['location'].lower())
                    platform_set = True
                if 'type' in platform_info:
                    ds.setncattr('platform_type', platform_info['type'])
                if 'location_keywords' in platform_info:
                    ds.setncattr('location_keywords', platform_info['location_keywords'])
                if 'deployment_mode' in platform_info:
                    ds.setncattr('deployment_mode', platform_info['deployment_mode'])
                ds.setncattr('platform_is_mobile', 'false')
            
            # Fallback to default if platform info not found in YAML
            if not platform_set:
                print("Using default platform info (cao)")
                ds.setncattr('platform', 'cao')
                ds.setncattr('platform_type', 'stationary_platform')
                ds.setncattr('location_keywords', 'Chilbolton, Hampshire, England')
                ds.setncattr('platform_is_mobile', 'false')
                ds.setncattr('deployment_mode', 'land')
            
            # Add platform_altitude - try multiple sources to ensure it's always set
            platform_altitude_set = False
            
            # First try: extract from NetCDF altitude variable
            if 'altitude' in ds.variables:
                try:
                    alt_var = ds.variables['altitude']
                    if len(alt_var) > 0:
                        platform_altitude = float(alt_var[0])
                        ds.setncattr('platform_altitude', f'{platform_altitude:.1f} m')
                        print(f"Set platform_altitude to {platform_altitude:.1f} m from altitude variable")
                        platform_altitude_set = True
                except Exception as e:
                    print(f"Warning: Could not extract platform_altitude from altitude variable: {e}")
            
            # Second try: get from YAML altitude value
            if not platform_altitude_set and project_instrument_info and 'altitude' in project_instrument_info:
                try:
                    alt_info = project_instrument_info['altitude']
                    if isinstance(alt_info, dict) and 'value' in alt_info:
                        alt_value = float(alt_info['value'])
                        alt_units = alt_info.get('units', 'm')
                        if alt_units == 'metres' or alt_units == 'm':
                            ds.setncattr('platform_altitude', f'{alt_value:.1f} m')
                            print(f"Set platform_altitude to {alt_value:.1f} m from YAML altitude")
                            platform_altitude_set = True
                except Exception as e:
                    print(f"Warning: Could not extract platform_altitude from YAML altitude: {e}")
            
            # Third try: get from YAML platform.altitude string
            if not platform_altitude_set and project_instrument_info and 'platform' in project_instrument_info:
                try:
                    platform_info = project_instrument_info['platform']
                    if isinstance(platform_info, dict) and 'altitude' in platform_info:
                        # Use the string as-is (e.g., "83 m (orthometric height above EGM2008 geoid)")
                        ds.setncattr('platform_altitude', platform_info['altitude'])
                        print(f"Set platform_altitude to '{platform_info['altitude']}' from YAML platform.altitude")
                        platform_altitude_set = True
                except Exception as e:
                    print(f"Warning: Could not extract platform_altitude from YAML platform: {e}")
            
            # Last resort: use default for Chilbolton
            if not platform_altitude_set:
                ds.setncattr('platform_altitude', '83.0 m')
                print("Set platform_altitude to default '83.0 m' for Chilbolton")
            
            # Processing information - try to get from YAML first
            if project_instrument_info and 'processing_software' in project_instrument_info:
                ps_info = project_instrument_info['processing_software']
                if 'url' in ps_info:
                    ds.setncattr('processing_software_url', ps_info['url'])
                    print(f"Set processing_software_url from YAML: {ps_info['url']}")
                else:
                    ds.setncattr('processing_software_url', 'https://github.com/longlostjames/kepler-radar-utils/releases/tag/v1.0.0')
                
                if 'version' in ps_info:
                    ds.setncattr('processing_software_version', ps_info['version'])
                    print(f"Set processing_software_version from YAML: {ps_info['version']}")
                else:
                    ds.setncattr('processing_software_version', 'v1.0.0')
            else:
                # Fallback to defaults
                ds.setncattr('processing_software_url', 'https://github.com/longlostjames/kepler-radar-utils/releases/tag/v1.0.0')
                ds.setncattr('processing_software_version', 'v1.0.0')
                print("Using default processing_software metadata")
            
            # Add instrument software info if available
            ds.setncattr('instrument_software', 'rx_client')
            ds.setncattr('instrument_software_version', '/home/vin/WORK/MBR3_XCRL_B3XC/ews/rx_client/Debug/rx_client 20161014-1951')
            
            # Add geospatial bounds - CALCULATE DYNAMICALLY
            try:
                geospatial_bounds = cfradial_get_bbox(cfradial_file)
                ds.setncattr('geospatial_bounds', geospatial_bounds)
                print(f"Added dynamic geospatial bounds: {geospatial_bounds}")
            except Exception as e:
                print(f"Warning: Could not calculate geospatial bounds: {e}")
                # Fallback to Chilbolton area
                ds.setncattr('geospatial_bounds', 'Bounding box: 51.15N -1.44E, 51.49N -1.37E')
    
            # Update timestamp
            current_time = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
            ds.setncattr('last_revised_date', current_time)
            ds.setncattr('date_created', current_time)

            # Write history attribute (required by NCAS checker; must be non-empty)
            existing_history = ''
            if hasattr(ds, 'history'):
                existing_history = ds.getncattr('history') or ''
            if not existing_history.strip():
                timestamp = datetime.utcnow().strftime('%a %b %d %H:%M:%S %Y')
                username = os.environ.get('USER', 'unknown')
                hostname = os.environ.get('HOSTNAME', socket.gethostname())
                new_entry = (
                    f"{timestamp} - user:{username} machine:{hostname} "
                    f"Convert mmclx to CF-Radial with NCAS metadata"
                )
                ds.setncattr('history', new_entry)

            # Update variable attributes for NCAS compliance
            _update_variable_attributes_for_ncas(ds, revised_northangle)
        
    except Exception as e:
        raise RuntimeError(f"Error modifying NetCDF file {cfradial_file}: {e}")
    
    print(f"Successfully added manual NCAS metadata to {cfradial_file}")

def _update_variable_attributes_for_ncas(ds, revised_northangle: float = 55.9):
    """
    Update variable attributes to be NCAS compliant.
    
    PyART automatically sets the 'type' attribute based on each variable's dtype,
    so we rely on that automatic behavior rather than explicitly setting types.
    We only update other attributes like standard_name, comment, units, etc.
    
    Args:
        ds: NetCDF4 Dataset object
        revised_northangle: North angle correction in degrees (default: 55.9 for COBALT)
    """

    def _s(val):
        """Encode str as bytes so netCDF4 writes NC_CHAR, not NC_STRING."""
        return val.encode('utf-8') if isinstance(val, str) else val

    # Update time variable
    if 'time' in ds.variables:
        time_var = ds.variables['time']
        time_var.setncattr('comment', _s(' '))
        time_var.setncattr('long_name', _s('time_since_time_reference'))
    
    # Update range variable  
    if 'range' in ds.variables:
        range_var = ds.variables['range']
        range_var.setncattr('comment', _s('Range to centre of each bin'))
        # projection_range_coordinate is not a CF standard name — use proposed_standard_name.
        if hasattr(range_var, 'standard_name'):
            range_var.delncattr('standard_name')
        range_var.setncattr('proposed_standard_name', _s('projection_range_coordinate'))
        if hasattr(range_var, 'meters_to_center_of_first_gate'):
            # Keep existing value
            pass
        else:
            range_var.setncattr('meters_to_center_of_first_gate', range_var[0])
    
    # Update azimuth variable
    if 'azimuth' in ds.variables:
        azim_var = ds.variables['azimuth']
        _default_azi_comment = 'Azimuth of antenna relative to true north'
        # Preserve any reconstruction/correction comment already written by the reader.
        # Guard against double-prefixing when this function is called more than once.
        _existing_azi = getattr(azim_var, 'comment', '')
        if _existing_azi.startswith(_default_azi_comment):
            pass  # Already prefixed — leave as is
        elif _existing_azi and _existing_azi != _default_azi_comment:
            azim_var.setncattr('comment', _s(_default_azi_comment + '. ' + _existing_azi))
        else:
            azim_var.setncattr('comment', _s(_default_azi_comment))
        # Remove any standard_name written by PyART (e.g. 'beam_azimuth_angle');
        # ray_azimuth_angle is an NCAS proposed name, not a CF standard name.
        if hasattr(azim_var, 'standard_name'):
            azim_var.delncattr('standard_name')
        azim_var.setncattr('proposed_standard_name', _s('ray_azimuth_angle'))
    
    # Update elevation variable
    if 'elevation' in ds.variables:
        elev_var = ds.variables['elevation']
        _default_elev_comment = 'Elevation of antenna relative to the horizontal plane'
        # Preserve any reconstruction/correction comment already written by the reader.
        # Guard against double-prefixing when this function is called more than once.
        _existing_elev = getattr(elev_var, 'comment', '')
        if _existing_elev.startswith(_default_elev_comment):
            pass  # Already prefixed — leave as is
        elif _existing_elev and _existing_elev != _default_elev_comment:
            elev_var.setncattr('comment', _s(_default_elev_comment + '. ' + _existing_elev))
        else:
            elev_var.setncattr('comment', _s(_default_elev_comment))
        # Remove any standard_name written by PyART (e.g. 'beam_elevation_angle');
        # ray_elevation_angle is an NCAS proposed name, not a CF standard name.
        if hasattr(elev_var, 'standard_name'):
            elev_var.delncattr('standard_name')
        elev_var.setncattr('proposed_standard_name', _s('ray_elevation_angle'))
    
    # Update latitude variable
    if 'latitude' in ds.variables:
        lat_var = ds.variables['latitude']
        # PyART automatically sets type based on variable dtype (float64 -> 'double')
    
    # Update longitude variable
    if 'longitude' in ds.variables:
        lon_var = ds.variables['longitude']
        # PyART automatically sets type based on variable dtype (float64 -> 'double')
    
    # Update altitude variable
    if 'altitude' in ds.variables:
        alt_var = ds.variables['altitude']
        alt_var.setncattr('comment', _s('Altitude of the centre of rotation of the antenna above the geoid using the WGS84 ellipsoid and EGM2008 geoid model'))
        alt_var.setncattr('units', _s('metres'))
        # PyART automatically sets type based on variable dtype (float64 -> 'double')
    
    # Update prt variable
    if 'prt' in ds.variables:
        prt_var = ds.variables['prt']
        # PyART automatically sets type based on variable dtype
    
    # Update radar_measured_transmit_power_h
    if 'radar_measured_transmit_power_h' in ds.variables:
        power_var = ds.variables['radar_measured_transmit_power_h']
        if not hasattr(power_var, 'units'):
            power_var.setncattr('units', _s('dBm'))
        # PyART automatically sets type based on variable dtype
    
    # Update sweep_mode variable
    if 'sweep_mode' in ds.variables:
        sweep_var = ds.variables['sweep_mode']
        sweep_var.setncattr('comment', _s('Options are: "sector", "coplane", "rhi", "vertical_pointing", "idle", "azimuth_surveillance", "elevation_surveillance", "sunscan", "pointing", "manual_ppi", "manual_rhi"'))
        sweep_var.setncattr('units', _s(''))
    
    # Update antenna_transition variable
    if 'antenna_transition' in ds.variables:
        ant_var = ds.variables['antenna_transition']
        ant_var.setncattr('comment', _s('1 if antenna is in transition, 0 otherwise'))
        ant_var.setncattr('units', _s(''))
    
    # Add frequency variable if missing
    if 'frequency' not in ds.variables and 'frequency' not in ds.dimensions:
        # Add frequency dimension and variable
        ds.createDimension('frequency', 1)
        freq_var = ds.createVariable('frequency', 'f4', ('frequency',))
        freq_var.setncattr('standard_name', _s('radiation_frequency'))
        freq_var.setncattr('long_name', _s('frequency of transmitted radiation'))
        freq_var.setncattr('units', _s('s-1'))
        freq_var.setncattr('meta_group', _s('instrument_parameters'))
        freq_var[:] = [35.5e9]  # 35.5 GHz for Ka-band
    
    # Add azimuth_correction variable if missing
    if 'azimuth_correction' not in ds.variables:
        azim_corr_var = ds.createVariable('azimuth_correction', 'f4')
        azim_corr_var.setncattr('long_name', _s('azimuth correction applied'))
        azim_corr_var.setncattr('units', _s('degrees'))
        azim_corr_var.setncattr('meta_group', _s('geometry_correction'))
        azim_corr_var.setncattr('comment', _s('Azimuth correction applied. North angle relative to instrument home azimuth.'))
        azim_corr_var[:] = revised_northangle

def multi_mmclx2cfrad(
   
    mmclxfiles: List[str],
    outdir: str,
    scan_name: str = 'RHI',
    gzip_flag: bool = True,
    azimuth_offset: float = 0.0,
    tracking_tag: str = '',
    campaign: str = '',
    revised_northangle: float = 55.9,
    data_version: str = "1.0.0",
    single_sweep: bool = False,
    yaml_project_file: str = None,
    yaml_instrument_file: str = None,
    logp_dir: str = None
) -> Optional[Radar]:
    """
    Convert multiple mmclx files to CF-Radial file(s).
    
    Args:
        mmclxfiles: List of mmclx file paths
        outdir: Output directory
        scan_name: Type of scan ('RHI', 'PPI', 'VPT')
        gzip_flag: Whether files are gzip compressed
        azimuth_offset: Azimuth offset correction
        tracking_tag: AMOF tracking tag
        campaign: Campaign name
        revised_northangle: North angle correction
        data_version: Data version string
        single_sweep: If True, create separate files for each sweep
        yaml_project_file: Path to project YAML file (required)
        yaml_instrument_file: Path to instrument YAML file (required)
        logp_dir: Optional path to axis-log directory (*.axis.gz). When provided,
            used to reconstruct angles for scans where mmclx angle data is missing.
        
    Returns:
        Combined radar object or None if processing fails
        
    Raises:
        ValueError: If required YAML file paths are not provided
        FileNotFoundError: If YAML files don't exist
    """
    if not mmclxfiles:
        print("No input files provided")
        return None
    
    # Validate required YAML file paths
    if yaml_project_file is None:
        raise ValueError("yaml_project_file must be provided - cannot use default paths")
    if yaml_instrument_file is None:
        raise ValueError("yaml_instrument_file must be provided - cannot use default paths")
    
    # Check that YAML files exist
    if not os.path.exists(yaml_project_file):
        raise FileNotFoundError(f"Project YAML file not found: {yaml_project_file}")
    if not os.path.exists(yaml_instrument_file):
        raise FileNotFoundError(f"Instrument YAML file not found: {yaml_instrument_file}")
    
    print(f"Processing {len(mmclxfiles)} files for {scan_name}")
    print(f"Single sweep mode: {single_sweep}")
    print(f"Data version: {data_version}")
    print(f"Using YAML files:")
    print(f"  Project: {yaml_project_file}")
    print(f"  Instrument: {yaml_instrument_file}")
    
    # Read all radar files
    radars = []
    for mmclx_file in mmclxfiles:
        try:
            radar = read_mira35_mmclx(
                mmclx_file, gzip_flag=gzip_flag,
                revised_northangle=revised_northangle,
                logp_dir=logp_dir
            )
            radars.append(radar)
        except Exception as e:
            print(f"Error reading {mmclx_file}: {e}")
            continue
    
    if not radars:
        print("No valid radar files could be read")
        return None
    
    # Determine location from project YAML file
    location = 'unknown'  # default
    try:
        with open(yaml_project_file, 'r') as f:
            projects = yaml.safe_load(f)
        
        # Find the project with the matching tracking_tag
        project = None
        for p in projects:
            if tracking_tag in p:
                project = p[tracking_tag]
                break
        
        _instr_key = next((k for k in ('chilbolton_instruments', 'ncas_instruments') if k in project), None)
        if project and _instr_key:
            # Find the radar instrument
            for instrument in project[_instr_key]:
                if any(name in instrument for name in ('ncas-mobile-ka-band-radar-1', 'reading-mobile-ka-band-radar-1')):
                    radar_name_key = next(name for name in ('reading-mobile-ka-band-radar-1', 'ncas-mobile-ka-band-radar-1') if name in instrument)
                    radar_info = instrument[radar_name_key]
                    if 'platform' in radar_info and 'location' in radar_info['platform']:
                        location = radar_info['platform']['location'].lower()
                        print(f"Found location from YAML: {location}")
                        break
        
        if location == 'unknown':
            print(f"Warning: Could not find platform location in {yaml_project_file}, using 'unknown'")
    
    except Exception as e:
        print(f"Warning: Error reading location from {yaml_project_file}: {e}")
        print("Using 'unknown' as location")
    
    scan_name_lower = scan_name.lower().replace('_', '-')
    
   
    
    if single_sweep:
        # Create separate file for each sweep
        created_files = []
        
        for i, radar in enumerate(radars):
            try:
                # Get the actual start time from the radar data
                first_time = cftime.num2pydate(radar.time['data'][0], radar.time['units'])
                dtstr = first_time.strftime('%Y%m%d-%H%M%S')
                
                # Simple filename with correct data version
                outfile = os.path.join(
                    outdir,
                    f'ncas-mobile-ka-band-radar-1_{location}_{dtstr}_{scan_name_lower}_l1_v{data_version}.nc'
                )
                
                print(f"Creating output file: {outfile}")
                
                # If multiple files would have the same name, add a sequence number
                if os.path.exists(outfile):
                    base_name = outfile.replace('.nc', '')
                    counter = 1
                    while os.path.exists(f'{base_name}_{counter:02d}.nc'):
                        counter += 1
                    outfile = f'{base_name}_{counter:02d}.nc'
                    print(f"File exists, using: {outfile}")
                
                # Write CF-Radial file
                pyart.io.write_cfradial(outfile, radar, format='NETCDF4', time_reference=True)
                
                # Add elevation and azimuth scan rates for MAN scans
                add_scan_rate_coordinates(outfile, radar)
                
                # Add time coverage attributes
                cfradial_add_time_coverage(outfile)
                
                # Add NCAS metadata
                cfradial_add_ncas_metadata(outfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version, revised_northangle)
                
                # Update history with correct information - USE "deg" INSTEAD OF DEGREE SYMBOL
                if scan_name == 'RHI':
                    angle_info = f"azimuth={radar.fixed_angle['data'][0]:.1f}deg"
                elif scan_name == 'PPI':
                    angle_info = f"elevation={radar.fixed_angle['data'][0]:.1f}deg"
                else:
                    angle_info = ""
                
                history_msg = f"{campaign.upper()}: Single-sweep {scan_name} conversion"
                if angle_info:
                    history_msg += f", {angle_info}"
                if revised_northangle:
                    history_msg += f", revised_northangle={revised_northangle}deg"
                
                update_history_attribute(outfile, history_msg)
                
                print(f"Created single-sweep file: {outfile}")
                created_files.append(outfile)
                
            except Exception as e:
                print(f"Error processing sweep {i}: {e}")
                continue
        
        print(f"Created {len(created_files)} single-sweep files")
        return radars[0] if radars else None
    
    else:
        # Original multi-sweep processing
        if len(radars) == 1:
            combined_radar = radars[0]
        else:
            try:
                # Combine multiple radar objects
                combined_radar = _combine_radars(radars, scan_name)
            except Exception as e:
                print(f"Error combining sweeps: {e}")
                return None
        
        # Generate output filename for multi-sweep file
        first_time = cftime.num2pydate(combined_radar.time['data'][0], combined_radar.time['units'])
        dtstr = first_time.strftime('%Y%m%d-%H%M%S')
        
        outfile = os.path.join(
            outdir,
            f'ncas-mobile-ka-band-radar-1_{location}_{dtstr}_{scan_name_lower}_l1_v{data_version}.nc'
        )
        
        print(f"Creating multi-sweep output file: {outfile}")
        
        # Write CF-Radial file
        pyart.io.write_cfradial(outfile, combined_radar, format='NETCDF4', time_reference=True)
        
        # Add elevation and azimuth scan rates for MAN scans
        add_scan_rate_coordinates(outfile, combined_radar)
        
        # Add time coverage attributes
        cfradial_add_time_coverage(outfile)
        
        # Add NCAS metadata
        cfradial_add_ncas_metadata(outfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version, revised_northangle)
        
        # Update history - USE "deg" INSTEAD OF DEGREE SYMBOL
        history_msg = f"{campaign.upper()}: Multi-sweep {scan_name} conversion"
        if revised_northangle:
            history_msg += f", revised_northangle={revised_northangle}deg"
        history_msg += f", {len(radars)} sweeps"
        
        update_history_attribute(outfile, history_msg)
        
        print(f"Created multi-sweep file: {outfile}")
        return combined_radar

def cfradial_add_time_coverage(cfradial_file: str) -> None:
    """
    Add time_coverage_start and time_coverage_end attributes to CF-Radial file.
    
    These attributes are required for NCAS compliance and help with data discovery.
    They are calculated from the actual time data in the file.
    
    Args:
        cfradial_file: Path to CF-Radial file to modify
    """
    try:
        with nc4.Dataset(cfradial_file, 'a') as ds:
            # Check if time variable exists
            if 'time' not in ds.variables:
                print(f"Warning: No time variable found in {cfradial_file}")
                return
            
            time_var = ds.variables['time']
            
            # Get time units for conversion
            if hasattr(time_var, 'units'):
                time_units = time_var.units
            else:
                print(f"Warning: No time units found in {cfradial_file}")
                return
            
            # Convert first and last time values to datetime
            try:
                start_time = cftime.num2pydate(time_var[0], time_units)
                end_time = cftime.num2pydate(time_var[-1], time_units)
                
                # Format as ISO 8601 strings
                time_coverage_start = start_time.strftime('%Y-%m-%dT%H:%M:%SZ')
                time_coverage_end = end_time.strftime('%Y-%m-%dT%H:%M:%SZ')
                
                # Set the attributes
                ds.setncattr('time_coverage_start', time_coverage_start)
                ds.setncattr('time_coverage_end', time_coverage_end)
                
                print(f"Added time coverage: {time_coverage_start} to {time_coverage_end}")
                
            except Exception as e:
                print(f"Error converting time values: {e}")
                return
                
    except Exception as e:
        print(f"Error adding time coverage to {cfradial_file}: {e}")

def _safe_parse_time(nc_dataset, time_index=0):
    """
    Safely parse time from netCDF dataset, handling microseconds and various time unit formats.
    
    Args:
        nc_dataset: Open netCDF4 Dataset
        time_index: Index of time value to parse (default: 0)

        
    Returns:
        datetime object or None if parsing fails
    """
    try:
        # First try the standard approach with microseconds (as in your old code)
        if 'time' in nc_dataset.variables and 'microsec' in nc_dataset.variables:
            time_var = nc_dataset.variables['time']
            microsec_var = nc_dataset.variables['microsec']
            
            # Parse the base time
            dtime = cftime.num2pydate(time_var[time_index], 'seconds since 1970-01-01 00:00:00')
            
            # Add microseconds
            dtime = dtime.replace(microsecond=int(microsec_var[time_index]))
            
            return dtime
            
    except (ValueError, OSError, IndexError) as e:
        print(f"Microsecond time parsing failed: {e}")
        
    try:
        # Fallback: Try standard time parsing without microseconds
        if 'time' in nc_dataset.variables:
            time_var = nc_dataset.variables['time']
            if hasattr(time_var, 'units'):
                # Clean up common problematic units
                units = time_var.units
                units = units.replace('seconds since 1970-1-1 0:0:0', 'seconds since 1970-01-01 00:00:00')
                units = units.replace('seconds since 1970-01-01 0:00:00', 'seconds since 1970-01-01 00:00:00')
                return cftime.num2pydate(time_var[time_index], units)
            else:
                return cftime.num2pydate(time_var[time_index], 'seconds since 1970-01-01 00:00:00')
    except (ValueError, OSError, IndexError) as e:
        print(f"Standard time parsing failed: {e}")
        
    try:
        # Method 3: Use raw time value and assume epoch
        if 'time' in nc_dataset.variables:
            time_val = nc_dataset.variables['time'][time_index]
            return datetime.datetime.utcfromtimestamp(time_val)
    except:
        pass
    
    # Method 4: Return a default time based on current processing
    print("Warning: Could not parse time from file, using current time")
    return datetime.datetime.utcnow()

def _combine_radars(radar_list: List[Radar], scan_name: str) -> Radar:
    """
    Combine multiple radar objects into a single radar object.
    
    Args:
        radar_list: List of PyART Radar objects to combine
        scan_name: Type of scan being combined
        
    Returns:
        Combined PyART Radar object
    """
    # Use scan-specific combination strategies
    scan_name_lower = scan_name.lower()
    
    if scan_name_lower == 'ppi':
        return _combine_ppi_sweeps(radar_list)
    elif scan_name_lower in ['vpt', 'vertical_pointing']:
        return _combine_vpt_sweeps(radar_list)
    elif scan_name_lower in ['rhi', 'manual_rhi']:
        return _combine_rhi_sweeps(radar_list)
    else:
        # Default to manual combination
        print(f"Using manual combination for scan type: {scan_name}")
        return _manual_combine_radars(radar_list, scan_name)

def _manual_combine_radars(radar_list: List[Radar], scan_name: str) -> Radar:
    """
    Manually combine radar objects by concatenating arrays.
    
    Args:
        radar_list: List of PyART Radar objects to combine
        scan_name: Type of scan being combined
        
    Returns:
        Combined PyART Radar object
    """
    base_radar = radar_list[0]
    
    # Calculate total dimensions
    total_nrays = sum(r.nrays for r in radar_list)
    nsweeps = len(radar_list)
    
    print(f"Combining {nsweeps} sweeps with {total_nrays} total rays")
    
    # Create new time array
    time_data = []
    for radar in radar_list:
        time_data.extend(radar.time['data'])
    
    combined_time = base_radar.time.copy()
    combined_time['data'] = np.array(time_data)
    
    # Create new angle arrays - explicitly use float32
    azimuth_data = []
    elevation_data = []
    for radar in radar_list:
        azimuth_data.extend(radar.azimuth['data'])
        elevation_data.extend(radar.elevation['data'])
    
    combined_azimuth = base_radar.azimuth.copy()
    combined_azimuth['data'] = np.array(azimuth_data, dtype=np.float32)
    
    combined_elevation = base_radar.elevation.copy()
    combined_elevation['data'] = np.array(elevation_data, dtype=np.float32)
    
    # Combine fields
    combined_fields = {}
    for field_name in base_radar.fields:
        field_data = []
        for radar in radar_list:
            if field_name in radar.fields:
                field_data.append(radar.fields[field_name]['data'])
        
        if field_data:
            combined_field = base_radar.fields[field_name].copy()
            combined_field['data'] = np.concatenate(field_data, axis=0)
            combined_fields[field_name] = combined_field
    
    # Create sweep indexing
    sweep_start_ray_index = {'data': np.zeros(nsweeps, dtype='int32')}
    sweep_end_ray_index = {'data': np.zeros(nsweeps, dtype='int32')}
    sweep_number = {'data': np.arange(nsweeps, dtype='int32')}
    fixed_angle = {'data': np.zeros(nsweeps, dtype='float32')}
    sweep_mode = {'data': np.array([scan_name.lower()] * nsweeps)}
    
    current_ray = 0
    for i, radar in enumerate(radar_list):
        sweep_start_ray_index['data'][i] = current_ray
        sweep_end_ray_index['data'][i] = current_ray + radar.nrays - 1
        fixed_angle['data'][i] = radar.fixed_angle['data'][0]
        current_ray += radar.nrays
    
    # Combine antenna_transition if present on any radar
    combined_antenna_transition = None
    if any(r.antenna_transition is not None for r in radar_list):
        at_data = []
        for radar in radar_list:
            if radar.antenna_transition is not None:
                at_data.append(radar.antenna_transition['data'])
            else:
                # Radar has no antenna_transition: treat all its rays as non-transition
                at_data.append(np.zeros(radar.nrays, dtype='int8'))
        combined_antenna_transition = base_radar.antenna_transition.copy() if base_radar.antenna_transition is not None else {}
        combined_antenna_transition['data'] = np.concatenate(at_data)

    # Create combined radar object
    combined_radar = Radar(
        combined_time, base_radar.range, combined_fields,
        base_radar.metadata, scan_name.lower(),
        base_radar.latitude, base_radar.longitude, base_radar.altitude,
        sweep_number, sweep_mode, fixed_angle,
        sweep_start_ray_index, sweep_end_ray_index,
        combined_azimuth, combined_elevation,
        antenna_transition=combined_antenna_transition,
        instrument_parameters=base_radar.instrument_parameters
    )
    
    return combined_radar

def split_monotonic_sequence(azimuth_data: np.ndarray, tolerance: float = 5.0) -> List[Tuple[int, int]]:
    """
    Split azimuth data into monotonic sequences for PPI processing.
    
    This function identifies breaks in monotonic azimuth sequences that indicate
    separate PPI sweeps, handling 360-degree wraparound.
    
    Args:
        azimuth_data: Array of azimuth angles in degrees
        tolerance: Maximum allowed jump in degrees before considering a new sequence
        
    Returns:
        List of tuples (start_index, end_index) for each monotonic sequence
        
    Example:
        >>> azimuth = np.array([0, 10, 20, 30, 350, 0, 10, 20])
        >>> sequences = split_monotonic_sequence(azimuth, tolerance=30.0)
        >>> print(sequences)
        [(0, 3), (4, 7)]
    """
    if len(azimuth_data) < 2:
        return [(0, len(azimuth_data) - 1)]
    
    sequences = []
    start_idx = 0
    
    for i in range(1, len(azimuth_data)):
        # Calculate azimuth difference, handling 360-degree wraparound
        diff = azimuth_data[i] - azimuth_data[i-1]
        
        # Handle wraparound cases
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360
        
        # Check if this represents a break in the monotonic sequence
        if abs(diff) > tolerance:
            # End the current sequence
            sequences.append((start_idx, i - 1))
            start_idx = i
    
    # Add the final sequence
    sequences.append((start_idx, len(azimuth_data) - 1))
    
    print(f"Found {len(sequences)} monotonic sequences in azimuth data")
    for i, (start, end) in enumerate(sequences):
        print(f"  Sequence {i+1}: indices {start}-{end} "
              f"(az: {azimuth_data[start]:.1f}° to {azimuth_data[end]:.1f}°)")
    
    return sequences

def _combine_ppi_sweeps(radar_list: List[Radar]) -> Radar:
    """
    Combine PPI radar sweeps, handling monotonic azimuth sequences.
    
    Args:
        radar_list: List of PyART Radar objects to combine
        
    Returns:
        Combined PyART Radar object with proper sweep indexing
    """
    if len(radar_list) == 1:
        return radar_list[0]
    
    base_radar = radar_list[0]
    
    # Collect all azimuth data to identify sweep boundaries
    all_azimuths = []
    all_times = []
    all_elevations = []
    ray_counts = []
    
    for radar in radar_list:
        all_azimuths.extend(radar.azimuth['data'])
        all_elevations.extend(radar.elevation['data'])
        all_times.extend(radar.time['data'])
        ray_counts.append(radar.nrays)
    
    all_azimuths = np.array(all_azimuths, dtype=np.float32)
    all_elevations = np.array(all_elevations, dtype=np.float32)
    all_times = np.array(all_times)
    
    # Identify monotonic azimuth sequences (individual sweeps)
    sequences = split_monotonic_sequence(all_azimuths, tolerance=10.0)
    
    nsweeps = len(sequences)
    total_nrays = len(all_azimuths)
    
    print(f"Combining {len(radar_list)} files into {nsweeps} PPI sweeps with {total_nrays} total rays")
    
    # Create combined time, azimuth, elevation arrays
    combined_time = base_radar.time.copy()
    combined_time['data'] = all_times
    
    combined_azimuth = base_radar.azimuth.copy() 
    combined_azimuth['data'] = all_azimuths
    
    combined_elevation = base_radar.elevation.copy()
    combined_elevation['data'] = all_elevations
    
    # Combine fields from all radars
    combined_fields = {}
    for field_name in base_radar.fields:
        field_data = []
        for radar in radar_list:
            if field_name in radar.fields:
                field_data.append(radar.fields[field_name]['data'])
        
        if field_data:
            combined_field = base_radar.fields[field_name].copy()
            combined_field['data'] = np.concatenate(field_data, axis=0)
            combined_fields[field_name] = combined_field
    
    # Create sweep indexing based on monotonic sequences
    sweep_start_ray_index = {'data': np.array([seq[0] for seq in sequences], dtype='int32')}
    sweep_end_ray_index = {'data': np.array([seq[1] for seq in sequences], dtype='int32')}
    sweep_number = {'data': np.arange(nsweeps, dtype='int32')}
    
    # Calculate fixed angle for each sweep (median elevation)
    fixed_angles = []
    for start_idx, end_idx in sequences:
        sweep_elevations = all_elevations[start_idx:end_idx+1]
        fixed_angles.append(np.median(sweep_elevations))
    
    fixed_angle = {'data': np.array(fixed_angles, dtype='float32')}
    sweep_mode = {'data': np.array(['ppi'] * nsweeps)}
    
    # Combine antenna_transition if present on any radar
    combined_ppi_antenna_transition = None
    if any(r.antenna_transition is not None for r in radar_list):
        at_data = []
        for radar in radar_list:
            if radar.antenna_transition is not None:
                at_data.append(radar.antenna_transition['data'])
            else:
                at_data.append(np.zeros(radar.nrays, dtype='int8'))
        combined_ppi_antenna_transition = base_radar.antenna_transition.copy() if base_radar.antenna_transition is not None else {}
        combined_ppi_antenna_transition['data'] = np.concatenate(at_data)

    # Create combined radar object
    combined_radar = Radar(
        combined_time, base_radar.range, combined_fields,
        base_radar.metadata, 'ppi',
        base_radar.latitude, base_radar.longitude, base_radar.altitude,
        sweep_number, sweep_mode, fixed_angle,
        sweep_start_ray_index, sweep_end_ray_index,
        combined_azimuth, combined_elevation,
        antenna_transition=combined_ppi_antenna_transition,
        instrument_parameters=base_radar.instrument_parameters
    )
    
    return combined_radar

def _combine_rhi_sweeps(radar_list: List[Radar]) -> Radar:
    """
    Combine RHI radar sweeps.
    
    Args:
        radar_list: List of PyART Radar objects to combine
        
    Returns:
        Combined PyART Radar object
    """
    # For RHI, use the general radar combination since each file is typically one sweep
    return _manual_combine_radars(radar_list, 'rhi')

def _combine_vpt_sweeps(radar_list: List[Radar]) -> Radar:
    """
    Combine VPT (Vertical Pointing) radar sweeps.
    
    VPT scans are typically continuous vertical pointing measurements
    that should be combined chronologically into a multi-sweep time series.
    
    Args:
        radar_list: List of PyART Radar objects to combine
        
    Returns:
        Combined PyART Radar object with chronological sweep organization
    """
    if len(radar_list) == 1:
        return radar_list[0]
    
    # Sort radars by time to ensure chronological order
    radar_list = sorted(radar_list, key=lambda r: r.time['data'][0])
    
    print(f"Combining {len(radar_list)} VPT sweeps chronologically")
    
    # Use manual combination for VPT since each file is typically one sweep
    return _manual_combine_radars(radar_list, 'vertical_pointing')

def find_mmclx_vad_files(
    start_time: str,
    end_time: str,
    elev_min: float,
    elev_max: float,
    inpath: str,
    gzip_flag: bool = False
) -> List[str]:
    """
    Find VAD (Volume Antenna Diagram) mmclx files within time and elevation ranges.
    
    Args:
        start_time: Start time as 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'
        end_time: End time as 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DDTHH:MM:SSZ'
        elev_min: Minimum elevation angle
        elev_max: Maximum elevation angle
        inpath: Directory path to search
        gzip_flag: Whether files are gzip compressed
        
    Returns:
        List of matching VAD file paths
    """
    # Handle both time formats
    def parse_time_string(time_str):
        """Parse time string in either format"""
        try:
            # Try ISO format first
            return datetime.datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ')
        except ValueError:
            try:
                # Try space-separated format
                return datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
            except ValueError:
                raise ValueError(f"Time string '{time_str}' does not match expected formats")
    
    # Convert the input times to datetime objects
    start_datetime = parse_time_string(start_time)
    end_datetime = parse_time_string(end_time)
    
    print(f"Searching for VAD files from {start_time} to {end_time}")
    print(f"Elevation range: {elev_min}° to {elev_max}°")
    
    # Get the date string from start time
    date_str = start_datetime.strftime('%Y%m%d')
    
    # Try date subdirectory first (standard campaigns); fall back to inpath directly (WOEST)
    inpath_obj = Path(inpath)
    date_subdir = inpath_obj / date_str
    if date_subdir.exists():
        search_path = date_subdir
        print(f"Looking in date directory: {search_path}")
    elif inpath_obj.exists():
        search_path = inpath_obj
        print(f"Using existing path directly (no date subdir): {search_path}")
    else:
        print(f"Path does not exist: {inpath_obj} (with or without date {date_str})")
        return []
    
    if not search_path.exists():
        print(f"Date directory does not exist: {search_path}")
        return []
    
    # Initialize matching files list
    matching_files = []
    
    # Look for VAD files (usually PPI scans at high elevation)
    if gzip_flag:
        vad_pattern = "*ppi*.mmclx.gz"  # VAD files are typically PPI scans
    else:
        vad_pattern = "*ppi*.mmclx"
    
    candidate_files = list(search_path.glob(vad_pattern))
    print(f"Found {len(candidate_files)} candidate VAD files")
    
    # Check each file's timestamp and elevation
    for file_path in candidate_files:
        try:
            if gzip_flag:
                with gzip.open(file_path, 'rb') as gz:
                    with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                        if 'time' in nc.dimensions and len(nc.dimensions['time']) > 0:
                            file_time = _safe_parse_time(nc, 0)
                            if file_time is None:
                                continue
                                
                            # Check elevation angle
                            if 'elv' in nc.variables:
                                elev = nc['elv'][0]  # First elevation angle
                                
                                if (start_datetime <= file_time <= end_datetime and 
                                    elev_min <= elev <= elev_max):
                                    matching_files.append(str(file_path))
                                    print(f"Added VAD file: {file_path.name} (elev: {elev:.1f}°)")
            else:
                with nc4.Dataset(file_path, 'r') as nc:
                    if 'time' in nc.dimensions and len(nc.dimensions['time']) > 0:
                        file_time = _safe_parse_time(nc, 0)
                        if file_time is None:
                            continue
                            
                        # Check elevation angle
                        if 'elv' in nc.variables:
                            elev = nc['elv'][0]  # First elevation angle
                            
                            if (start_datetime <= file_time <= end_datetime and 
                                elev_min <= elev <= elev_max):
                                matching_files.append(str(file_path))
                                print(f"Added VAD file: {file_path.name} (elev: {elev:.1f}°)")
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            continue
    
    print(f"Found {len(matching_files)} matching VAD files")
    return sorted(matching_files)





def export_sweep_metadata_to_csv(
    cfradial_file: str,
    sweep_csv: str = None,
    antenna_transition_csv: str = None
) -> None:
    """
    Export sweep metadata and antenna_transition flags from CF-Radial file to CSV.
    
    Creates two CSV files for manual editing:
    1. Sweep metadata (sweep_number, start_idx, end_idx, start_time, end_time, n_rays, sweep_mode, phase, fixed_angle, start_az, end_az, start_el, end_el)
    2. Antenna transition flags (ray_index, time, azimuth, elevation, antenna_transition)
    
    The CSV includes both:
    - sweep_mode: CF-Radial standard (manual_rhi, vertical_pointing, etc.) - preserved in NetCDF
    - phase: Campaign-specific descriptions (upward_rhi, downward_rhi, tracking, dwell, etc.)
      For PICASSO: from phase_sequence metadata
      For other campaigns: defaults to sweep_mode if no phase_sequence exists
    
    Args:
        cfradial_file: Path to CF-Radial NetCDF file
        sweep_csv: Output path for sweep metadata CSV (default: <file>_sweeps.csv)
        antenna_transition_csv: Output path for antenna transition CSV (default: <file>_antenna_transition.csv)
        
    Example:
        >>> export_sweep_metadata_to_csv('myfile.nc')
        Exported sweep metadata to myfile_sweeps.csv (77 sweeps)
        Exported antenna transition flags to myfile_antenna_transition.csv (1234 rays, 56 transitions)
    """
    import csv
    import netCDF4 as nc4
    from pathlib import Path
    
    # Generate default output paths if not provided
    if sweep_csv is None:
        sweep_csv = str(Path(cfradial_file).with_suffix('')) + '_sweeps.csv'
    if antenna_transition_csv is None:
        antenna_transition_csv = str(Path(cfradial_file).with_suffix('')) + '_antenna_transition.csv'
    
    with nc4.Dataset(cfradial_file, 'r') as ds:
        # Export sweep metadata
        nsweeps = len(ds.dimensions['sweep'])
        
        # Read time variable and convert to datetime
        time_var = ds.variables['time']
        time_data = time_var[:]
        time_units = time_var.units
        
        # Convert time to datetime objects
        import cftime
        time_datetimes = cftime.num2pydate(time_data, time_units)
        
        # Read azimuth and elevation data
        azimuth_data = ds.variables['azimuth'][:]
        elevation_data = ds.variables['elevation'][:]
        
        # Read sweep_mode variable (CF-Radial standard)
        sweep_mode_var = ds.variables['sweep_mode']
        if sweep_mode_var.dtype.char == 'S':
            # Character array - use netCDF4's chartostring conversion
            sweep_modes = nc4.chartostring(sweep_mode_var[:])
            sweep_modes = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in sweep_modes]
        else:
            # Already strings or other format
            sweep_modes = [str(sweep_mode_var[i]).strip() for i in range(nsweeps)]
        
        # Read phase_sequence metadata (campaign-specific phase descriptions)
        phase_list = []
        if 'phase_sequence' in ds.ncattrs():
            phase_sequence = ds.getncattr('phase_sequence')
            phase_list = [p.strip() for p in phase_sequence.split(',')]
        
        # If no phase_sequence or mismatch, use sweep_mode as phase
        if len(phase_list) != nsweeps:
            if phase_list:
                print(f"  INFO: phase_sequence has {len(phase_list)} phases but file has {nsweeps} sweeps, using sweep_mode for phase column")
            phase_list = sweep_modes.copy()
        
        # Get total number of rays for bounds checking
        n_rays_total = len(time_data)
        
        with open(sweep_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['sweep_number', 'start_idx', 'end_idx', 'start_time', 'end_time', 'n_rays', 'sweep_mode', 'phase', 'fixed_angle', 'start_az', 'end_az', 'start_el', 'end_el'])
            
            for i in range(nsweeps):
                sweep_num = int(ds.variables['sweep_number'][i])
                start_idx = int(ds.variables['sweep_start_ray_index'][i])
                end_idx = int(ds.variables['sweep_end_ray_index'][i])
                
                # Validate indices are within bounds
                if start_idx < 0 or start_idx >= n_rays_total:
                    print(f"  WARNING: Sweep {i} start_idx {start_idx} out of bounds (0-{n_rays_total-1}), skipping")
                    continue
                if end_idx < 0 or end_idx >= n_rays_total:
                    print(f"  WARNING: Sweep {i} end_idx {end_idx} out of bounds (0-{n_rays_total-1}), clamping to {n_rays_total-1}")
                    end_idx = n_rays_total - 1
                
                n_rays = end_idx - start_idx + 1
                
                # Get timestamps for start and end rays
                start_time = time_datetimes[start_idx].strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]  # milliseconds
                end_time = time_datetimes[end_idx].strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
                
                # Get azimuth and elevation at start and end of sweep
                start_az = float(azimuth_data[start_idx])
                end_az = float(azimuth_data[end_idx])
                start_el = float(elevation_data[start_idx])
                end_el = float(elevation_data[end_idx])
                
                # Get sweep_mode and phase for this sweep
                sweep_mode = sweep_modes[i].strip()
                phase = phase_list[i].strip() if i < len(phase_list) else sweep_mode
                
                fixed_angle = float(ds.variables['fixed_angle'][i])
                
                writer.writerow([sweep_num, start_idx, end_idx, start_time, end_time, n_rays, sweep_mode, phase, f'{fixed_angle:.2f}', f'{start_az:.2f}', f'{end_az:.2f}', f'{start_el:.2f}', f'{end_el:.2f}'])
        
        print(f"Exported sweep metadata to {sweep_csv} ({nsweeps} sweeps)")
        
        # Export antenna_transition flags
        if 'antenna_transition' in ds.variables:
            antenna_transition = ds.variables['antenna_transition'][:]
        else:
            # Variable absent (e.g. file pre-dates pipeline support): derive from elevation.
            # Rays below 85° are considered transitioning (not yet at vertical pointing).
            antenna_transition = np.where(elevation_data < 85.0, 1, 0).astype('int8')
            print(f"  No antenna_transition variable in file – deriving from elevation (< 85° → 1)")

        n_rays = len(antenna_transition)
        n_transitions = int(np.sum(antenna_transition == 1))

        with open(antenna_transition_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['ray_index', 'time', 'azimuth', 'elevation', 'antenna_transition'])

            for i in range(n_rays):
                ray_time = time_datetimes[i].strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]  # milliseconds
                azimuth = float(azimuth_data[i])
                elevation = float(elevation_data[i])
                writer.writerow([i, ray_time, f'{azimuth:.2f}', f'{elevation:.2f}', int(antenna_transition[i])])

        print(f"Exported antenna transition flags to {antenna_transition_csv} ({n_rays} rays, {n_transitions} transitions)")


def import_sweep_metadata_from_csv(
    cfradial_file: str,
    sweep_csv: str = None,
    antenna_transition_csv: str = None,
    backup: bool = True,
    write_phase_sequence: bool = True
) -> None:
    """
    Import edited sweep metadata and antenna_transition flags from CSV into CF-Radial file.
    
    Updates the NetCDF file with edited values from CSV. Creates a backup by default.
    
    IMPORTANT: This function can only modify existing sweep boundaries and antenna_transition
    flags. It CANNOT change the number of sweeps (e.g., splitting one sweep into two).
    If you need to change the number of sweeps, you must regenerate the file from raw data
    with the modified sweep structure.
    
    Args:
        cfradial_file: Path to CF-Radial NetCDF file to modify
        sweep_csv: Input path for sweep metadata CSV (default: <file>_sweeps.csv)
        antenna_transition_csv: Input path for antenna transition CSV (default: <file>_antenna_transition.csv)
        backup: Create backup file before modifying (default: True)
        write_phase_sequence: Write phase_sequence global attribute from the phase column
            of the sweep CSV (default: True). Set to False for campaigns that do not use
            phase_sequence metadata (e.g. COBALT).
        
    Example:
        >>> import_sweep_metadata_from_csv('myfile.nc')
        Created backup: myfile.nc.bak
        Updated sweep metadata (77 sweeps)
        Updated antenna transition flags (1234 rays, 62 transitions)
        
    Note:
        To split or merge sweeps (changing the total number), regenerate the file rather
        than using this import function.
    """
    import csv
    import netCDF4 as nc4
    from pathlib import Path
    import shutil
    
    # Generate default input paths if not provided
    if sweep_csv is None:
        sweep_csv = str(Path(cfradial_file).with_suffix('')) + '_sweeps.csv'
    if antenna_transition_csv is None:
        antenna_transition_csv = str(Path(cfradial_file).with_suffix('')) + '_antenna_transition.csv'
    
    # Create backup if requested
    if backup:
        backup_file = cfradial_file + '.bak'
        shutil.copy2(cfradial_file, backup_file)
        print(f"Created backup: {backup_file}")
    
    # Read sweep metadata from CSV
    sweep_data = []
    phases = []
    sweep_modes = []
    if Path(sweep_csv).exists():
        with open(sweep_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                sweep_data.append({
                    'sweep_number': int(row['sweep_number']),
                    'start_idx': int(row['start_idx']),
                    'end_idx': int(row['end_idx']),
                    'fixed_angle': float(row['fixed_angle'])
                })
                # Read phase if it exists in CSV (campaign-specific)
                if 'phase' in row:
                    phases.append(row['phase'].strip())
                # Read sweep_mode if it exists in CSV (CF-Radial standard)
                if 'sweep_mode' in row:
                    sweep_modes.append(row['sweep_mode'].strip())
    else:
        print(f"Sweep CSV not found: {sweep_csv}")
    
    # Read antenna_transition from CSV
    antenna_transition_data = {}
    if Path(antenna_transition_csv).exists():
        with open(antenna_transition_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                ray_idx = int(row['ray_index'])
                antenna_transition_data[ray_idx] = int(row['antenna_transition'])
    else:
        print(f"Antenna transition CSV not found: {antenna_transition_csv}")
    
    # Update NetCDF file
    with nc4.Dataset(cfradial_file, 'a') as ds:
        # Update sweep metadata
        sweep_updated = False
        if sweep_data:
            nsweeps_csv = len(sweep_data)
            nsweeps_nc = len(ds.dimensions['sweep'])
            
            if nsweeps_csv != nsweeps_nc:
                print(f"\nERROR: Cannot change the number of sweeps!")
                print(f"  Current file:  {nsweeps_nc} sweeps")
                print(f"  CSV file:      {nsweeps_csv} sweeps")
                print(f"  Difference:    {nsweeps_csv - nsweeps_nc:+d} sweeps")
                print(f"\nNetCDF dimension sizes are fixed and cannot be changed after file creation.")
                print(f"\nOptions:")
                print(f"  1. Use rebuild_file_from_csv() instead to recreate the file:")
                print(f"     >>> from kepler_utils import rebuild_file_from_csv")
                print(f"     >>> rebuild_file_from_csv('{cfradial_file}')")
                print(f"  2. Revert your CSV to have {nsweeps_nc} sweeps")
                print(f"  3. Regenerate from raw data with modified phase detection parameters")
                return
            
            for i, sweep in enumerate(sweep_data):
                ds.variables['sweep_number'][i] = sweep['sweep_number']
                ds.variables['sweep_start_ray_index'][i] = sweep['start_idx']
                ds.variables['sweep_end_ray_index'][i] = sweep['end_idx']
                ds.variables['fixed_angle'][i] = sweep['fixed_angle']
            
            # Update sweep_mode variable if values were provided in CSV
            if sweep_modes and len(sweep_modes) == nsweeps_csv and 'sweep_mode' in ds.variables:
                sweep_mode_var = ds.variables['sweep_mode']
                # Handle character array (old-style NetCDF3 strings)
                if sweep_mode_var.dtype.char == 'S':
                    # Convert list of strings to character array
                    sweep_mode_array = np.array(sweep_modes, dtype='S32')
                    char_array = nc4.stringtochar(sweep_mode_array)
                    sweep_mode_var[:] = char_array
                else:
                    # Direct string assignment for NetCDF4 string types
                    for i, mode in enumerate(sweep_modes):
                        sweep_mode_var[i] = mode
                print(f"Updated sweep_mode variable for {nsweeps_csv} sweeps")
            
            # Update phase_sequence metadata if phases were provided and requested
            if write_phase_sequence and phases and len(phases) == nsweeps_csv:
                phase_sequence_str = ', '.join(phases)
                ds.setncattr('phase_sequence', phase_sequence_str)
                print(f"Updated phase_sequence metadata: {phase_sequence_str}")
            
            sweep_updated = True
            print(f"Updated sweep metadata ({nsweeps_csv} sweeps)")
        
        # Update antenna_transition
        transition_updated = False
        if antenna_transition_data:
            n_rays_total = len(ds.variables['time'])

            if 'antenna_transition' not in ds.variables:
                # Create the variable – this file predates pipeline support
                at_var = ds.createVariable('antenna_transition', 'i1', ('time',),
                                           fill_value=np.int8(0))
                at_var.long_name = "antenna is in transition between sweeps"
                at_var.comment = (
                    "1 if elevation < 85° (antenna not yet at vertical pointing position), "
                    "0 otherwise"
                )
                at_var[:] = np.zeros(n_rays_total, dtype='int8')
                print(f"  Created antenna_transition variable ({n_rays_total} rays)")

            antenna_transition_var = ds.variables['antenna_transition']
            n_rays = len(antenna_transition_var)

            # Write values from CSV
            for ray_idx, value in antenna_transition_data.items():
                if 0 <= ray_idx < n_rays:
                    antenna_transition_var[ray_idx] = value

            n_transitions = int(np.sum(antenna_transition_var[:] == 1))
            transition_updated = True
            print(f"Updated antenna transition flags ({n_rays} rays, {n_transitions} transitions)")
        
        # Update last_revised_date and history
        if sweep_updated or transition_updated:
            import datetime
            current_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
            ds.setncattr('last_revised_date', current_time)

            history_parts = []
            if sweep_updated:
                history_parts.append("sweep metadata")
            if transition_updated:
                history_parts.append("antenna_transition flags")
            history_msg = f"Applied manual edits from CSV ({', '.join(history_parts)})"

            timestamp = datetime.datetime.utcnow().strftime('%a %b %d %H:%M:%S %Y')
            username = os.environ.get('USER', 'unknown')
            hostname = os.environ.get('HOSTNAME', 'unknown')
            new_entry = f"{timestamp} - user:{username} machine:{hostname} {history_msg}"

            current_history = ""
            if 'history' in ds.ncattrs():
                current_history = ds.getncattr('history')
            updated_history = f"{new_entry}\n{current_history}" if current_history else new_entry
            ds.setncattr('history', updated_history)

            print(f"Updated last_revised_date and history")

    print(f"Successfully updated {cfradial_file}")


def export_sweep_metadata_for_date(
    datestr: str,
    data_dir: str,
    pattern: str = '*.nc',
    overwrite: bool = False
) -> None:
    """
    Export sweep metadata and antenna_transition flags from all NetCDF files for a given date.
    
    Finds all NetCDF files matching the pattern for the specified date and exports
    their sweep metadata and antenna transition flags to CSV files for manual editing.
    
    Args:
        datestr: Date string in YYYYMMDD format
        data_dir: Base directory containing date subdirectories (e.g., /path/to/output)
        pattern: Filename pattern to match (default: '*.nc')
        overwrite: Overwrite existing CSV files (default: False)
        
    Example:
        >>> export_sweep_metadata_for_date('20171213', '/path/to/picasso/L1_v1.0.0', pattern='*_man_*.nc')
        Processing date: 20171213
        Found 5 files to process
        
        Processing: ncas-mobile-ka-band-radar-1_cao_20171213-212527_man_l1_v1.0.0.nc
        Exported sweep metadata to ...sweeps.csv (12 sweeps)
        Exported antenna transition flags to ...antenna_transition.csv (1234 rays, 56 transitions)
        ...
        
        Completed: 5 files processed, 0 skipped
    """
    import glob
    import os
    from pathlib import Path
    
    # Construct the date directory path
    date_path = os.path.join(data_dir, datestr)
    
    if not os.path.exists(date_path):
        print(f"ERROR: Date directory not found: {date_path}")
        return
    
    # Find all matching files
    search_pattern = os.path.join(date_path, pattern)
    nc_files = sorted(glob.glob(search_pattern))
    
    if not nc_files:
        print(f"No files found matching pattern: {search_pattern}")
        return
    
    print(f"Processing date: {datestr}")
    print(f"Found {len(nc_files)} file(s) to process\n")
    
    processed = 0
    skipped = 0
    
    for nc_file in nc_files:
        print(f"Processing: {os.path.basename(nc_file)}")
        
        # Generate expected CSV filenames
        base_name = str(Path(nc_file).with_suffix(''))
        sweep_csv = base_name + '_sweeps.csv'
        antenna_csv = base_name + '_antenna_transition.csv'
        
        # Check if CSV files already exist
        if not overwrite and (os.path.exists(sweep_csv) or os.path.exists(antenna_csv)):
            print(f"  CSV files already exist - skipping (use --overwrite to replace)")
            skipped += 1
            continue
        
        try:
            export_sweep_metadata_to_csv(
                cfradial_file=nc_file,
                sweep_csv=sweep_csv,
                antenna_transition_csv=antenna_csv
            )
            processed += 1
        except Exception as e:
            print(f"  ERROR: Failed to export {os.path.basename(nc_file)}: {e}")
            import traceback
            traceback.print_exc()
            skipped += 1
        
        print()  # Blank line between files
    
    print("="*70)
    print(f"Completed: {processed} file(s) processed, {skipped} skipped")
    print("="*70)


def apply_csv_edits_for_date(
    datestr: str,
    data_dir: str,
    csv_dir: str = None,
    backup: bool = True,
    write_phase_sequence: bool = True,
) -> None:
    """
    Apply CSV edits to NetCDF files for a given date.

    Discovers which NetCDF files need updating by scanning the CSV directory for
    files whose names end in ``_sweeps.csv`` or ``_antenna_transition.csv``.  The
    corresponding NetCDF filename is derived by stripping that suffix and appending
    ``.nc``.  This means any scan type (vpt, man, ppi, …) can be edited without
    requiring a fixed filename pattern.

    Args:
        datestr: Date string in YYYYMMDD format
        data_dir: Base directory containing date subdirectories with NetCDF files
        csv_dir: Directory containing CSV files in YYYYMMDD subdirectories
                 (default: same as data_dir)
        backup: Create backup files before modifying (default: True)
        write_phase_sequence: Write phase_sequence global attribute from the phase column
            of the sweep CSV (default: True). Set to False for campaigns that do not use
            phase_sequence metadata (e.g. COBALT).

    Example:
        >>> apply_csv_edits_for_date('20171216', '/path/to/picasso/output')
        Processing date: 20171216
        Found 1 NetCDF file(s) to process (from 1 CSV file(s))

        Processing: ncas-mobile-ka-band-radar-1_cao_20171216-111356_vpt_l1_v1.0.0.nc
        Applying antenna transition edits...
        ...

        Completed: 1 file(s) processed, 0 skipped

    Note:
        When the CSV has a different number of sweeps than the NetCDF file,
        the entire file is rebuilt. When the sweep count is unchanged, only
        the specific fields are updated (faster operation).
    """
    import glob
    import os

    # If csv_dir not specified, use data_dir
    if csv_dir is None:
        csv_dir = data_dir

    # Construct the date directory paths
    data_date_path = os.path.join(data_dir, datestr)
    csv_date_path = os.path.join(csv_dir, datestr)

    if not os.path.exists(data_date_path):
        print(f"ERROR: Data directory not found: {data_date_path}")
        return

    if not os.path.exists(csv_date_path):
        print(f"ERROR: CSV directory not found: {csv_date_path}")
        return

    # Find all relevant CSV files and build a dict: nc_stem -> {sweep_csv, antenna_csv}
    csv_files = (
        sorted(glob.glob(os.path.join(csv_date_path, '*_antenna_transition.csv'))) +
        sorted(glob.glob(os.path.join(csv_date_path, '*_sweeps.csv')))
    )

    nc_map = {}
    for csv_file in csv_files:
        stem = os.path.basename(csv_file)
        if stem.endswith('_antenna_transition.csv'):
            nc_stem = stem[:-len('_antenna_transition.csv')]
            key = 'antenna_csv'
        elif stem.endswith('_sweeps.csv'):
            nc_stem = stem[:-len('_sweeps.csv')]
            key = 'sweep_csv'
        else:
            continue
        if nc_stem not in nc_map:
            nc_map[nc_stem] = {}
        nc_map[nc_stem][key] = csv_file

    if not nc_map:
        print(f"No CSV files found in: {csv_date_path}")
        return

    print(f"Processing date: {datestr}")
    print(f"Data directory: {data_date_path}")
    print(f"CSV directory:  {csv_date_path}")
    print(f"Found {len(nc_map)} NetCDF file(s) to process (from {len(csv_files)} CSV file(s))\n")

    processed = 0
    skipped = 0

    for nc_stem in sorted(nc_map.keys()):
        nc_file = os.path.join(data_date_path, nc_stem + '.nc')
        csv_info = nc_map[nc_stem]
        sweep_csv = csv_info.get('sweep_csv')
        antenna_csv = csv_info.get('antenna_csv')

        print(f"Processing: {nc_stem}.nc")

        if not os.path.exists(nc_file):
            print(f"  NetCDF file not found: {nc_file} - skipping")
            skipped += 1
            continue

        # Determine if we need to rebuild (sweep count changed) or just update
        needs_rebuild = False
        if sweep_csv is not None:
            import csv
            import netCDF4 as nc4

            # Count sweeps in CSV
            with open(sweep_csv, 'r') as f:
                reader = csv.DictReader(f)
                nsweeps_csv = sum(1 for _ in reader)

            # Count sweeps in NetCDF
            with nc4.Dataset(nc_file, 'r') as ds:
                nsweeps_nc = len(ds.dimensions['sweep'])

            if nsweeps_csv != nsweeps_nc:
                needs_rebuild = True
                print(f"  Sweep count changed ({nsweeps_nc} -> {nsweeps_csv}), rebuilding file...")

        # Apply the edits using appropriate method
        try:
            if needs_rebuild:
                # Rebuild entire file with new sweep structure
                rebuild_file_from_csv(
                    cfradial_file=nc_file,
                    sweep_csv=sweep_csv,
                    backup=backup
                )
                # Apply antenna transition edits after rebuild if CSV exists
                if antenna_csv is not None:
                    print("  Applying antenna transition edits...")
                    import_sweep_metadata_from_csv(
                        cfradial_file=nc_file,
                        sweep_csv=None,  # Don't update sweeps again
                        antenna_transition_csv=antenna_csv,
                        backup=False,  # Already backed up
                        write_phase_sequence=write_phase_sequence
                    )
            else:
                # Just update existing file (same sweep count)
                import_sweep_metadata_from_csv(
                    cfradial_file=nc_file,
                    sweep_csv=sweep_csv,
                    antenna_transition_csv=antenna_csv,
                    backup=backup,
                    write_phase_sequence=write_phase_sequence
                )
            processed += 1
        except Exception as e:
            print(f"  ERROR: Failed to process {nc_stem}.nc: {e}")
            import traceback
            traceback.print_exc()
            skipped += 1
        
        print()  # Blank line between files
    
    print(f"Completed: {processed} file(s) processed, {skipped} skipped")


def rebuild_file_from_csv(
    cfradial_file: str,
    sweep_csv: str = None,
    backup: bool = True
) -> None:
    """
    Rebuild a CF-Radial file with new sweep structure from edited CSV.
    
    This function recreates the NetCDF file with a new sweep structure based on
    the CSV file. Use this when you need to split or merge sweeps (change the
    number of sweeps), which cannot be done with import_sweep_metadata_from_csv.
    
    Args:
        cfradial_file: Path to CF-Radial NetCDF file to rebuild
        sweep_csv: Input path for sweep metadata CSV (default: <file>_sweeps.csv)
        backup: Create backup file before rebuilding (default: True)
        
    Example:
        >>> rebuild_file_from_csv('myfile.nc')
        Created backup: myfile.nc.bak
        Original file: 5 sweeps
        CSV file: 7 sweeps (2 sweeps added)
        Reading radar data...
        Applying new sweep structure...
        Writing updated file...
        Successfully rebuilt myfile.nc with 7 sweeps
        
    Note:
        This recreates the entire file, which takes longer than import_sweep_metadata_from_csv
        but allows changing the number of sweeps.
    """
    import csv
    import netCDF4 as nc4
    from pathlib import Path
    import shutil
    import pyart
    from campaign_processing import apply_phase_sweeps_to_radar
    
    # Generate default input path if not provided
    if sweep_csv is None:
        sweep_csv = str(Path(cfradial_file).with_suffix('')) + '_sweeps.csv'
    
    if not Path(sweep_csv).exists():
        print(f"ERROR: Sweep CSV not found: {sweep_csv}")
        return
    
    # Read sweep metadata from CSV
    sweep_data = []
    phases_list = []
    sweep_modes_list = []
    with open(sweep_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sweep_data.append({
                'sweep_number': int(row['sweep_number']),
                'start_idx': int(row['start_idx']),
                'end_idx': int(row['end_idx']),
                'fixed_angle': float(row['fixed_angle'])
            })
            # Read phase if it exists (campaign-specific)
            if 'phase' in row:
                phases_list.append(row['phase'].strip())
            # Read sweep_mode if it exists (CF-Radial standard)
            if 'sweep_mode' in row:
                sweep_modes_list.append(row['sweep_mode'].strip())
            # If neither exists, use default
            if 'phase' not in row and 'sweep_mode' not in row:
                phases_list.append('tracking')  # Default fallback
    
    # Check original file sweep count
    with nc4.Dataset(cfradial_file, 'r') as ds:
        nsweeps_original = len(ds.dimensions['sweep'])
    
    nsweeps_csv = len(sweep_data)
    print(f"Original file: {nsweeps_original} sweeps")
    print(f"CSV file: {nsweeps_csv} sweeps", end='')
    if nsweeps_csv != nsweeps_original:
        diff = nsweeps_csv - nsweeps_original
        print(f" ({diff:+d} sweeps {'added' if diff > 0 else 'removed'})")
    else:
        print()
    
    # Create backup if requested
    if backup:
        backup_file = cfradial_file + '.bak'
        shutil.copy2(cfradial_file, backup_file)
        print(f"Created backup: {backup_file}")
    
    # Read radar data
    print("Reading radar data...")
    try:
        radar = pyart.io.read_cfradial(cfradial_file)
    except Exception as e:
        print(f"ERROR: Failed to read radar file: {e}")
        return
    
    # Build phases list from CSV (format expected by apply_phase_sweeps_to_radar)
    print("Applying new sweep structure...")
    phases = []
    for i, sweep in enumerate(sweep_data):
        # Use phase directly from CSV
        phase_type = phases_list[i]
        
        # Calculate mean azimuth and elevation from the ray indices
        start_idx = sweep['start_idx']
        end_idx = sweep['end_idx']
        mean_az = float(np.mean(radar.azimuth['data'][start_idx:end_idx+1]))
        mean_el = float(np.mean(radar.elevation['data'][start_idx:end_idx+1]))
        
        phases.append({
            'phase': phase_type,
            'start_idx': start_idx,
            'end_idx': end_idx,
            'n_rays': end_idx - start_idx + 1,
            'mean_az': mean_az,
            'mean_el': mean_el,
            'transition_ray_indices': []  # Empty for manual CSV structure
        })
    
    # Apply the new sweep structure
    try:
        radar = apply_phase_sweeps_to_radar(radar, phases)
    except Exception as e:
        print(f"ERROR: Failed to apply sweep structure: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Write the file
    print("Writing updated file...")
    try:
        pyart.io.write_cfradial(cfradial_file, radar, format='NETCDF4')
    except Exception as e:
        print(f"ERROR: Failed to write file: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Update sweep_mode variable if values were provided in CSV
    if sweep_modes_list and len(sweep_modes_list) == nsweeps_csv:
        with nc4.Dataset(cfradial_file, 'a') as ds:
            if 'sweep_mode' in ds.variables:
                sweep_mode_var = ds.variables['sweep_mode']
                # Handle character array (old-style NetCDF3 strings)
                if sweep_mode_var.dtype.char == 'S':
                    # Convert list of strings to character array
                    sweep_mode_array = np.array(sweep_modes_list, dtype='S32')
                    char_array = nc4.stringtochar(sweep_mode_array)
                    sweep_mode_var[:] = char_array
                else:
                    # Direct string assignment for NetCDF4 string types
                    for i, mode in enumerate(sweep_modes_list):
                        sweep_mode_var[i] = mode
                print(f"Updated sweep_mode variable from CSV for {nsweeps_csv} sweeps")
    
    # Update time coverage attributes
    cfradial_add_time_coverage(cfradial_file)

    # Update history
    history_msg = f"Rebuilt file with new sweep structure from CSV ({nsweeps_csv} sweeps)"
    update_history_attribute(cfradial_file, history_msg)

    # Update last_revised_date
    import datetime
    with nc4.Dataset(cfradial_file, 'a') as ds:
        current_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
        ds.setncattr('last_revised_date', current_time)

    print(f"Successfully rebuilt {cfradial_file} with {nsweeps_csv} sweeps")


def is_sweep_all_antenna_transition(radar: Radar, sweep_idx: int) -> bool:
    """
    Check if a sweep consists entirely of antenna transition rays.
    
    Args:
        radar: PyART Radar object
        sweep_idx: Sweep index to check
        
    Returns:
        True if all rays in the sweep have antenna_transition=1, False otherwise
        
    Example:
        >>> radar = pyart.io.read_cfradial('myfile.nc')
        >>> is_transition = is_sweep_all_antenna_transition(radar, 3)
        >>> if is_transition:
        ...     print("Sweep 3 is all antenna transitions")
    """
    # Check if antenna_transition coordinate exists (stored as radar attribute, not field)
    if not hasattr(radar, 'antenna_transition') or radar.antenna_transition is None:
        return False
    
    if 'data' not in radar.antenna_transition:
        return False
    
    # Get ray indices for this sweep
    sweep_slice = radar.get_slice(sweep_idx)
    antenna_transition_data = radar.antenna_transition['data'][sweep_slice]
    
    # Check if all values are 1 (transition)
    return np.all(antenna_transition_data == 1)


def get_valid_sweep_indices(
    radar: Radar, 
    skip_all_transition: bool = False
) -> List[int]:
    """
    Get list of valid sweep indices, optionally filtering out all-transition sweeps.
    
    Args:
        radar: PyART Radar object
        skip_all_transition: If True, exclude sweeps where 100% of rays are antenna_transition=1
        
    Returns:
        List of valid sweep indices
        
    Example:
        >>> radar = pyart.io.read_cfradial('myfile.nc')
        >>> # Get all sweeps
        >>> all_sweeps = get_valid_sweep_indices(radar)
        >>> # Get only sweeps with actual data (not pure transitions)
        >>> data_sweeps = get_valid_sweep_indices(radar, skip_all_transition=True)
        >>> print(f"Total sweeps: {len(all_sweeps)}, Data sweeps: {len(data_sweeps)}")
    """
    nsweeps = radar.nsweeps
    
    if not skip_all_transition:
        return list(range(nsweeps))
    
    valid_indices = []
    for sweep_idx in range(nsweeps):
        if not is_sweep_all_antenna_transition(radar, sweep_idx):
            valid_indices.append(sweep_idx)
    
    return valid_indices
