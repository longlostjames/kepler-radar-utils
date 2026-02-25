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
        azimuth, elevation, scan_rate, antenna_transition, target_scan_rate = _extract_angle_info(
            ncobj, scan_name, time, revised_northangle, filemetadata
        )
        
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
    
    # Parse latitude
    lat_str = ncobj.getncattr('Latitude')
    z = StringIO(lat_str)
    lat_data = np.genfromtxt(z, dtype=None, names=['lat', 'direction'])
    
    if lat_data['direction'] == b'S' and lat_data['lat'] > 0:
        latitude['data'] = np.array([-lat_data['lat']], dtype='f4')
    else:
        latitude['data'] = np.array([lat_data['lat']], dtype='f4')
    
    # Parse longitude  
    lon_str = ncobj.getncattr('Longitude')
    z = StringIO(lon_str)
    lon_data = np.genfromtxt(z, dtype=None, names=['lon', 'direction'])
    
    if lon_data['direction'] == b'W' and lon_data['lon'] > 0:
        longitude['data'] = np.array([-lon_data['lon']], dtype='f4')
    else:
        longitude['data'] = np.array([lon_data['lon']], dtype='f4')
    
    # Parse altitude
    alt_str = ncobj.getncattr('Altitude')
    z = StringIO(alt_str)
    alt_data = np.genfromtxt(z, dtype=None, names=['alt', 'units'])
    altitude['data'] = np.array([alt_data['alt']], dtype='f4')
    
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
    elif scan_name in ['ppi', 'vertical_pointing']:
        # For PPI/VPT, fixed angle is the elevation
        fixed_angle_value = np.round(ncvars['elv'][i4fixed_angle], 2)
        fixed_angle["data"] = np.array([fixed_angle_value], dtype='f')
    
    print(f"Detected scan type: {scan_name}, fixed angle: {fixed_angle['data'][0]}")
    
    return scan_name, sweep_mode, fixed_angle

def _extract_angle_info(
    ncobj: nc4.Dataset, 
    scan_name: str, 
    time: Dict, 
    revised_northangle: float,
    filemetadata: FileMetadata
) -> Tuple[Dict, Dict, Optional[Dict], Optional[Dict], Optional[Dict]]:
    """
    Extract azimuth, elevation and scan rate information.
    
    Args:
        ncobj: Open netCDF4 Dataset
        scan_name: Type of scan (ppi, rhi, etc.)
        time: Time information dictionary
        revised_northangle: North angle correction in degrees
        filemetadata: PyART FileMetadata object
        
    Returns:
        Tuple of (azimuth, elevation, scan_rate, antenna_transition, target_scan_rate)
    """
    ncvars = ncobj.variables
    
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    
    # Basic angle assignments
    azimuth['data'] = (ncvars['azi'][:] + revised_northangle) % 360
    elevation['data'] = ncvars['elv'][:]
    
    # Calculate ray durations for scan rate correction
    ray_duration = np.diff(time['data'])
    
    if scan_name in ['ppi', 'rhi']:
        # For scanning modes, calculate scan rates and antenna transitions
        scan_rate = filemetadata("scan_rate")
        antenna_transition = filemetadata("antenna_transition") 
        target_scan_rate = filemetadata("target_scan_rate")
        
        # Determine primary scan rate variable based on scan type
        if scan_name == 'ppi':
            scan_rate['data'] = ncvars['aziv'][:]
        else:  # rhi
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
        
        # Apply timing corrections to angles
        _apply_timing_corrections(azimuth, elevation, ncvars, ray_duration, revised_northangle)
        
    else:
        # For non-scanning modes (VPT), no scan rate information
        scan_rate = None
        antenna_transition = None  
        target_scan_rate = None
    
    return azimuth, elevation, scan_rate, antenna_transition, target_scan_rate

def _apply_timing_corrections(
    azimuth: Dict, 
    elevation: Dict, 
    ncvars: Dict, 
    ray_duration: np.ndarray,
    revised_northangle: float
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
    """
    # Calculate target ray duration and identify anomalies
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
    azimuth['proposed_standard_name'] = "sensor_to_target_azimuth_angle"
    azimuth['long_name'] = "sensor to target azimuth angle"
    
    # Special handling for long-duration rays
    for idx in long_duration:
        if idx < len(azimuth['data']):
            azimuth['data'][idx] = (ncvars['azi'][idx] + revised_northangle) % 360
            azimuth['data'][idx] -= 0.5 * target_ray_duration * ncvars['aziv'][idx]
    
    # Apply corrections for elevation
    elevation['data'] -= 0.5 * ray_duration_extended * ncvars['elvv'][:]
    elevation['units'] = "degrees"
    elevation['proposed_standard_name'] = "sensor_to_target_elevation_angle" 
    elevation['long_name'] = "sensor to target elevation angle"
    
    # Special handling for long-duration rays
    for idx in long_duration:
        if idx < len(elevation['data']):
            elevation['data'][idx] = ncvars['elv'][idx]
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
                # Convert linear to log scale
                field_dict['data'] = 10.0 * np.log10(ncvars[mmclx_name][:])
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
        prt = 1.0 / ncvars["prf"][:]
        prt_dict["data"] = np.ones((nrays,), dtype="float32") * prt
        instrument_parameters["prt"] = prt_dict
    
    # Parse hardware parameters from 'hrd' attribute
    hrd_variables = {}
    if 'hrd' in ncobj.ncattrs():
        hrd_attribute_value = ncobj.getncattr('hrd')
        hrd_variables = _parse_hrd_attribute(hrd_attribute_value)
    
    # Pulse width
    if "PULSE_WIDTH" in hrd_variables:
        pw_dict = filemetadata("pulse_width")
        pulse_width = hrd_variables["PULSE_WIDTH"]
        pw_dict["data"] = np.ones((nrays,), dtype="float32") * pulse_width
        instrument_parameters["pulse_width"] = pw_dict
    
    # Transmit power
    if "tpow" in ncvars:
        tp_dict = filemetadata("radar_measured_transmit_power_h")
        txpower = 10.0 * np.log10(ncvars["tpow"][:]) + 30  # Convert to dBm
        tp_dict["data"] = txpower
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
    for n in project["ncas_instruments"]:
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
    
    # Add NCAS metadata
    cfradial_add_ncas_metadata(outfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version)
    
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
            time_var.long_name = "time in seconds since time_reference"

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
    
    # Get the date string from start time to look in the right subdirectory
    date_str = start_datetime.strftime('%Y%m%d')
    search_path = Path(inpath) / date_str
    
    print(f"Looking in date directory: {search_path}")
    
    if not search_path.exists():
        print(f"Date directory does not exist: {search_path}")
        return []
    
    # Define search patterns based on sweep type
    if gzip_flag:
        patterns = {
            'rhi': "*rhi*.mmclx.gz",
            'ppi': "*ppi*.mmclx.gz", 
            'vert': "*vert*.mmclx.gz",
            'vad': "*ppi*.mmclx.gz"  # VAD files are typically PPI scans
        }
    else:
        patterns = {
            'rhi': "*rhi*.mmclx",
            'ppi': "*ppi*.mmclx",
            'vert': "*vert*.mmclx", 
            'vad': "*ppi*.mmclx"  # VAD files are typically PPI scans
        }
    
    pattern = patterns.get(sweep_type.lower(), f"*{sweep_type}*.mmclx{'gz' if gzip_flag else ''}")
    
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
    search_path = Path(inpath) / date_str
    
    print(f"Looking in date directory: {search_path}")
    
    if not search_path.exists():
        print(f"Date directory does not exist: {search_path}")
        return []
    
    # Define a list to store the found files
    matching_files = []
    
    az_search_offset = -8.0  # Before July
    
    # Look for RHI files
    if gzip_flag:
        rhi_pattern = "*rhi*.mmclx.gz"
    else:
        rhi_pattern = "*rhi*.mmclx"
    
    candidate_files = list(search_path.glob(rhi_pattern))
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
    
    # Get the date string from start time to look in the right subdirectory
    date_str = start_datetime.strftime('%Y%m%d')
    search_path = Path(inpath) / date_str
    
    print(f"Looking in date directory: {search_path}")
    
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
    data_version: str
) -> None:
    """
    Add NCAS metadata to a CF-Radial file.
    
    Args:
        cfradial_file: Path to CF-Radial file to modify
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file  
        tracking_tag: AMOF tracking tag
        data_version: Data version string
        
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
    
    try:
        # Try to use the official NCAS metadata library
        from ncas_amof_netcdf_template import util
        
        util.add_ncas_metadata_to_netcdf_file(
            cfradial_file, 
            yaml_project_file, 
            yaml_instrument_file, 
            tracking_tag,
            data_version
        )
        
        print(f"Successfully added NCAS metadata using official library")
        return
        
    except ImportError as e:
        print(f"NCAS metadata library not available: {e}")
        print("Falling back to manual metadata addition")
    except Exception as e:
        print(f"Error using NCAS metadata library: {e}")
        print("Falling back to manual metadata addition")
    
    # Fallback: Add NCAS metadata manually (but still require YAML files)
    try:
        _add_ncas_metadata_manually(cfradial_file, yaml_project_file, yaml_instrument_file, tracking_tag, data_version)
    except Exception as e:
        raise RuntimeError(f"Error adding NCAS metadata to {cfradial_file}: {e}")

def _add_ncas_metadata_manually(
    cfradial_file: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    tracking_tag: str,
    data_version: str
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
    
    # Also look for instrument info in the project YAML file under ncas_instruments
    project_instrument_info = None
    
    # Check if project_info has ncas_instruments
    if project_info and isinstance(project_info, dict) and 'ncas_instruments' in project_info:
        instruments_list = project_info['ncas_instruments']
        if isinstance(instruments_list, list):
            for instrument in instruments_list:
                if isinstance(instrument, dict):
                    # Handle dict with instrument name as key
                    for inst_name, inst_data in instrument.items():
                        # Check for exact match or partial match (case insensitive)
                        if inst_name == 'ncas-mobile-ka-band-radar-1' or 'ncas-mobile-ka-band-radar-1' in inst_name.lower():
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
            
            # Update core attributes
            ds.setncattr('title', 'Moment data from NCAS Mobile Ka-band Radar (Kepler)')
            ds.setncattr('institution', 'National Centre for Atmospheric Science (NCAS) and Science and Technology Facilities Council (STFC) as part of UK Research and Innovation (UKRI)')
            ds.setncattr('source', 'NCAS Mobile Ka-band Radar unit 1')
            ds.setncattr('references', 'www.metek.de')
            ds.setncattr('product_version', f'v{data_version}')
            ds.setncattr('processing_level', '1')
            ds.setncattr('licence', 'This dataset is released for use under the Open Government Licence, OGL-UK-3.0 (see https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).')
            
            # Add creator information
            ds.setncattr('creator_name', 'Chris Walden')
            ds.setncattr('creator_email', 'chris.walden@ncas.ac.uk')
            ds.setncattr('creator_url', 'https://orcid.org/0000-0002-5718-466X')
            
            # Add instrument metadata from YAML if available
            if instrument_info:
                print("Adding instrument metadata from YAML")
                if 'manufacturer' in instrument_info:
                    ds.setncattr('instrument_manufacturer', instrument_info['manufacturer'])
                if 'model' in instrument_info:
                    ds.setncattr('instrument_model', instrument_info['model'])
                if 'serial_number' in instrument_info:
                    ds.setncattr('instrument_serial_number', instrument_info['serial_number'])
                if 'pid' in instrument_info:
                    ds.setncattr('instrument_pid', instrument_info['pid'])
            else:
                print("Using default instrument metadata")
                # Default instrument metadata
                ds.setncattr('instrument_manufacturer', 'Meteorologische Messtechnik (Metek) GmbH')
                ds.setncattr('instrument_model', 'MIRA-35S')
                ds.setncattr('instrument_serial_number', 'BX3C')
                ds.setncattr('instrument_pid', 'https://hdl.handle.net/21.12132/3.b7a8298b7a54405f')
            
            # Add project metadata from YAML if available
            if project_info:
                print("Adding project metadata from YAML")
                if 'title' in project_info:
                    ds.setncattr('project', project_info['title'])
                
                if 'principal_investigator' in project_info:
                    pi = project_info['principal_investigator']
                    if isinstance(pi, dict):
                        if 'name' in pi:
                            ds.setncattr('project_principal_investigator', pi['name'])
                        if 'email' in pi:
                            ds.setncattr('project_principal_investigator_email', pi['email'])
                        if 'orcid' in pi:
                            ds.setncattr('project_principal_investigator_url', pi['orcid'])
                
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
                # Default project metadata for COBALT
                ds.setncattr('project', 'Contrail Observations and Lifecycle Tracking (COBALT)')
                ds.setncattr('project_principal_investigator', 'Edward Gryspeerdt')
                ds.setncattr('project_principal_investigator_email', 'e.gryspeerdt@imperial.ac.uk')
                ds.setncattr('project_principal_investigator_url', 'https://orcid.org/0000-0002-3815-4756')
                ds.setncattr('acknowledgement', 'This dataset was developed as part of the activity "Contrail Observations and Lifecycle Tracking (COBALT)", funded by\nNatural Environment Research Council (NERC) Grant NE/Z503794/1.\nIt uses instrumentation provided by the Atmospheric Measurement and Observation Facility (AMOF), part of NERC National Capability.\nUsers should acknowledge the National Centre for Atmospheric Science (NCAS) as the data provider.')
            
            # Platform information
            ds.setncattr('platform', 'CAO')
            ds.setncattr('platform_type', 'stationary_platform')
            ds.setncattr('location_keywords', 'Chilbolton, Hampshire, England')
            ds.setncattr('platform_is_mobile', 'false')
            ds.setncattr('deployment_mode', 'land')
            
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
            
            # Update variable attributes for NCAS compliance
            _update_variable_attributes_for_ncas(ds)
        
    except Exception as e:
        raise RuntimeError(f"Error modifying NetCDF file {cfradial_file}: {e}")
    
    print(f"Successfully added manual NCAS metadata to {cfradial_file}")

def _update_variable_attributes_for_ncas(ds):
    """Update variable attributes to be NCAS compliant."""
    
    # Update time variable
    if 'time' in ds.variables:
        time_var = ds.variables['time']
        time_var.setncattr('comment', '')
        time_var.setncattr('long_name', 'time in seconds since time_reference')
    
    # Update range variable  
    if 'range' in ds.variables:
        range_var = ds.variables['range']
        range_var.setncattr('comment', 'Range to centre of each bin')
        if hasattr(range_var, 'meters_to_center_of_first_gate'):
            # Keep existing value
            pass
        else:
            range_var.setncattr('meters_to_center_of_first_gate', range_var[0])
    
    # Update azimuth variable
    if 'azimuth' in ds.variables:
        azim_var = ds.variables['azimuth']
        azim_var.setncattr('comment', 'Azimuth of antenna relative to true north')
    
    # Update elevation variable
    if 'elevation' in ds.variables:
        elev_var = ds.variables['elevation']
        elev_var.setncattr('comment', 'Elevation of antenna relative to the horizontal plane')
    
    # Update altitude variable
    if 'altitude' in ds.variables:
        alt_var = ds.variables['altitude']
        alt_var.setncattr('comment', 'Altitude of the centre of rotation of the antenna above the geoid using the WGS84 ellipsoid and EGM2008 geoid model')
        alt_var.setncattr('units', 'metres')
    
    # Update sweep_mode variable
    if 'sweep_mode' in ds.variables:
        sweep_var = ds.variables['sweep_mode']
        sweep_var.setncattr('comment', 'Options are: "sector", "coplane", "rhi", "vertical_pointing", "idle", "azimuth_surveillance", "elevation_surveillance", "sunscan", "pointing", "manual_ppi", "manual_rhi"')
        sweep_var.setncattr('units', '')
    
    # Update antenna_transition variable
    if 'antenna_transition' in ds.variables:
        ant_var = ds.variables['antenna_transition']
        ant_var.setncattr('comment', '1 if antenna is in transition, 0 otherwise')
        ant_var.setncattr('units', '')
    
    # Add frequency variable if missing
    if 'frequency' not in ds.variables and 'frequency' not in ds.dimensions:
        # Add frequency dimension and variable
        ds.createDimension('frequency', 1)
        freq_var = ds.createVariable('frequency', 'f4', ('frequency',))
        freq_var.setncattr('standard_name', 'radiation_frequency')
        freq_var.setncattr('long_name', 'frequency of transmitted radiation')
        freq_var.setncattr('units', 's-1')
        freq_var.setncattr('meta_group', 'instrument_parameters')
        freq_var[:] = [35.5e9]  # 35.5 GHz for Ka-band
    
    # Add azimuth_correction variable if missing
    if 'azimuth_correction' not in ds.variables:
        azim_corr_var = ds.createVariable('azimuth_correction', 'f4')
        azim_corr_var.setncattr('long_name', 'azimuth correction applied')
        azim_corr_var.setncattr('units', 'degrees')
        azim_corr_var.setncattr('meta_group', 'geometry_correction')
        azim_corr_var.setncattr('comment', 'Azimuth correction applied. North angle relative to instrument home azimuth.')
        azim_corr_var[:] = 55.9  # Default value for COBALT

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
    yaml_instrument_file: str = None
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
                revised_northangle=revised_northangle
            )
            radars.append(radar)
        except Exception as e:
            print(f"Error reading {mmclx_file}: {e}")
            continue
    
    if not radars:
        print("No valid radar files could be read")
        return None
    
    # Determine location from campaign
    location_map = {
        'woest': 'lyneham',
        'ccrest': 'cao',
        'ccrest-m': 'cao',
        'coalesc3': 'wao',
        'cobalt': 'cao',
        'kasbex': 'cao'
    }
    location = location_map.get(campaign.lower(), 'unknown')
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
                
                # Add time coverage attributes
                cfradial_add_time_coverage(outfile)
                
                # Add NCAS metadata
                cfradial_add_ncas_metadata(outfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version)
                
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
        
        # Add time coverage attributes
        cfradial_add_time_coverage(outfile)
        
        # Add NCAS metadata
        cfradial_add_ncas_metadata(outfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version)
        
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
    
    # Create new angle arrays
    azimuth_data = []
    elevation_data = []
    for radar in radar_list:
        azimuth_data.extend(radar.azimuth['data'])
        elevation_data.extend(radar.elevation['data'])
    
    combined_azimuth = base_radar.azimuth.copy()
    combined_azimuth['data'] = np.array(azimuth_data)
    
    combined_elevation = base_radar.elevation.copy()
    combined_elevation['data'] = np.array(elevation_data)
    
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
    
    # Create combined radar object
    combined_radar = Radar(
        combined_time, base_radar.range, combined_fields,
        base_radar.metadata, scan_name.lower(),
        base_radar.latitude, base_radar.longitude, base_radar.altitude,
        sweep_number, sweep_mode, fixed_angle,
        sweep_start_ray_index, sweep_end_ray_index,
        combined_azimuth, combined_elevation,
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
    
    all_azimuths = np.array(all_azimuths)
    all_elevations = np.array(all_elevations)
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
    
    # Create combined radar object
    combined_radar = Radar(
        combined_time, base_radar.range, combined_fields,
        base_radar.metadata, 'ppi',
        base_radar.latitude, base_radar.longitude, base_radar.altitude,
        sweep_number, sweep_mode, fixed_angle,
        sweep_start_ray_index, sweep_end_ray_index,
        combined_azimuth, combined_elevation,
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
    
    # Get the date string from start time to look in the right subdirectory
    date_str = start_datetime.strftime('%Y%m%d')
    search_path = Path(inpath) / date_str
    
    print(f"Looking in date directory: {search_path}")
    
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



