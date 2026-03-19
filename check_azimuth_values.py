#!/usr/bin/env python
"""
Quick diagnostic script to check azimuth values and corrections in NetCDF files.
Usage: python check_azimuth_values.py <netcdf_file>
"""

import sys
import netCDF4 as nc4
import numpy as np

if len(sys.argv) < 2:
    print("Usage: python check_azimuth_values.py <netcdf_file>")
    sys.exit(1)

ncfile = sys.argv[1]

print(f"\nChecking azimuth information in: {ncfile}\n")
print("=" * 70)

with nc4.Dataset(ncfile, 'r') as ds:
    # Check for azimuth_correction variable
    if 'azimuth_correction' in ds.variables:
        azimuth_corr = ds.variables['azimuth_correction'][:]
        print(f"✓ azimuth_correction variable found: {azimuth_corr}")
        if hasattr(ds.variables['azimuth_correction'], 'comment'):
            print(f"  Comment: {ds.variables['azimuth_correction'].comment}")
    else:
        print("✗ azimuth_correction variable NOT found")
    
    print()
    
    # Check azimuth values
    if 'azimuth' in ds.variables:
        azimuth_data = ds.variables['azimuth'][:]
        print(f"Azimuth data:")
        print(f"  Min: {np.min(azimuth_data):.2f}°")
        print(f"  Max: {np.max(azimuth_data):.2f}°")
        print(f"  Mean: {np.mean(azimuth_data):.2f}°")
        print(f"  First 10 values: {azimuth_data[:10]}")
        
        if hasattr(ds.variables['azimuth'], 'long_name'):
            print(f"  long_name: {ds.variables['azimuth'].long_name}")
    else:
        print("✗ azimuth variable NOT found")
    
    print()
    
    # Check history for north_angle information
    if hasattr(ds, 'history'):
        history = ds.history
        print("History attribute:")
        # Look for north_angle mentions
        for line in history.split('\n'):
            if 'north' in line.lower() or 'azimuth' in line.lower():
                print(f"  {line}")
    
    print()
    
    # Check if this is a PPI scan
    if 'sweep_mode' in ds.variables:
        sweep_mode = ds.variables['sweep_mode'][:]
        print(f"Sweep mode: {sweep_mode}")
    
    print()
    
    # Check scan_type
    if hasattr(ds, 'scan_type'):
        print(f"Scan type: {ds.scan_type}")

print("=" * 70)
