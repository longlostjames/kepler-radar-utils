"""
fix_empty_history_l1.py
-----------------------
One-off fix: write a non-empty 'history' global attribute to every CCREST-M L1
NetCDF file that currently has history = '' (empty string).

The NCAS checker requires that history is non-empty.

Usage:
    conda run -n cao_3_11 python3 fix_empty_history_l1.py [--dry-run]
"""

import datetime
import glob
import os
import socket
import sys

import netCDF4 as nc4

BASE = (
    '/gws/pw/j07/ncas_obs_vol2/cao/processing/'
    'ncas-mobile-ka-band-radar-1/ccrest-m/L1_v1.0.1'
)

DRY_RUN = '--dry-run' in sys.argv

timestamp = datetime.datetime.utcnow().strftime('%a %b %d %H:%M:%S %Y')
username = os.environ.get('USER', 'unknown')
hostname = os.environ.get('HOSTNAME', socket.gethostname())
HISTORY_ENTRY = (
    f"{timestamp} - user:{username} machine:{hostname} "
    f"Convert mmclx to CF-Radial with NCAS metadata"
)

files = sorted(glob.glob(os.path.join(BASE, '*', '*_l1_*.nc')))
print(f"Found {len(files)} L1 files")

fixed = 0
skipped = 0
errors = 0

for f in files:
    try:
        with nc4.Dataset(f, 'r') as ds:
            current = ds.getncattr('history') if hasattr(ds, 'history') else None

        if current and current.strip():
            skipped += 1
            continue

        if DRY_RUN:
            print(f"[DRY RUN] Would fix: {os.path.basename(f)}")
            fixed += 1
            continue

        with nc4.Dataset(f, 'a') as ds:
            ds.setncattr('history', HISTORY_ENTRY)

        fixed += 1
        print(f"Fixed: {os.path.basename(f)}")

    except Exception as e:
        print(f"ERROR: {os.path.basename(f)}: {e}")
        errors += 1

print()
print(f"Fixed  : {fixed}")
print(f"Skipped: {skipped} (already had history)")
print(f"Errors : {errors}")
