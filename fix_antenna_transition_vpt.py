"""
fix_antenna_transition_vpt.py
------------------------------
Set antenna_transition=1 for any ray in a VPT L1 CF-Radial file whose
elevation angle is below a threshold (default 80°).

This catches radar rays that were recorded during antenna motion between
interleaved VPT stares and RHI scans — e.g. when the MIRA SIRTA scheduler
alternates scan modes — but were not flagged by the upstream processor.

Only rays that are both low-elevation AND currently unflagged (=0) are
modified.  The file is edited in-place; a dry-run mode is available.

Usage
-----
    python fix_antenna_transition_vpt.py [options] <file> [<file> ...]

Options
-------
    -t, --threshold  Elevation threshold in degrees (default: 80.0)
    -n, --dry-run    Print what would change without modifying any file
    -h, --help       Show this message

Examples
--------
    # Fix a single file
    python fix_antenna_transition_vpt.py ncas-..._vpt_l1_v1.0.3.nc

    # Fix all VPT files for a range of dates
    python fix_antenna_transition_vpt.py /path/to/L1/202412*/ncas-*_vpt_l1*.nc

    # Check but don't modify
    python fix_antenna_transition_vpt.py --dry-run /path/to/L1/202412*/ncas-*_vpt_l1*.nc
"""

import datetime
import getopt
import glob
import os
import socket
import sys

import netCDF4 as nc4
import numpy as np


ELEVATION_THRESHOLD_DEFAULT = 80.0


def _history_entry(msg):
    """Return a history line: 'DDD Mon DD HH:MM:SS YYYY - user:<u> machine:<h> <msg>'"""
    now = datetime.datetime.utcnow()
    user = os.environ.get('USER', 'unknown')
    host = socket.gethostname()
    timestamp = now.strftime('%a %b %d %H:%M:%S %Y')
    return f'{timestamp} - user:{user} machine:{host} {msg}'


def fix_file(filepath, threshold=ELEVATION_THRESHOLD_DEFAULT, dry_run=False):
    """
    Inspect *filepath* and set antenna_transition=1 for rays with
    elevation < *threshold* degrees that are currently flagged as 0.

    Parameters
    ----------
    filepath : str
        Path to a VPT L1 CF-Radial NetCDF file.
    threshold : float
        Elevation angle (degrees) below which a ray is considered a
        transition ray.
    dry_run : bool
        If True, report what would change but do not modify the file.

    Returns
    -------
    int
        Number of rays updated (0 if the file was already correct).
    """
    if not os.path.isfile(filepath):
        print(f"ERROR: file not found: {filepath}", file=sys.stderr)
        return 0

    with nc4.Dataset(filepath, 'r') as ds:
        if 'elevation' not in ds.variables:
            print(f"SKIP: no elevation variable in {filepath}", file=sys.stderr)
            return 0
        if 'antenna_transition' not in ds.variables:
            print(f"SKIP: no antenna_transition variable in {filepath}", file=sys.stderr)
            return 0

        elev  = ds.variables['elevation'][:]
        anttr = ds.variables['antenna_transition'][:]

    # Rays that need fixing: low-elevation AND currently NOT flagged
    to_fix = np.where((elev < threshold) & (anttr == 0))[0]
    n_fix = len(to_fix)
    already_flagged = int((elev < threshold).sum()) - n_fix

    if n_fix == 0:
        status = "already correct" if already_flagged == 0 else f"all {already_flagged} transition ray(s) already flagged"
        print(f"OK  {os.path.basename(filepath)} — {status}")
        return 0

    print(f"{'(dry-run) ' if dry_run else ''}FIX {os.path.basename(filepath)} "
          f"— setting antenna_transition=1 for {n_fix} ray(s) with elevation < {threshold}°")

    if dry_run:
        return n_fix

    # Apply in-place edit
    with nc4.Dataset(filepath, 'r+') as ds:
        ds.variables['antenna_transition'][to_fix] = np.int8(1)

        # Update history
        msg = (f"fix_antenna_transition_vpt.py: set antenna_transition=1 for "
               f"{n_fix} ray(s) with elevation < {threshold}°")
        entry = _history_entry(msg)
        existing_history = getattr(ds, 'history', '') or ''
        ds.setncattr('history', f"{existing_history}\n{entry}" if existing_history else entry)

        # Update last_modified
        now_iso = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
        ds.setncattr('last_modified', now_iso)

    return n_fix


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'ht:n', ['help', 'threshold=', 'dry-run'])
    except getopt.GetoptError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    threshold = ELEVATION_THRESHOLD_DEFAULT
    dry_run = False

    for opt, val in opts:
        if opt in ('-h', '--help'):
            print(__doc__)
            sys.exit(0)
        elif opt in ('-t', '--threshold'):
            threshold = float(val)
        elif opt in ('-n', '--dry-run'):
            dry_run = True

    if not args:
        print("Error: no input files specified.", file=sys.stderr)
        print("Usage: python fix_antenna_transition_vpt.py [options] <file> [<file> ...]",
              file=sys.stderr)
        sys.exit(1)

    # Expand any globs (helpful when called from a shell that doesn't expand them)
    files = []
    for pattern in args:
        expanded = sorted(glob.glob(pattern))
        files.extend(expanded if expanded else [pattern])

    total_fixed = 0
    for filepath in files:
        total_fixed += fix_file(filepath, threshold=threshold, dry_run=dry_run)

    print(f"\n{'Would update' if dry_run else 'Updated'} {total_fixed} ray(s) across {len(files)} file(s).")


if __name__ == '__main__':
    main()
