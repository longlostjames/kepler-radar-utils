#!/usr/bin/env python3
"""
Fix antenna_transition flags in a COBALT VPT L1 file where mmclx elevation
fill values (elv = -1000) were silently replaced by 90 degrees during L1
processing.

Strategy:
  1. Scan all *vert.mmclx.gz files for the given date.
  2. Collect Unix timestamps of rays where elv < -999 (the sentinel fill value).
  3. Convert those timestamps to the L1 file's time units.
  4. Flag every L1 ray whose time falls within ±TOLERANCE seconds of the
     fill-value period as antenna_transition = 1.

Usage:
    python fix_vpt_transitions_from_mmclx.py \\
        -d 20241216 \\
        --mmclx-dir /gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-mobile-ka-band-radar-1/data/campaign/cobalt/mom/20241216 \\
        --l1-file /gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/cobalt/L1_v1.0.3/20241216/ncas-mobile-ka-band-radar-1_cao_20241216-000007_vpt_l1_v1.0.3.nc \\
        [--dry-run] [--no-backup]
"""

import argparse
import datetime
import glob
import gzip
import os
import shutil
import tempfile

import netCDF4 as nc4
import numpy as np


# Rays with mmclx elv below this value are treated as fill / unknown elevation.
ELV_FILL_THRESHOLD = -999.0

# Maximum mismatch (seconds) between an mmclx ray timestamp and the L1 time
# array when deciding whether a contiguous fill-value block boundary applies.
TOLERANCE_S = 3.0


# ---------------------------------------------------------------------------
# Step 1: scan mmclx files
# ---------------------------------------------------------------------------

def scan_mmclx_fill_rays(mmclx_dir, datestr):
    """
    Scan every *vert.mmclx.gz file in *mmclx_dir* for the given *datestr*
    (YYYYMMDD) and return a sorted NumPy array of Unix timestamps (seconds)
    for rays whose elv variable is below ELV_FILL_THRESHOLD.
    """
    pattern = os.path.join(mmclx_dir, f'{datestr}_*.vert.mmclx.gz')
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f'No vert mmclx.gz files found matching: {pattern}')

    fill_unix_times = []
    n_files_with_fill = 0

    for fpath in files:
        with gzip.open(fpath) as gz:
            tmp = tempfile.NamedTemporaryFile(suffix='.nc', delete=False)
            tmp.write(gz.read())
            tmp.close()
        try:
            with nc4.Dataset(tmp.name) as ds:
                elv = np.array(ds.variables['elv'][:])
                t   = np.array(ds.variables['time'][:])   # Unix seconds
                bad = elv < ELV_FILL_THRESHOLD
                if np.any(bad):
                    fill_unix_times.extend(t[bad].tolist())
                    n_files_with_fill += 1
        finally:
            os.unlink(tmp.name)

    print(f'  {len(files)} vert mmclx files scanned; '
          f'{n_files_with_fill} contained elv fill values.')
    return np.array(sorted(fill_unix_times), dtype=float)


# ---------------------------------------------------------------------------
# Step 2: parse L1 time epoch
# ---------------------------------------------------------------------------

def l1_epoch_unix(time_units_str):
    """
    Parse a CF-convention time units string such as
    'seconds since 2024-12-16T00:00:00Z' and return the corresponding
    Unix timestamp (float seconds since 1970-01-01T00:00:00Z).
    """
    prefix = 'seconds since '
    if not time_units_str.startswith(prefix):
        raise ValueError(f'Unexpected time units: {time_units_str!r}')
    epoch_str = time_units_str[len(prefix):].strip()
    for fmt in ('%Y-%m-%dT%H:%M:%SZ', '%Y-%m-%dT%H:%M:%S',
                '%Y-%m-%d %H:%M:%S', '%Y-%m-%d'):
        try:
            dt = datetime.datetime.strptime(epoch_str, fmt).replace(
                tzinfo=datetime.timezone.utc)
            return dt.timestamp()
        except ValueError:
            continue
    raise ValueError(f'Cannot parse epoch from: {epoch_str!r}')


# ---------------------------------------------------------------------------
# Step 3: apply the fix
# ---------------------------------------------------------------------------

def apply_fix(l1_file, fill_unix_times, dry_run=False, no_backup=False):
    """
    Open *l1_file*, identify L1 rays that correspond to the mmclx fill-value
    period, set antenna_transition=1 for those rays, and update history.
    """
    with nc4.Dataset(l1_file) as ds:
        tu       = ds.variables['time'].units
        l1_times = np.array(ds.variables['time'][:], dtype=float)
        at_orig  = np.array(ds.variables['antenna_transition'][:])

    epoch_unix = l1_epoch_unix(tu)

    # Convert mmclx Unix timestamps to L1 time coordinates
    fill_l1 = fill_unix_times - epoch_unix

    t_min = fill_l1.min()
    t_max = fill_l1.max()

    now_utc = datetime.datetime.utcfromtimestamp

    print(f'  Fill-value period in L1 time coordinates: '
          f'{t_min:.0f} – {t_max:.0f} s  '
          f'({now_utc(fill_unix_times.min()).strftime("%H:%M:%S")} – '
          f'{now_utc(fill_unix_times.max()).strftime("%H:%M:%S")} UTC)')

    # Flag all L1 rays within [t_min - TOL, t_max + TOL]
    window = (l1_times >= t_min - TOLERANCE_S) & (l1_times <= t_max + TOLERANCE_S)

    n_window       = int(window.sum())
    n_already      = int(np.sum((at_orig == 1) & window))
    n_to_flag      = int(np.sum((at_orig == 0) & window))
    n_mmclx_rays   = len(fill_unix_times)

    print(f'  mmclx fill-value rays:        {n_mmclx_rays}')
    print(f'  L1 rays in window:            {n_window}')
    print(f'  already antenna_transition=1: {n_already}')
    print(f'  to be set to 1:               {n_to_flag}')

    if n_to_flag == 0:
        print('Nothing to change.')
        return

    if dry_run:
        print('[DRY RUN] No changes written to file.')
        return

    if not no_backup:
        backup = l1_file + '.bak'
        shutil.copy2(l1_file, backup)
        print(f'  Backup written: {backup}')

    with nc4.Dataset(l1_file, 'r+') as ds:
        at_new = np.array(ds.variables['antenna_transition'][:])
        at_new[window] = 1
        ds.variables['antenna_transition'][:] = at_new

        ts = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
        hist = getattr(ds, 'history', '').strip()
        entry = (
            f'{ts}: Set antenna_transition=1 for {n_window} rays '
            f'(L1 time {t_min:.0f}\u2013{t_max:.0f} s) where mmclx elv '
            f'contained fill value (-1000); fix applied by '
            f'fix_vpt_transitions_from_mmclx.py'
        )
        ds.history = f'{hist}\n{entry}' if hist else entry
        ds.last_modified = ts

    print(f'Done. {n_window} rays now have antenna_transition=1.')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            'Set antenna_transition=1 in a COBALT VPT L1 file for rays '
            'where the source mmclx elevation was a fill value (-1000).'
        )
    )
    ap.add_argument('-d', '--date', required=True,
                    help='Date string YYYYMMDD')
    ap.add_argument('--mmclx-dir', required=True,
                    help='Directory containing *vert.mmclx.gz files for DATE')
    ap.add_argument('--l1-file', required=True,
                    help='Full path to the VPT L1 NetCDF file to correct')
    ap.add_argument('--dry-run', action='store_true',
                    help='Report what would change but do not write anything')
    ap.add_argument('--no-backup', action='store_true',
                    help='Skip creating a .bak backup of the L1 file')
    return ap.parse_args()


def main():
    args = parse_args()

    print(f'Scanning mmclx fill-value rays for {args.date} ...')
    fill_times = scan_mmclx_fill_rays(args.mmclx_dir, args.date)

    if len(fill_times) == 0:
        print('No fill-value rays found in mmclx files. Nothing to do.')
        return

    print(f'Applying fix to: {args.l1_file}')
    apply_fix(
        args.l1_file,
        fill_times,
        dry_run=args.dry_run,
        no_backup=args.no_backup,
    )


if __name__ == '__main__':
    main()
