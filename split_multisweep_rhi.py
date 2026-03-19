#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
split_multisweep_rhi.py

Split CF-Radial RHI files that contain more than one sweep into individual
single-sweep files, one per sweep.

Each output file:
  - Contains a single sweep extracted from the source file
  - Is named using the sweep's own start time
  - Has time_coverage_start / time_coverage_end updated to match the sweep
  - Has phase_sequence updated to the single phase for that sweep

Usage:
    python split_multisweep_rhi.py -d DIRECTORY [-n] [--delete-originals]

Arguments:
    -d, --directory     Directory to search for RHI NetCDF files (searched
                        recursively for *_rhi_*.nc)
    -n, --dry-run       Print what would be done without writing any files
    --delete-originals  After successfully writing all split files, delete
                        the original multi-sweep files (default: keep them)

Author: Chris Walden, NCAS
"""

import getopt
import sys
import os
import glob
import datetime

import netCDF4 as nc4
import cftime
import numpy as np
import pyart


# ---------------------------------------------------------------------------
# Command-line parsing
# ---------------------------------------------------------------------------

def parse_args():
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "d:n",
                                ["directory=", "dry-run", "delete-originals"])
    except getopt.GetoptError as err:
        print(f"Error: {err}")
        print(__doc__)
        sys.exit(2)

    directory = None
    dry_run = False
    delete_originals = False

    for option, argument in opts:
        if option in ("-d", "--directory"):
            directory = argument
        elif option in ("-n", "--dry-run"):
            dry_run = True
        elif option == "--delete-originals":
            delete_originals = True

    if directory is None:
        print("Error: -d/--directory is required.")
        print(__doc__)
        sys.exit(2)

    if not os.path.isdir(directory):
        print(f"Error: directory does not exist: {directory}")
        sys.exit(2)

    return directory, dry_run, delete_originals


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_multisweep_rhi_files(directory):
    """Return list of *_rhi_*.nc files in directory tree that have >1 sweep."""
    all_rhi = sorted(glob.glob(os.path.join(directory, "**", "*_rhi_*.nc"),
                               recursive=True))
    multi = []
    for path in all_rhi:
        try:
            with nc4.Dataset(path, 'r') as ds:
                nsweeps = len(ds.variables['sweep_number'])
            if nsweeps > 1:
                multi.append((path, nsweeps))
        except Exception as e:
            print(f"  Warning: could not read {path}: {e}")
    return multi


def build_output_filename(source_path, sweep_start_time):
    """
    Derive an output filename from source_path by substituting the sweep
    start timestamp.

    Source pattern:
        <instrument>_cao_YYYYMMDD-HHMMSS_rhi_l1_<version>.nc

    Output:
        <instrument>_cao_YYYYMMDD-HHMMSS_rhi_l1_<version>.nc
    using the sweep's own start time.
    """
    basename = os.path.basename(source_path)
    parts = basename.split('_')

    # Find the position of the timestamp field (YYYYMMDD-HHMMSS)
    ts_idx = None
    for i, p in enumerate(parts):
        if len(p) == 15 and p[8] == '-' and p[:8].isdigit() and p[9:].isdigit():
            ts_idx = i
            break

    if ts_idx is None:
        # Fallback: insert timestamp before 'rhi'
        new_ts = sweep_start_time.strftime('%Y%m%d-%H%M%S')
        new_basename = basename.replace('_rhi_', f'_{new_ts}_rhi_')
    else:
        new_ts = sweep_start_time.strftime('%Y%m%d-%H%M%S')
        parts[ts_idx] = new_ts
        new_basename = '_'.join(parts)

    return os.path.join(os.path.dirname(source_path), new_basename)


def update_netcdf_time_metadata(path, sweep_start, sweep_end, phase_str=None):
    """
    Fix up time_coverage_start / time_coverage_end global attributes and
    string variables in-place after pyart has written the file.

    Also updates phase_sequence if phase_str is supplied.
    """
    fmt = "%Y-%m-%dT%H:%M:%SZ"
    t_start_str = sweep_start.strftime(fmt)
    t_end_str   = sweep_end.strftime(fmt)

    with nc4.Dataset(path, 'a') as ds:
        # Global attributes
        ds.setncattr('time_coverage_start', t_start_str)
        ds.setncattr('time_coverage_end',   t_end_str)
        if phase_str is not None:
            ds.setncattr('phase_sequence', phase_str)

        # Character array variables (CF-Radial convention)
        for var_name, value in [('time_coverage_start', t_start_str),
                                 ('time_coverage_end',   t_end_str)]:
            if var_name in ds.variables:
                var = ds.variables[var_name]
                encoded = np.array(list(value), dtype='S1')
                # Pad or truncate to match existing variable length
                n = var.shape[0]
                if len(encoded) < n:
                    pad = np.array([b''] * (n - len(encoded)), dtype='S1')
                    encoded = np.concatenate([encoded, pad])
                else:
                    encoded = encoded[:n]
                var[:] = encoded


# ---------------------------------------------------------------------------
# Main splitting logic
# ---------------------------------------------------------------------------

def split_file(source_path, nsweeps, dry_run=False):
    """
    Split source_path into nsweeps single-sweep files.

    Returns list of output paths written (empty in dry-run mode).
    """
    print(f"\n  Source ({nsweeps} sweeps): {os.path.basename(source_path)}")

    # Read phase_sequence from source file before loading with pyart
    with nc4.Dataset(source_path, 'r') as ds:
        phase_sequence = None
        if 'phase_sequence' in ds.ncattrs():
            raw = ds.getncattr('phase_sequence')
            # May be comma-separated string or a single value
            phase_sequence = [p.strip() for p in str(raw).split(',')]

    radar = pyart.io.read_cfradial(source_path)

    output_paths = []

    for sweep_idx in range(nsweeps):
        # Extract single sweep
        sweep_radar = radar.extract_sweeps([sweep_idx])

        # Determine sweep time range
        t_start = cftime.num2pydate(sweep_radar.time['data'][0],
                                    sweep_radar.time['units'])
        t_end   = cftime.num2pydate(sweep_radar.time['data'][-1],
                                    sweep_radar.time['units'])

        # Update metadata dict so pyart writes correct global attrs
        sweep_radar.metadata['time_coverage_start'] = \
            t_start.strftime("%Y-%m-%dT%H:%M:%SZ")
        sweep_radar.metadata['time_coverage_end'] = \
            t_end.strftime("%Y-%m-%dT%H:%M:%SZ")

        # Single-phase phase_sequence for this sweep
        sweep_phase = None
        if phase_sequence is not None:
            if sweep_idx < len(phase_sequence):
                sweep_phase = phase_sequence[sweep_idx]
                sweep_radar.metadata['phase_sequence'] = sweep_phase
            else:
                sweep_phase = phase_sequence[-1]
                sweep_radar.metadata['phase_sequence'] = sweep_phase

        # Build output filename
        out_path = build_output_filename(source_path, t_start)

        fixed_angle = float(sweep_radar.fixed_angle['data'][0])
        print(f"    Sweep {sweep_idx}: az={fixed_angle:.2f}°  "
              f"{t_start.strftime('%H:%M:%S')} -> {t_end.strftime('%H:%M:%S')}  "
              f"-> {os.path.basename(out_path)}")

        if dry_run:
            continue

        if os.path.exists(out_path):
            print(f"      Warning: output already exists, overwriting: {out_path}")

        # Cast sweep index variables to int32 - CF-Radial spec requires 'int'
        # (NC_INT / 32-bit); pyart's extract_sweeps can produce int64.
        for key in ('sweep_number', 'sweep_start_ray_index', 'sweep_end_ray_index'):
            if hasattr(sweep_radar, key):
                obj = getattr(sweep_radar, key)
                if 'data' in obj:
                    obj['data'] = obj['data'].astype(np.int32)

        # Write file
        pyart.io.write_cfradial(out_path, sweep_radar)

        # Fix up time metadata that pyart may not update correctly
        update_netcdf_time_metadata(out_path, t_start, t_end,
                                    phase_str=sweep_phase)

        output_paths.append(out_path)
        print(f"      Written.")

    return output_paths


def main():
    directory, dry_run, delete_originals = parse_args()

    if dry_run:
        print("DRY RUN — no files will be written or deleted.\n")

    print(f"Searching for multi-sweep RHI files in: {directory}")
    multi_files = find_multisweep_rhi_files(directory)

    if not multi_files:
        print("No multi-sweep RHI files found.")
        return

    print(f"Found {len(multi_files)} multi-sweep RHI file(s):")
    for path, n in multi_files:
        print(f"  {n} sweeps  {path}")

    all_written = []
    sources_to_delete = []

    for source_path, nsweeps in multi_files:
        written = split_file(source_path, nsweeps, dry_run=dry_run)
        if len(written) == nsweeps:
            all_written.extend(written)
            sources_to_delete.append(source_path)
        elif not dry_run:
            print(f"  Warning: only {len(written)}/{nsweeps} sweeps written "
                  f"for {source_path} — original will not be deleted.")

    if not dry_run:
        print(f"\nWrote {len(all_written)} file(s) total.")

        if delete_originals and sources_to_delete:
            print(f"Deleting {len(sources_to_delete)} original file(s):")
            for path in sources_to_delete:
                os.remove(path)
                print(f"  Deleted: {path}")
        elif not delete_originals:
            print("Original multi-sweep files kept "
                  "(use --delete-originals to remove them).")
    else:
        print(f"\nDry run complete. {len(multi_files)} file(s) would be split.")


if __name__ == '__main__':
    main()
