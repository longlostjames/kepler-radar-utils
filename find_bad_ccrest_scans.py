#!/usr/bin/env python
"""
find_bad_ccrest_scans.py
------------------------
Scan CCREST-M L1 NetCDF files for a given day and flag files / rays that
contain known data quality problems:

  1. Sentinel az or el values (< --sentinel-threshold, default -900)
  2. Files where ALL rays are sentinel (entire file unusable)
  3. Files with very few rays (< --min-rays, default 15)
  4. [optional] Az or el disagreement with the axis log beyond
     --az-tol / --el-tol after removing the systematic az offset

Results are printed as a table and optionally written to a CSV.

Usage
-----
    python find_bad_ccrest_scans.py YYYYMMDD L1_PATH [options]

    # With axis-log cross-check:
    python find_bad_ccrest_scans.py YYYYMMDD L1_PATH \\
        --logp-path /path/to/logp \\
        --az-tol 5 --el-tol 2

Positional arguments
--------------------
    YYYYMMDD   Date string, e.g. 20240307
    L1_PATH    Root L1 path; NC files expected under L1_PATH/YYYYMMDD/*.nc

Optional arguments
------------------
    --logp-path PATH        Directory containing *.axis.gz log files.
                            If supplied, enables axis-log cross-check.
    --sentinel-threshold T  Az or el values below this are sentinel/fill
                            (default: -900)
    --min-rays N            Files with fewer than N rays are flagged as short
                            (default: 15)
    --az-tol DEGREES        Max |corrected az residual| vs log to flag
                            (default: 5, only used with --logp-path)
    --el-tol DEGREES        Max |el residual| vs log to flag
                            (default: 2, only used with --logp-path)
    --output CSV            Write results to this CSV file
    --all                   Also report files with no issues (default: bad only)
"""

import argparse
import csv
import datetime
import glob
import gzip
import os
import sys

import numpy as np

try:
    import netCDF4 as nc4
    import cftime
except ImportError:
    sys.exit('Error: netCDF4 is required.  conda install netCDF4')


# ---------------------------------------------------------------------------
# Axis log helpers
# ---------------------------------------------------------------------------

def _load_axis_log(logp_dir, datestr):
    """Load all log entries for *datestr*. Returns (times, az, el) arrays."""
    day_start = datetime.datetime.strptime(datestr, '%Y%m%d').replace(
        tzinfo=datetime.timezone.utc).timestamp()
    day_end = day_start + 86400.0
    day_start_dt = datetime.datetime.utcfromtimestamp(day_start).replace(
        tzinfo=datetime.timezone.utc)
    day_end_dt   = datetime.datetime.utcfromtimestamp(day_end).replace(
        tzinfo=datetime.timezone.utc)

    all_gz = sorted(glob.glob(os.path.join(logp_dir, '*.axis.gz')))
    matching = []
    for path in all_gz:
        stem = os.path.basename(path).replace('.axis.gz', '')
        try:
            fdt = datetime.datetime.strptime(stem, '%Y%m%d%H%M').replace(
                tzinfo=datetime.timezone.utc)
        except ValueError:
            continue
        if fdt < day_end_dt and fdt >= day_start_dt - datetime.timedelta(hours=1):
            matching.append(path)

    if not matching:
        return None, None, None

    times_list, az_list, el_list = [], [], []
    for path in matching:
        try:
            with gzip.open(path, 'rt') as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    t = float(parts[0])
                    if t < day_start or t >= day_end:
                        continue
                    times_list.append(t)
                    az_list.append(float(parts[1]))
                    el_list.append(float(parts[2]))
        except Exception as exc:
            print(f'  Warning: could not read {os.path.basename(path)}: {exc}')

    if not times_list:
        return None, None, None

    times = np.array(times_list)
    az    = np.array(az_list)
    el    = np.array(el_list)
    order = np.argsort(times)
    return times[order], az[order], el[order]


def _az_diff(a, b):
    """Signed azimuth difference a − b wrapped to (−180, +180]."""
    return (a - b + 180.0) % 360.0 - 180.0


def _interp_log(log_t, log_az, log_el, query_unix):
    """Interpolate log az/el (with unwrapped az) to query_unix times."""
    log_az_uw = np.unwrap(np.deg2rad(log_az))
    interp_az = np.rad2deg(np.interp(query_unix, log_t, log_az_uw)) % 360.0
    interp_el = np.interp(query_unix, log_t, log_el)
    return interp_az, interp_el


# ---------------------------------------------------------------------------
# Per-file check
# ---------------------------------------------------------------------------

def check_file(ncfile, sentinel_threshold, min_rays,
               log_t=None, log_az=None, log_el=None,
               az_tol=5.0, el_tol=2.0):
    """
    Return a dict with quality information for one NC file.

    Keys
    ----
    filename, scan_type, n_rays,
    n_sentinel_az, n_sentinel_el, n_sentinel_any,
    all_sentinel,
    short_file,
    log_az_offset,          # median NC−log az residual (non-sentinel)
    n_az_disagree,          # rays where corrected |Δaz| > az_tol
    n_el_disagree,          # rays where |Δel| > el_tol (non-sentinel)
    log_el_at_sentinel,     # median log el during sentinel rays (if log available)
    flags,                  # list of plain-English flag strings
    bad                     # True if any flag is set
    """
    basename = os.path.basename(ncfile)
    result = dict(
        filename=basename,
        scan_type='unknown',
        n_rays=0,
        n_sentinel_az=0, n_sentinel_el=0, n_sentinel_any=0,
        all_sentinel=False,
        short_file=False,
        log_az_offset=float('nan'),
        n_az_disagree=0, n_el_disagree=0,
        log_el_at_sentinel='',
        flags=[], bad=False,
    )

    # Infer scan type from filename
    b = basename.lower()
    for t in ('ppi', 'rhi', 'vpt', 'pointing'):
        if f'_{t}_' in b:
            result['scan_type'] = t
            break

    try:
        with nc4.Dataset(ncfile) as ds:
            t_data  = ds['time'][:]
            t_units = ds['time'].units
            t_cal   = getattr(ds['time'], 'calendar', 'standard')
            nc_az   = ds['azimuth'][:].data.astype(float)
            nc_el   = ds['elevation'][:].data.astype(float)
    except Exception as exc:
        result['flags'].append(f'read_error: {exc}')
        result['bad'] = True
        return result

    n = len(nc_az)
    result['n_rays'] = n

    # --- Short file check ---
    if n < min_rays:
        result['short_file'] = True
        result['flags'].append(f'short_file ({n} rays < {min_rays})')

    # --- Sentinel checks ---
    sent_az  = nc_az < sentinel_threshold
    sent_el  = nc_el < sentinel_threshold
    sent_any = sent_az | sent_el

    result['n_sentinel_az']  = int(sent_az.sum())
    result['n_sentinel_el']  = int(sent_el.sum())
    result['n_sentinel_any'] = int(sent_any.sum())

    if sent_az.any():
        pct = 100.0 * sent_az.sum() / n
        result['flags'].append(
            f'sentinel_az ({sent_az.sum()}/{n} rays, {pct:.0f}%)')
    if sent_el.any():
        pct = 100.0 * sent_el.sum() / n
        result['flags'].append(
            f'sentinel_el ({sent_el.sum()}/{n} rays, {pct:.0f}%)')
    if sent_any.all():
        result['all_sentinel'] = True
        result['flags'].append('all_sentinel')

    # --- Axis log cross-check ---
    if log_t is not None:
        try:
            dts  = cftime.num2pydate(t_data, t_units, calendar=t_cal)
            unix = np.array([d.replace(tzinfo=datetime.timezone.utc).timestamp()
                             for d in dts])
        except Exception as exc:
            result['flags'].append(f'time_conv_error: {exc}')
            result['bad'] = bool(result['flags'])
            return result

        in_range = (unix >= log_t[0]) & (unix <= log_t[-1])
        good = in_range & ~sent_any

        if good.any():
            iaz, iel = _interp_log(log_t, log_az, log_el, unix[good])
            daz = _az_diff(nc_az[good], iaz)
            # Auto-detect and store the systematic az offset
            az_offset = float(np.median(daz))
            result['log_az_offset'] = az_offset
            daz_corr = daz - az_offset
            del_corr = nc_el[good] - iel

            n_az_dis = int((np.abs(daz_corr) > az_tol).sum())
            n_el_dis = int((np.abs(del_corr) > el_tol).sum())
            result['n_az_disagree'] = n_az_dis
            result['n_el_disagree'] = n_el_dis
            if n_az_dis:
                result['flags'].append(
                    f'az_disagree_with_log ({n_az_dis} rays, tol {az_tol}°)')
            if n_el_dis:
                result['flags'].append(
                    f'el_disagree_with_log ({n_el_dis} rays, tol {el_tol}°)')

        # What was the log el doing during sentinel rays?
        sent_in_range = in_range & sent_any
        if sent_in_range.any():
            _, iel_sent = _interp_log(log_t, log_az, log_el, unix[sent_in_range])
            result['log_el_at_sentinel'] = f'{np.median(iel_sent):.1f}'

    result['bad'] = bool(result['flags'])
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Find bad CCREST-M L1 scans for a given day.')
    parser.add_argument('datestr',  help='Date string YYYYMMDD')
    parser.add_argument('l1_path',  help='Root L1 path (NC under L1_PATH/YYYYMMDD/)')
    parser.add_argument('--logp-path',          default=None,
                        help='Directory with *.axis.gz log files (enables log cross-check)')
    parser.add_argument('--sentinel-threshold', type=float, default=-900.0,
                        help='Az or el values below this are sentinel/fill (default: -900)')
    parser.add_argument('--min-rays',           type=int,   default=15,
                        help='Files with fewer rays than this are flagged (default: 15)')
    parser.add_argument('--az-tol',             type=float, default=5.0,
                        help='Corrected az residual threshold vs log in deg (default: 5)')
    parser.add_argument('--el-tol',             type=float, default=2.0,
                        help='El residual threshold vs log in deg (default: 2)')
    parser.add_argument('--output',             default=None,
                        help='Write results to this CSV file')
    parser.add_argument('--all',  action='store_true',
                        help='Report all files, not just bad ones')
    args = parser.parse_args()

    inpath_date = os.path.join(args.l1_path, args.datestr)
    if not os.path.isdir(inpath_date):
        sys.exit(f'Directory not found: {inpath_date}')

    nc_files = sorted(glob.glob(os.path.join(inpath_date, '*.nc')))
    if not nc_files:
        sys.exit(f'No NC files found in {inpath_date}')

    # Load axis log if requested
    log_t = log_az = log_el = None
    if args.logp_path:
        print(f'Loading axis log for {args.datestr}...')
        log_t, log_az, log_el = _load_axis_log(args.logp_path, args.datestr)
        if log_t is None:
            print('  Warning: no log data found; skipping log cross-check')
        else:
            print(f'  Loaded {len(log_t):,} log samples')

    print(f'\nChecking {len(nc_files)} NC files...\n')

    results = []
    for ncfile in nc_files:
        r = check_file(
            ncfile,
            sentinel_threshold=args.sentinel_threshold,
            min_rays=args.min_rays,
            log_t=log_t, log_az=log_az, log_el=log_el,
            az_tol=args.az_tol, el_tol=args.el_tol,
        )
        results.append(r)

    # --- Summary ---
    bad_results = [r for r in results if r['bad']]
    print(f'Files checked : {len(results)}')
    print(f'Files flagged : {len(bad_results)}')

    # Counts by flag type
    flag_types = ['short_file', 'sentinel_az', 'sentinel_el', 'all_sentinel',
                  'az_disagree', 'el_disagree']
    for ft in flag_types:
        n = sum(1 for r in results if any(ft in fl for fl in r['flags']))
        if n:
            print(f'  {ft:30s}: {n}')

    # --- Print table ---
    to_print = results if args.all else bad_results
    if to_print:
        print(f'\n{"Filename":<65s} {"Type":8s} {"Rays":>6s} '
              f'{"SentAz":>7s} {"SentEl":>7s} {"LogEl@Sent":>11s}  Flags')
        print('-' * 130)
        for r in to_print:
            print(
                f'{r["filename"]:<65s} {r["scan_type"]:8s} {r["n_rays"]:>6d} '
                f'{r["n_sentinel_az"]:>7d} {r["n_sentinel_el"]:>7d} '
                f'{r["log_el_at_sentinel"]:>11s}  '
                + ', '.join(r["flags"])
            )
    else:
        print('\nNo bad files found.')

    # --- CSV output ---
    if args.output:
        fieldnames = [
            'filename', 'scan_type', 'n_rays',
            'n_sentinel_az', 'n_sentinel_el', 'n_sentinel_any',
            'all_sentinel', 'short_file',
            'log_az_offset', 'n_az_disagree', 'n_el_disagree',
            'log_el_at_sentinel', 'flags',
        ]
        with open(args.output, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            for r in results if args.all else bad_results:
                row = dict(r)
                row['flags'] = '; '.join(r['flags'])
                writer.writerow(row)
        print(f'\nResults written to {args.output}')


if __name__ == '__main__':
    main()
