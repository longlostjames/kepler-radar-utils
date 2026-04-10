#!/usr/bin/env python
"""
compare_axis_log_vs_netcdf.py
------------------------------
Compare antenna azimuth / elevation from the raw axis log files against the
values stored in the L1 NetCDF files for a given day.

The axis logs are sampled at ~2 Hz; the NetCDF files contain one value per
radar ray (may be much faster). For each NetCDF ray the log az/el is linearly
interpolated to the ray timestamp and the residual (NC – log) is computed.

Sentinel rays (az=0, el=90 or any other physically suspicious values written
during antenna parking / errors) are highlighted in the output plot.

Usage
-----
    python compare_axis_log_vs_netcdf.py YYYYMMDD L1_PATH LOGP_PATH [options]

Positional arguments
--------------------
    YYYYMMDD    Date string, e.g. 20240203
    L1_PATH     Root L1 path; NC files are expected under L1_PATH/YYYYMMDD/*.nc
    LOGP_PATH   Directory containing *.axis.gz log files

Optional arguments
------------------
    --output OUTPUT.png    Output plot file (default: axis_vs_nc_YYYYMMDD.png)
    --downsample N         Keep 1-in-N NC rays for speed (default: 10)
    --az-offset DEGREES    Known systematic azimuth offset to subtract from
                           (NC_az − log_az) before computing residuals
                           (default: auto-detect from median residual)
    --sentinel-threshold T Az or el values below this threshold are treated as
                           sentinel / fill values (default: -900)
"""

import argparse
import glob
import gzip
import os
import sys
import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

try:
    import netCDF4 as nc4
    import cftime
except ImportError:
    sys.exit('Error: netCDF4 is required.  conda install netCDF4')

# ---------------------------------------------------------------------------
# Scan-type colour scheme
# ---------------------------------------------------------------------------
SCAN_COLORS = {
    'ppi':      '#e6194b',   # red
    'rhi':      '#3cb44b',   # green
    'vpt':      '#4363d8',   # blue
    'pointing': '#f58231',   # orange
    'other':    '#a9a9a9',   # gray
}


# ---------------------------------------------------------------------------
# Log loading
# ---------------------------------------------------------------------------

def _log_files_for_day(logp_dir: str, datestr: str):
    """Return all axis log files whose name starts with YYYYMMDD or that
    cover the day (files are named YYYYMMDDHHMM.axis.gz)."""
    # The midnight file covers the whole day, but there may also be a file
    # from the previous day that rolls over into this one.
    all_gz = sorted(glob.glob(os.path.join(logp_dir, '*.axis.gz')))
    if not all_gz:
        return []

    day_start = datetime.datetime.strptime(datestr, '%Y%m%d').replace(
        tzinfo=datetime.timezone.utc)
    day_end = day_start + datetime.timedelta(days=1)

    matching = []
    for path in all_gz:
        stem = os.path.basename(path).replace('.axis.gz', '')
        try:
            file_dt = datetime.datetime.strptime(stem, '%Y%m%d%H%M').replace(
                tzinfo=datetime.timezone.utc)
        except ValueError:
            continue
        # Include if the file starts within the day *or* starts on the
        # previous day (could span midnight) — we will filter by timestamp
        # when loading.
        if file_dt < day_end and file_dt >= day_start - datetime.timedelta(hours=1):
            matching.append(path)

    return matching


def load_axis_logs(logp_dir: str, datestr: str):
    """Load all log entries for *datestr* from the axis log directory.

    Returns arrays (unix_times, azimuths, elevations), all float64, sorted
    by ascending unix time.  Only rows whose timestamp falls within the UTC
    day are kept.
    """
    day_start = datetime.datetime.strptime(datestr, '%Y%m%d').replace(
        tzinfo=datetime.timezone.utc).timestamp()
    day_end = day_start + 86400.0

    log_files = _log_files_for_day(logp_dir, datestr)
    if not log_files:
        sys.exit(f'No axis log files found for {datestr} in {logp_dir}')

    print(f'Loading {len(log_files)} axis log file(s)...')

    times_list, az_list, el_list = [], [], []
    for path in log_files:
        print(f'  {os.path.basename(path)}')
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
            print(f'  WARNING: could not read {os.path.basename(path)}: {exc}')

    if not times_list:
        sys.exit(f'No log entries found within UTC day {datestr}')

    times = np.array(times_list, dtype=np.float64)
    az    = np.array(az_list,    dtype=np.float64)
    el    = np.array(el_list,    dtype=np.float64)

    # Sort by time
    order = np.argsort(times)
    print(f'Loaded {len(times):,} log samples')
    return times[order], az[order], el[order]


# ---------------------------------------------------------------------------
# NetCDF loading
# ---------------------------------------------------------------------------

def _scan_type_from_filename(filename: str) -> str:
    b = os.path.basename(filename).lower()
    for t in ('ppi', 'rhi', 'vpt', 'pointing'):
        if f'_{t}_' in b:
            return t
    return 'other'


def load_netcdf_rays(l1_path: str, datestr: str, downsample: int = 10):
    """Load az/el/time from all NC files for the day.

    Returns:
        unix_times  : float64 array, seconds since epoch
        azimuths    : float64 array, degrees
        elevations  : float64 array, degrees
        scan_types  : list of str, one per ray
        filenames   : list of str, one per ray (basename without path)
    """
    inpath_date = os.path.join(l1_path, datestr)
    if not os.path.isdir(inpath_date):
        sys.exit(f'L1 directory not found: {inpath_date}')

    nc_files = sorted(glob.glob(os.path.join(inpath_date, '*.nc')))
    if not nc_files:
        sys.exit(f'No NC files found in {inpath_date}')

    print(f'Loading {len(nc_files)} NC files (downsample 1-in-{downsample})...')

    all_times, all_az, all_el, all_types, all_files = [], [], [], [], []

    for ncfile in nc_files:
        scan_type = _scan_type_from_filename(ncfile)
        try:
            with nc4.Dataset(ncfile) as ds:
                t_data   = ds['time'][:]
                t_units  = ds['time'].units
                t_cal    = getattr(ds['time'], 'calendar', 'standard')
                az_data  = ds['azimuth'][:]
                el_data  = ds['elevation'][:]
        except Exception as exc:
            print(f'  WARNING: {os.path.basename(ncfile)}: {exc}')
            continue

        try:
            py_dates = cftime.num2pydate(t_data, t_units, calendar=t_cal)
            unix = np.array([d.replace(tzinfo=datetime.timezone.utc).timestamp()
                             for d in py_dates], dtype=np.float64)
        except Exception as exc:
            print(f'  WARNING: time conversion failed for {os.path.basename(ncfile)}: {exc}')
            continue

        # Downsample
        idx = np.arange(0, len(unix), downsample)
        all_times.append(unix[idx])
        all_az.append(np.asarray(az_data[idx], dtype=np.float64))
        all_el.append(np.asarray(el_data[idx], dtype=np.float64))
        all_types.extend([scan_type] * len(idx))
        all_files.extend([os.path.basename(ncfile)] * len(idx))

    unix_times = np.concatenate(all_times)
    azimuths   = np.concatenate(all_az)
    elevations = np.concatenate(all_el)
    order      = np.argsort(unix_times)
    print(f'Loaded {len(unix_times):,} NC rays (after downsampling)')
    return (unix_times[order], azimuths[order], elevations[order],
            [all_types[i] for i in order],
            [all_files[i] for i in order])


# ---------------------------------------------------------------------------
# Residual computation with azimuth wrap handling
# ---------------------------------------------------------------------------

def _az_diff(a, b):
    """Signed difference a − b in (−180, +180]."""
    d = (a - b + 180.0) % 360.0 - 180.0
    return d


def interpolate_log(log_times, log_az, log_el, query_times):
    """Linearly interpolate log az/el to *query_times*.

    Azimuth interpolation uses unwrapped values to handle the 0/360 boundary.
    Returns (interp_az, interp_el).
    """
    log_az_unwrap = np.unwrap(np.deg2rad(log_az))
    interp_az_rad = np.interp(query_times, log_times, log_az_unwrap)
    interp_az = np.rad2deg(interp_az_rad) % 360.0

    interp_el = np.interp(query_times, log_times, log_el)

    return interp_az, interp_el


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _per_scan_type_stats(az_resid_raw, el_resid, scan_types):
    """Return a list of (stype, n, mean_daz, std_daz, mean_del, std_del) tuples."""
    rows = []
    for stype in ('ppi', 'rhi', 'vpt', 'pointing', 'other'):
        mask = np.array([s == stype for s in scan_types])
        if not np.any(mask):
            continue
        da = az_resid_raw[mask]
        de = el_resid[mask]
        rows.append((stype, int(mask.sum()),
                     float(da.mean()), float(da.std()),
                     float(de.mean()), float(de.std())))
    return rows


def make_plot(datestr, log_times, log_az, log_el,
              nc_times, nc_az, nc_el,
              interp_az, interp_el,
              scan_types, filenames,
              sentinel_threshold, az_offset,
              output_path):

    dt_log = np.array([datetime.datetime.utcfromtimestamp(t) for t in log_times])
    dt_nc  = np.array([datetime.datetime.utcfromtimestamp(t) for t in nc_times])

    # Sentinel rays: NC az or el value below the sentinel threshold (fill values)
    az_sentinel = nc_az < sentinel_threshold
    el_sentinel = nc_el < sentinel_threshold
    sentinel = az_sentinel | el_sentinel
    n_az_sent = int(np.sum(az_sentinel))
    n_el_sent = int(np.sum(el_sentinel))
    n_sent    = int(np.sum(sentinel))
    print(f'Sentinel rays (NC az < {sentinel_threshold}): {n_az_sent}')
    print(f'Sentinel rays (NC el < {sentinel_threshold}): {n_el_sent}')

    # Compute residuals only for non-sentinel rays
    good = ~sentinel
    az_resid_raw = np.full(len(nc_az), np.nan)
    el_resid     = np.full(len(nc_el), np.nan)
    az_resid_raw[good] = _az_diff(nc_az[good], interp_az[good])
    el_resid[good]     = nc_el[good] - interp_el[good]

    # Per-scan-type stats (non-sentinel rays, raw residual so the offset is visible)
    stats_rows = _per_scan_type_stats(az_resid_raw, el_resid, scan_types)

    print(f'\nPer-scan-type statistics (NC − log, non-sentinel rays, before offset correction):')
    print(f'{"Scan":10s} {"Rays":>7s}  {"MeanΔaz":>8s} {"StdΔaz":>8s}  {"MeanΔel":>8s} {"StdΔel":>7s}')
    print('-' * 62)
    for stype, n, mda, sda, mde, sde in stats_rows:
        print(f'{stype:10s} {n:>7d}  {mda:>+8.2f}° {sda:>8.2f}°  {mde:>+8.2f}° {sde:>7.2f}°')

    # Auto-detect az offset from median of non-sentinel residuals
    valid_daz = az_resid_raw[good & np.isfinite(az_resid_raw)]
    if az_offset is None:
        az_offset = float(np.median(valid_daz)) if len(valid_daz) else 0.0
        print(f'\nAuto-detected az offset (median NC−log, non-sentinel): {az_offset:+.2f}°')
    else:
        print(f'\nUsing supplied az offset: {az_offset:+.2f}°')

    az_resid = az_resid_raw - az_offset   # corrected residual (NaN for sentinels)

    fig, axes = plt.subplots(4, 1, figsize=(16, 14), sharex=True,
                             gridspec_kw={'height_ratios': [2, 2, 1.2, 1.2]})
    fig.suptitle(
        f'Axis Log vs NetCDF Comparison \u2014 {datestr[:4]}-{datestr[4:6]}-{datestr[6:]}\n'
        f'Log {len(log_times):,} samples  |  NC {len(nc_times):,} rays (downsampled)  '
        f'|  Sentinel threshold: <{sentinel_threshold}  '
        f'|  Az offset applied: {az_offset:+.2f}\u00b0',
        fontsize=12, y=0.99,
    )

    date_fmt = mdates.DateFormatter('%H:%M')
    xmin = datetime.datetime(int(datestr[:4]), int(datestr[4:6]), int(datestr[6:]),  0, 0)
    xmax = xmin + datetime.timedelta(days=1)

    # ---- Azimuth panel ----
    ax0 = axes[0]
    ax0.plot(dt_log, log_az, color='silver', lw=0.5, zorder=1, label='Axis log')
    for stype, colour in SCAN_COLORS.items():
        mask = np.array([s == stype for s in scan_types])
        if not np.any(mask):
            continue
        ax0.scatter(dt_nc[mask], nc_az[mask], color=colour, s=2, alpha=0.6,
                    zorder=2, label=stype)
    # Highlight sentinel rays (value < threshold)
    if n_sent:
        ax0.scatter(dt_nc[sentinel], nc_az[sentinel],
                    marker='x', color='black', s=20, zorder=3,
                    label=f'Sentinel (<{sentinel_threshold}): {n_sent}')
    ax0.set_ylabel('Azimuth (°)')
    ax0.set_ylim(-5, 365)
    ax0.legend(loc='upper right', fontsize=7, markerscale=2)
    ax0.set_title('Azimuth', fontsize=10)
    ax0.grid(axis='x', lw=0.3, color='#dddddd')

    # ---- Elevation panel ----
    ax1 = axes[1]
    ax1.plot(dt_log, log_el, color='silver', lw=0.5, zorder=1, label='Axis log')
    for stype, colour in SCAN_COLORS.items():
        mask = np.array([s == stype for s in scan_types])
        if not np.any(mask):
            continue
        ax1.scatter(dt_nc[mask], nc_el[mask], color=colour, s=2, alpha=0.6,
                    zorder=2, label=stype)
    if n_sent:
        ax1.scatter(dt_nc[sentinel], nc_el[sentinel],
                    marker='x', color='black', s=20, zorder=3,
                    label=f'Sentinel (<{sentinel_threshold}): {n_sent}')
    ax1.set_ylabel('Elevation (°)')
    ax1.legend(loc='upper right', fontsize=7, markerscale=2)
    ax1.set_title('Elevation', fontsize=10)
    ax1.grid(axis='x', lw=0.3, color='#dddddd')

    # ---- Az residual panel (offset-corrected, NaN for sentinel rays) ----
    ax2 = axes[2]
    ax2.axhline(0, color='silver', lw=0.8)
    good_mask = np.isfinite(az_resid)
    ax2.scatter(dt_nc[good_mask], az_resid[good_mask], c='steelblue', s=1, alpha=0.4)
    if n_sent:
        # Plot sentinel rays at y=0 as crosses to show when they occur
        ax2.scatter(dt_nc[sentinel], np.zeros(n_sent),
                    marker='x', color='red', s=15, zorder=3,
                    label=f'Sentinel rays: {n_sent}')
    ax2.set_ylabel(f'\u0394 Az (NC\u2212log\u2212{az_offset:.1f}\u00b0)')
    ax2.legend(loc='upper right', fontsize=7)
    ax2.grid(axis='x', lw=0.3, color='#dddddd')

    # ---- Statistics table inset ----
    col_labels = ['Scan', 'Rays', 'Mean\u0394az', 'Std\u0394az', 'Mean\u0394el', 'Std\u0394el']
    table_data = [[s, str(n), f'{mda:+.2f}\u00b0', f'{sda:.2f}\u00b0', f'{mde:+.2f}\u00b0', f'{sde:.2f}\u00b0']
                  for s, n, mda, sda, mde, sde in stats_rows]
    tbl = ax2.table(cellText=table_data, colLabels=col_labels,
                    loc='upper left', cellLoc='right', bbox=[0.0, -0.05, 0.38, 0.95])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7)

    # ---- El residual panel ----
    ax3 = axes[3]
    ax3.axhline(0, color='silver', lw=0.8)
    good_el_mask = np.isfinite(el_resid)
    ax3.scatter(dt_nc[good_el_mask], el_resid[good_el_mask], c='steelblue', s=1, alpha=0.4)
    if n_sent:
        ax3.scatter(dt_nc[sentinel], np.zeros(n_sent),
                    marker='x', color='red', s=15, zorder=3,
                    label=f'Sentinel rays: {n_sent}')
    ax3.set_ylabel('\u0394 El (NC\u2212log) \u00b0')
    ax3.legend(loc='upper right', fontsize=7)
    ax3.xaxis.set_major_formatter(date_fmt)
    ax3.set_xlabel(f'Time UTC ({datestr[:4]}-{datestr[4:6]}-{datestr[6:]})')
    ax3.grid(axis='x', lw=0.3, color='#dddddd')

    for ax in axes:
        ax.set_xlim(xmin, xmax)

    fig.autofmt_xdate(rotation=30, ha='right')
    fig.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'\nPlot saved: {output_path}')
    print(f'Total sentinel rays (NC az or el < {sentinel_threshold}): {n_sent}')
    if n_sent:
        print('\nSentinel rays (first 20):')
        sent_idx = np.where(sentinel)[0]
        for i in sent_idx[:20]:
            print(
                f'  {dt_nc[i].strftime("%H:%M:%S")}  {filenames[i]}'
                f'  NC az={nc_az[i]:.2f}  el={nc_el[i]:.2f}'
            )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compare axis log az/el vs NetCDF ray az/el for a given day.'
    )
    parser.add_argument('datestr',  help='Date string YYYYMMDD')
    parser.add_argument('l1_path',  help='Root L1 path (NC files under L1_PATH/YYYYMMDD/)')
    parser.add_argument('logp_path', help='Directory containing *.axis.gz log files')
    parser.add_argument('--output', default=None,
                        help='Output PNG file (default: axis_vs_nc_YYYYMMDD.png)')
    parser.add_argument('--downsample', type=int, default=10,
                        help='Keep 1-in-N NC rays for speed (default: 10)')
    parser.add_argument('--az-offset', type=float, default=None,
                        help='Systematic az offset to subtract from NC-log residual before '
                             'flagging sentinels (default: auto-detect from median)')
    parser.add_argument('--sentinel-threshold', type=float, default=-900.0,
                        help='Az or el values below this are treated as sentinel/fill '
                             '(default: -900)')
    args = parser.parse_args()

    output = args.output or f'axis_vs_nc_{args.datestr}.png'

    log_times, log_az, log_el = load_axis_logs(args.logp_path, args.datestr)
    nc_times, nc_az, nc_el, scan_types, filenames = load_netcdf_rays(
        args.l1_path, args.datestr, downsample=args.downsample)

    # Clip NC times to log coverage (beyond log range, interp would extrapolate)
    in_range = (nc_times >= log_times[0]) & (nc_times <= log_times[-1])
    frac_out = 1.0 - in_range.mean()
    if frac_out > 0.01:
        print(f'Warning: {frac_out*100:.1f}% of NC rays fall outside log time range and will be excluded')
    nc_times   = nc_times[in_range]
    nc_az      = nc_az[in_range]
    nc_el      = nc_el[in_range]
    scan_types = [scan_types[i] for i in np.where(in_range)[0]]
    filenames  = [filenames[i]  for i in np.where(in_range)[0]]

    print('Interpolating log to NC ray times...')
    interp_az, interp_el = interpolate_log(log_times, log_az, log_el, nc_times)

    make_plot(
        args.datestr,
        log_times, log_az, log_el,
        nc_times, nc_az, nc_el,
        interp_az, interp_el,
        scan_types, filenames,
        args.sentinel_threshold,
        args.az_offset,
        output,
    )


if __name__ == '__main__':
    main()
