#!/usr/bin/env python3
"""
fix_midnight_crossing_vpt.py
-----------------------------
General-purpose tool for fixing VPT CF-Radial files whose rays spill past
midnight into the following day.

Operations
----------
1. Trim the source file so it only contains rays on its own calendar day
   (up to 23:59:59 UTC).  A backup (*_original.nc) is made first.

2. Prepend the spillover rays to the next-day file as a new sweep 0.
   If no next-day file exists the spillover is written as a standalone file.

Usage
-----
  python fix_midnight_crossing_vpt.py SOURCE_FILE [NEXT_DAY_FILE]

  SOURCE_FILE    VPT L1 CF-Radial file that crosses midnight.
  NEXT_DAY_FILE  (optional) existing VPT file for the following day.
                 If omitted, the script searches the same directory for a
                 file whose date token is SOURCE_DATE+1 day.

Flags
-----
  --dry-run    Print what would be done without modifying any files.
  --no-backup  Skip creating *_original.nc backups.
  --force      Overwrite existing backup files.
"""

import argparse
import datetime
import os
import re
import shutil
import socket
import sys

import netCDF4 as nc4
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _history_entry(msg: str) -> str:
    now = datetime.datetime.utcnow()
    user = os.environ.get('USER', 'unknown')
    host = socket.gethostname()
    return f"{now.strftime('%a %b %d %H:%M:%S %Y')} - user:{user} machine:{host} {msg}"


def _copy_var_attrs(src, dst):
    for attr in src.ncattrs():
        if attr != '_FillValue':
            setattr(dst, attr, getattr(src, attr))


def _copy_global_attrs(src_ds, dst_ds, overrides=None):
    overrides = overrides or {}
    for attr in src_ds.ncattrs():
        setattr(dst_ds, attr, overrides.get(attr, getattr(src_ds, attr)))
    for k, v in overrides.items():
        if k not in src_ds.ncattrs():
            setattr(dst_ds, k, v)


def _str_to_char_arr(s: str, length: int) -> np.ndarray:
    arr = np.zeros(length, dtype='S1')
    for i, b in enumerate(s.encode('ascii', errors='replace')):
        arr[i] = bytes([b])
    return arr


def _write_char_var(ds, vname: str, s: str, length: int, src_var=None):
    dst = ds.createVariable(vname, 'S1', ('string_length',))
    if src_var is not None:
        _copy_var_attrs(src_var, dst)
    dst[:] = _str_to_char_arr(s, length)
    return dst


def _find_split(time_var) -> int:
    """Return the index of the first ray on the *next* calendar day.

    Raises ValueError if no midnight crossing is found.
    """
    times = nc4.num2date(time_var[:], time_var.units,
                         only_use_cftime_datetimes=False)
    day0 = times[0].day
    for i, t in enumerate(times):
        if t.day != day0:
            return i
    raise ValueError("No midnight crossing detected in this file.")


def _spillover_already_present(next_day_file: str, src_date: datetime.date) -> bool:
    """Return True if next_day_file's first ray already belongs to src_date.

    This indicates the spillover was prepended in a previous run, so we must
    not prepend again when reprocessing src_date.
    """
    with nc4.Dataset(next_day_file, 'r') as ds:
        t = ds.variables['time']
        first = nc4.num2date(t[:1], t.units, only_use_cftime_datetimes=False)[0]
    return (first.year == src_date.year and
            first.month == src_date.month and
            first.day == src_date.day)


def _detect_next_day_file(src_path: str, src_date: datetime.date) -> str | None:
    """Search the parent directory for a VPT file whose date token is src_date+1."""
    next_date = src_date + datetime.timedelta(days=1)
    date_token = next_date.strftime('%Y%m%d')
    parent = os.path.dirname(src_path)
    candidates = []
    for fname in os.listdir(parent):
        if date_token in fname and 'vpt' in fname.lower() and fname.endswith('.nc'):
            candidates.append(os.path.join(parent, fname))
    if not candidates:
        # Also look one directory up / in a sibling date directory
        grandparent = os.path.dirname(parent)
        sibling = os.path.join(grandparent, date_token)
        if os.path.isdir(sibling):
            for fname in os.listdir(sibling):
                if date_token in fname and 'vpt' in fname.lower() and fname.endswith('.nc'):
                    candidates.append(os.path.join(sibling, fname))
    if not candidates:
        return None
    # Prefer the file that sorts first (earliest start time)
    candidates.sort()
    # Exclude any *_original.nc backups
    candidates = [c for c in candidates if '_original' not in c]
    return candidates[0] if candidates else None


def _output_path_for_spillover(src_path: str, first_spill_time: datetime.datetime,
                                next_day_file: str | None) -> str:
    """Generate the output file path for the merged/standalone next-day file."""
    if next_day_file:
        # Replace the time token in the next-day filename
        time_token = first_spill_time.strftime('%Y%m%d-%H%M%S')
        base = os.path.basename(next_day_file)
        # Replace existing YYYYMMDD-HHMMSS or YYYYMMDD token
        new_base = re.sub(r'\d{8}-\d{6}|\d{8}', time_token, base, count=1)
        return os.path.join(os.path.dirname(next_day_file), new_base)
    else:
        # Write standalone spillover file next to the source
        time_token = first_spill_time.strftime('%Y%m%d-%H%M%S')
        base = os.path.basename(src_path)
        new_base = re.sub(r'\d{8}-\d{6}|\d{8}', time_token, base, count=1)
        return os.path.join(os.path.dirname(src_path), new_base)


def _src_date(ncfile: str) -> datetime.date:
    """Extract the calendar date of the source file from its time variable."""
    with nc4.Dataset(ncfile) as ds:
        t = ds.variables['time']
        times = nc4.num2date(t[:1], t.units, only_use_cftime_datetimes=False)
    return datetime.date(times[0].year, times[0].month, times[0].day)


# ---------------------------------------------------------------------------
# Step 1: trim the source file in-place
# ---------------------------------------------------------------------------

def trim_source(src_path: str, split: int, dry_run: bool = False,
                no_backup: bool = False, force: bool = False) -> None:
    """Overwrite src_path so it only contains rays 0 .. split-1."""
    bak = src_path.replace('.nc', '_original.nc')
    if not no_backup:
        if os.path.exists(bak) and not force:
            print(f"  Backup already exists (use --force to overwrite): {bak}")
        elif not dry_run:
            shutil.copy2(src_path, bak)
            print(f"  Backup: {bak}")

    with nc4.Dataset(src_path, 'r') as ds:
        t = ds.variables['time']
        times = nc4.num2date(t[:], t.units, only_use_cftime_datetimes=False)
        tce_str = times[split - 1].strftime('%Y-%m-%dT%H:%M:%SZ')
        sl = ds.dimensions['string_length'].size

        print(f"  Trimming from {ds.dimensions['time'].size} to {split} rays")
        print(f"  New time_coverage_end: {tce_str}")

        if dry_run:
            return

        tmp = src_path + '.trimtmp'
        with nc4.Dataset(tmp, 'w', format='NETCDF4') as out:
            out.createDimension('time', None)
            out.createDimension('range', ds.dimensions['range'].size)
            out.createDimension('sweep', ds.dimensions['sweep'].size)
            out.createDimension('string_length', sl)
            for dim in ('frequency',):
                if dim in ds.dimensions:
                    out.createDimension(dim, ds.dimensions[dim].size)

            _copy_global_attrs(ds, out, overrides={'time_coverage_end': tce_str})

            # Ray-indexed variables
            for vname, src in ds.variables.items():
                if not src.dimensions or src.dimensions[0] != 'time':
                    continue
                fv = src._FillValue if hasattr(src, '_FillValue') else None
                zlib_kw = {} if vname == 'time' else dict(zlib=True, complevel=4)
                dst = out.createVariable(vname, src.dtype, src.dimensions,
                                         fill_value=fv, **zlib_kw)
                _copy_var_attrs(src, dst)
                dst[:] = src[:split] if src.ndim == 1 else src[:split, :]

            # Range
            src = ds.variables['range']
            dst = out.createVariable('range', src.dtype, src.dimensions,
                                     zlib=True, complevel=4)
            _copy_var_attrs(src, dst)
            dst[:] = src[:]

            # Sweep variables — adjust end indices for sweeps that cross midnight
            ss_arr = ds.variables['sweep_start_ray_index'][:]
            se_arr = ds.variables['sweep_end_ray_index'][:]
            # Only keep sweeps whose start is before split
            valid = ss_arr < split
            n_sweeps_new = int(valid.sum())

            # Resize sweep dimension if needed
            if n_sweeps_new < ds.dimensions['sweep'].size:
                # Recreate with smaller sweep dimension
                del out.dimensions  # not deletable; we'll rely on them being unlimited
                # netCDF4 doesn't let us resize fixed dims; recreate the file
                out.close()
                os.unlink(tmp)
                with nc4.Dataset(tmp, 'w', format='NETCDF4') as out2:
                    out2.createDimension('time', None)
                    out2.createDimension('range', ds.dimensions['range'].size)
                    out2.createDimension('sweep', n_sweeps_new)
                    out2.createDimension('string_length', sl)
                    for dim in ('frequency',):
                        if dim in ds.dimensions:
                            out2.createDimension(dim, ds.dimensions[dim].size)

                    _copy_global_attrs(ds, out2, overrides={'time_coverage_end': tce_str})

                    for vname, src in ds.variables.items():
                        if not src.dimensions or src.dimensions[0] != 'time':
                            continue
                        fv = src._FillValue if hasattr(src, '_FillValue') else None
                        zlib_kw = {} if vname == 'time' else dict(zlib=True, complevel=4)
                        dst = out2.createVariable(vname, src.dtype, src.dimensions,
                                                  fill_value=fv, **zlib_kw)
                        _copy_var_attrs(src, dst)
                        dst[:] = src[:split] if src.ndim == 1 else src[:split, :]

                    src = ds.variables['range']
                    dst = out2.createVariable('range', src.dtype, src.dimensions,
                                              zlib=True, complevel=4)
                    _copy_var_attrs(src, dst)
                    dst[:] = src[:]

                    new_ss = ss_arr[valid]
                    new_se = np.minimum(se_arr[valid], split - 1)
                    for vname, data in [
                        ('sweep_start_ray_index', new_ss),
                        ('sweep_end_ray_index',   new_se),
                        ('sweep_number',           np.arange(n_sweeps_new,
                                                  dtype=ds.variables['sweep_number'].dtype)),
                        ('fixed_angle',            ds.variables['fixed_angle'][:][valid]),
                        ('sweep_mode',             ds.variables['sweep_mode'][:][valid]),
                    ]:
                        src = ds.variables[vname]
                        dst = out2.createVariable(vname, src.dtype, src.dimensions,
                                                  zlib=True, complevel=4)
                        _copy_var_attrs(src, dst)
                        dst[:] = data

                    for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                                  'frequency', 'azimuth_correction'):
                        if vname not in ds.variables:
                            continue
                        src = ds.variables[vname]
                        dst = out2.createVariable(vname, src.dtype, src.dimensions)
                        _copy_var_attrs(src, dst)
                        dst[:] = src[:]

                    for vname in ('time_coverage_start', 'time_reference'):
                        if vname not in ds.variables:
                            continue
                        src = ds.variables[vname]
                        dst = out2.createVariable(vname, src.dtype, src.dimensions)
                        _copy_var_attrs(src, dst)
                        dst[:] = src[:]

                    _write_char_var(out2, 'time_coverage_end', tce_str, sl,
                                    src_var=ds.variables.get('time_coverage_end'))

                    existing_hist = getattr(out2, 'history', '')
                    new_entry = _history_entry(
                        f'Trimmed to {split} rays (removed post-midnight spillover '
                        f'to following day); {n_sweeps_new} sweep(s) retained'
                    )
                    out2.history = (new_entry + '\n' + existing_hist).strip()
                    out2.last_revised_date = datetime.datetime.utcnow().strftime(
                        '%Y-%m-%dT%H:%M:%SZ')
                return  # early return; file written above via context manager

            # Simple case: all sweeps start before split, just cap end indices
            new_se = np.minimum(se_arr, split - 1)
            for vname, data in [
                ('sweep_start_ray_index', ss_arr),
                ('sweep_end_ray_index',   new_se),
            ]:
                src = ds.variables[vname]
                dst = out.createVariable(vname, src.dtype, src.dimensions,
                                         zlib=True, complevel=4)
                _copy_var_attrs(src, dst)
                dst[:] = data

            for vname in ('sweep_number', 'fixed_angle', 'sweep_mode'):
                src = ds.variables[vname]
                dst = out.createVariable(vname, src.dtype, src.dimensions,
                                         zlib=True, complevel=4)
                _copy_var_attrs(src, dst)
                dst[:] = src[:]

            for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                          'frequency', 'azimuth_correction'):
                if vname not in ds.variables:
                    continue
                src = ds.variables[vname]
                dst = out.createVariable(vname, src.dtype, src.dimensions)
                _copy_var_attrs(src, dst)
                dst[:] = src[:]

            for vname in ('time_coverage_start', 'time_reference'):
                if vname not in ds.variables:
                    continue
                src = ds.variables[vname]
                dst = out.createVariable(vname, src.dtype, src.dimensions)
                _copy_var_attrs(src, dst)
                dst[:] = src[:]

            _write_char_var(out, 'time_coverage_end', tce_str, sl,
                            src_var=ds.variables.get('time_coverage_end'))

            existing_hist = getattr(out, 'history', '')
            new_entry = _history_entry(
                f'Trimmed to {split} rays (removed post-midnight spillover to following day)'
            )
            out.history = (new_entry + '\n' + existing_hist).strip()
            out.last_revised_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

    os.replace(tmp, src_path)
    print(f"  Written (trimmed): {src_path}")


# ---------------------------------------------------------------------------
# Step 2: prepend spillover to the next-day file
# ---------------------------------------------------------------------------

def prepend_spillover(src_bak: str, split: int, next_day_file: str | None,
                      out_path: str, dry_run: bool = False,
                      no_backup: bool = False, force: bool = False) -> None:
    """Write out_path as: spillover rays from src_bak prepended before next_day_file."""

    with nc4.Dataset(src_bak, 'r') as ds_src:
        t_src = ds_src.variables['time']
        times_src = nc4.num2date(t_src[:], t_src.units, only_use_cftime_datetimes=False)
        n_spill = ds_src.dimensions['time'].size - split

        # Times of spillover rays re-anchored to midnight of the next day
        next_midnight = datetime.datetime(
            times_src[split].year, times_src[split].month, times_src[split].day,
            tzinfo=datetime.timezone.utc
        )
        new_time_units = next_midnight.strftime('seconds since %Y-%m-%dT%H:%M:%SZ')
        midnight_sec = (
            datetime.datetime(times_src[0].year, times_src[0].month, times_src[0].day)
            - datetime.datetime(times_src[split].year, times_src[split].month, times_src[split].day)
        ).total_seconds()  # will be negative (−86400 for a 1-day cross)
        spill_times = t_src[split:] + midnight_sec   # convert to new epoch

        tcs_str = times_src[split].strftime('%Y-%m-%dT%H:%M:%SZ')

        if next_day_file:
            with nc4.Dataset(next_day_file, 'r') as ds_next:
                t_next = ds_next.variables['time']
                times_next = nc4.num2date(t_next[:], t_next.units,
                                          only_use_cftime_datetimes=False)
                n_next = ds_next.dimensions['time'].size
                n_total = n_spill + n_next
                sl = ds_next.dimensions['string_length'].size

                # Existing time values: rebase to new_time_units
                exist_times = nc4.date2num(
                    [datetime.datetime(t.year, t.month, t.day, t.hour, t.minute, t.second)
                     for t in times_next],
                    new_time_units
                )
                tce_str = times_next[-1].strftime('%Y-%m-%dT%H:%M:%SZ')

                # Sweep structure
                ss_src = ds_src.variables['sweep_start_ray_index'][:]
                se_src = ds_src.variables['sweep_end_ray_index'][:]
                fa_src = ds_src.variables['fixed_angle'][:]
                sm_src = ds_src.variables['sweep_mode'][:]
                # Only the sweeps whose rays are entirely in the spillover portion
                # (start >= split); for the last pre-midnight sweep that straddles,
                # treat the whole spillover as one new sweep.
                n_sweeps_next = ds_next.dimensions['sweep'].size
                ss_next = ds_next.variables['sweep_start_ray_index'][:]
                se_next = ds_next.variables['sweep_end_ray_index'][:]
                fa_next = ds_next.variables['fixed_angle'][:]
                sm_next = ds_next.variables['sweep_mode'][:]

                n_sweeps_new = 1 + n_sweeps_next
                new_ss = np.empty(n_sweeps_new, dtype=ss_next.dtype)
                new_se = np.empty(n_sweeps_new, dtype=se_next.dtype)
                new_fa = np.empty(n_sweeps_new, dtype=fa_next.dtype)
                new_sm = np.empty((n_sweeps_new, sl), dtype=sm_next.dtype)

                new_ss[0] = 0
                new_se[0] = n_spill - 1
                new_fa[0] = fa_src[0]
                new_sm[0] = sm_src[0]
                for i in range(n_sweeps_next):
                    new_ss[1 + i] = int(ss_next[i]) + n_spill
                    new_se[1 + i] = int(se_next[i]) + n_spill
                    new_fa[1 + i] = fa_next[i]
                    new_sm[1 + i] = sm_next[i]
                new_sn = np.arange(n_sweeps_new, dtype=ds_next.variables['sweep_number'].dtype)

                print(f"  Spillover rays:      {n_spill}")
                print(f"  Existing next-day rays: {n_next}")
                print(f"  Total merged rays:   {n_total}, {n_sweeps_new} sweep(s)")
                print(f"  time_coverage_start: {tcs_str}")
                print(f"  time_coverage_end:   {tce_str}")
                print(f"  Output: {out_path}")

                if dry_run:
                    return

                if next_day_file != out_path and not no_backup:
                    bak = next_day_file.replace('.nc', '_original.nc')
                    if os.path.exists(bak) and not force:
                        print(f"  Backup already exists (use --force to overwrite): {bak}")
                    else:
                        shutil.copy2(next_day_file, bak)
                        print(f"  Backup: {bak}")

                with nc4.Dataset(out_path, 'w', format='NETCDF4') as out:
                    out.createDimension('time', None)
                    out.createDimension('range', ds_next.dimensions['range'].size)
                    out.createDimension('sweep', n_sweeps_new)
                    out.createDimension('string_length', sl)
                    for dim in ('frequency',):
                        if dim in ds_next.dimensions:
                            out.createDimension(dim, ds_next.dimensions[dim].size)

                    _copy_global_attrs(ds_next, out, overrides={
                        'time_coverage_start': tcs_str,
                        'time_coverage_end':   tce_str,
                    })

                    # Time
                    src_t = ds_next.variables['time']
                    dst_t = out.createVariable('time', src_t.dtype, ('time',))
                    _copy_var_attrs(src_t, dst_t)
                    dst_t.units = new_time_units
                    dst_t[:n_spill] = spill_times
                    dst_t[n_spill:] = exist_times

                    # Other ray-indexed variables
                    ray_vars = [v for v in ds_next.variables
                                if ds_next.variables[v].dimensions
                                and ds_next.variables[v].dimensions[0] == 'time'
                                and v != 'time']
                    for vname in ray_vars:
                        src2 = ds_src.variables.get(vname)
                        src3 = ds_next.variables[vname]
                        fv = src3._FillValue if hasattr(src3, '_FillValue') else None
                        dst = out.createVariable(vname, src3.dtype, src3.dimensions,
                                                 fill_value=fv, zlib=True, complevel=4)
                        _copy_var_attrs(src3, dst)
                        if src2 is not None:
                            if src3.ndim == 1:
                                dst[:n_spill] = src2[split:]
                                dst[n_spill:] = src3[:]
                            else:
                                dst[:n_spill] = src2[split:, :]
                                dst[n_spill:] = src3[:]
                        else:
                            fill = fv if fv is not None else 0
                            if src3.ndim == 1:
                                dst[:n_spill] = fill
                                dst[n_spill:] = src3[:]
                            else:
                                dst[:n_spill] = fill
                                dst[n_spill:] = src3[:]

                    # Range
                    src = ds_next.variables['range']
                    dst = out.createVariable('range', src.dtype, src.dimensions,
                                             zlib=True, complevel=4)
                    _copy_var_attrs(src, dst)
                    dst[:] = src[:]

                    # Sweep variables
                    for vname, data in [
                        ('sweep_start_ray_index', new_ss),
                        ('sweep_end_ray_index',   new_se),
                        ('sweep_number',           new_sn),
                        ('fixed_angle',            new_fa),
                        ('sweep_mode',             new_sm),
                    ]:
                        src = ds_next.variables[vname]
                        dst = out.createVariable(vname, src.dtype, src.dimensions,
                                                 zlib=True, complevel=4)
                        _copy_var_attrs(src, dst)
                        dst[:] = data

                    for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                                  'frequency', 'azimuth_correction'):
                        if vname not in ds_next.variables:
                            continue
                        src = ds_next.variables[vname]
                        dst = out.createVariable(vname, src.dtype, src.dimensions)
                        _copy_var_attrs(src, dst)
                        dst[:] = src[:]

                    for vname in ('time_reference',):
                        if vname not in ds_next.variables:
                            continue
                        src = ds_next.variables[vname]
                        dst = out.createVariable(vname, src.dtype, src.dimensions)
                        _copy_var_attrs(src, dst)
                        dst[:] = src[:]

                    _write_char_var(out, 'time_coverage_start', tcs_str, sl,
                                    src_var=ds_next.variables.get('time_coverage_start'))
                    _write_char_var(out, 'time_coverage_end', tce_str, sl,
                                    src_var=ds_next.variables.get('time_coverage_end'))

                    existing_hist = getattr(out, 'history', '')
                    new_entry = _history_entry(
                        f'Prepended {n_spill} post-midnight spillover ray(s) from '
                        f'{os.path.basename(src_bak)} as new sweep 0; '
                        f'existing {n_sweeps_next} sweep(s) renumbered 1-{n_sweeps_next}'
                    )
                    out.history = (new_entry + '\n' + existing_hist).strip()
                    out.last_revised_date = datetime.datetime.utcnow().strftime(
                        '%Y-%m-%dT%H:%M:%SZ')

                # Remove the stale next-day file now that the correctly-named
                # out_path has been written (Linux allows unlink while open).
                if next_day_file != out_path:
                    os.remove(next_day_file)
                    print(f"  Removed stale file: {os.path.basename(next_day_file)}")

        else:
            # No next-day file — write spillover as standalone
            tce_str = times_src[-1].strftime('%Y-%m-%dT%H:%M:%SZ')
            sl = ds_src.dimensions['string_length'].size
            n_total = n_spill
            n_sweeps_new = 1

            print(f"  Spillover rays:      {n_spill} (no existing next-day file)")
            print(f"  time_coverage_start: {tcs_str}")
            print(f"  time_coverage_end:   {tce_str}")
            print(f"  Output: {out_path}")

            if dry_run:
                return

            with nc4.Dataset(out_path, 'w', format='NETCDF4') as out:
                out.createDimension('time', None)
                out.createDimension('range', ds_src.dimensions['range'].size)
                out.createDimension('sweep', 1)
                out.createDimension('string_length', sl)
                for dim in ('frequency',):
                    if dim in ds_src.dimensions:
                        out.createDimension(dim, ds_src.dimensions[dim].size)

                _copy_global_attrs(ds_src, out, overrides={
                    'time_coverage_start': tcs_str,
                    'time_coverage_end':   tce_str,
                })

                t_src_var = ds_src.variables['time']
                dst_t = out.createVariable('time', t_src_var.dtype, ('time',))
                _copy_var_attrs(t_src_var, dst_t)
                dst_t.units = new_time_units
                dst_t[:] = spill_times

                for vname, src in ds_src.variables.items():
                    if not src.dimensions or src.dimensions[0] != 'time' or vname == 'time':
                        continue
                    fv = src._FillValue if hasattr(src, '_FillValue') else None
                    dst = out.createVariable(vname, src.dtype, src.dimensions,
                                             fill_value=fv, zlib=True, complevel=4)
                    _copy_var_attrs(src, dst)
                    dst[:] = src[split:] if src.ndim == 1 else src[split:, :]

                src = ds_src.variables['range']
                dst = out.createVariable('range', src.dtype, src.dimensions,
                                         zlib=True, complevel=4)
                _copy_var_attrs(src, dst)
                dst[:] = src[:]

                ss_src = ds_src.variables['sweep_start_ray_index']
                se_src = ds_src.variables['sweep_end_ray_index']
                for vname, data in [
                    ('sweep_start_ray_index', np.array([0], dtype=ss_src.dtype)),
                    ('sweep_end_ray_index',   np.array([n_spill - 1], dtype=se_src.dtype)),
                    ('sweep_number',           np.array([0], dtype=ds_src.variables['sweep_number'].dtype)),
                    ('fixed_angle',            ds_src.variables['fixed_angle'][:1]),
                    ('sweep_mode',             ds_src.variables['sweep_mode'][:1]),
                ]:
                    src = ds_src.variables[vname]
                    dst = out.createVariable(vname, src.dtype, src.dimensions,
                                             zlib=True, complevel=4)
                    _copy_var_attrs(src, dst)
                    dst[:] = data

                for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                              'frequency', 'azimuth_correction'):
                    if vname not in ds_src.variables:
                        continue
                    src = ds_src.variables[vname]
                    dst = out.createVariable(vname, src.dtype, src.dimensions)
                    _copy_var_attrs(src, dst)
                    dst[:] = src[:]

                for vname in ('time_reference',):
                    if vname not in ds_src.variables:
                        continue
                    src = ds_src.variables[vname]
                    dst = out.createVariable(vname, src.dtype, src.dimensions)
                    _copy_var_attrs(src, dst)
                    dst[:] = src[:]

                _write_char_var(out, 'time_coverage_start', tcs_str, sl,
                                src_var=ds_src.variables.get('time_coverage_start'))
                _write_char_var(out, 'time_coverage_end', tce_str, sl,
                                src_var=ds_src.variables.get('time_coverage_end'))

                existing_hist = getattr(out, 'history', '')
                new_entry = _history_entry(
                    f'Extracted {n_spill} post-midnight spillover ray(s) from '
                    f'{os.path.basename(src_bak)} as standalone next-day file'
                )
                out.history = (new_entry + '\n' + existing_hist).strip()
                out.last_revised_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

    print(f"  Written: {out_path}")


# ---------------------------------------------------------------------------
# Step 3 (optional): merge a standalone spillover into a subsequently-written
#                    next-day file.
# Used when the standalone was written before the next-day file existed and
# the next-day file has since been created by the batch processor.
# ---------------------------------------------------------------------------

def merge_spillover_into_day(spillover_file: str, main_vpt_file: str,
                              dry_run: bool = False, no_backup: bool = False,
                              force: bool = False) -> None:
    """Prepend all rays from *spillover_file* before the rays in *main_vpt_file*.

    *main_vpt_file* is overwritten in-place (via a tmp file).  A backup of the
    original main file is created as ``*_original.nc`` unless *no_backup* is set.
    The *spillover_file* is removed once it has been successfully integrated.
    """
    with nc4.Dataset(spillover_file, 'r') as ds_spill, \
         nc4.Dataset(main_vpt_file, 'r') as ds_main:

        t_spill = ds_spill.variables['time']
        t_main  = ds_main.variables['time']
        times_spill = nc4.num2date(t_spill[:], t_spill.units,
                                   only_use_cftime_datetimes=False)
        times_main  = nc4.num2date(t_main[:],  t_main.units,
                                   only_use_cftime_datetimes=False)

        n_spill = ds_spill.dimensions['time'].size
        n_main  = ds_main.dimensions['time'].size
        sl      = ds_main.dimensions['string_length'].size

        # Anchor all times to midnight of the spillover day (= start of D+1)
        new_time_units = (
            f"seconds since "
            f"{times_spill[0].year:04d}-{times_spill[0].month:02d}-"
            f"{times_spill[0].day:02d}T00:00:00Z"
        )
        spill_seconds = nc4.date2num(
            [datetime.datetime(t.year, t.month, t.day, t.hour, t.minute, t.second)
             for t in times_spill],
            new_time_units
        )
        main_seconds = nc4.date2num(
            [datetime.datetime(t.year, t.month, t.day, t.hour, t.minute, t.second)
             for t in times_main],
            new_time_units
        )

        tcs_str = times_spill[0].strftime('%Y-%m-%dT%H:%M:%SZ')
        tce_str = times_main[-1].strftime('%Y-%m-%dT%H:%M:%SZ')

        # Build merged sweep structure
        n_sweeps_spill = ds_spill.dimensions['sweep'].size
        n_sweeps_main  = ds_main.dimensions['sweep'].size
        n_sweeps_new   = n_sweeps_spill + n_sweeps_main

        ss_spill = ds_spill.variables['sweep_start_ray_index'][:]
        se_spill = ds_spill.variables['sweep_end_ray_index'][:]
        fa_spill = ds_spill.variables['fixed_angle'][:]
        sm_spill = ds_spill.variables['sweep_mode'][:]

        ss_main = ds_main.variables['sweep_start_ray_index'][:]
        se_main = ds_main.variables['sweep_end_ray_index'][:]
        fa_main = ds_main.variables['fixed_angle'][:]
        sm_main = ds_main.variables['sweep_mode'][:]

        new_ss = np.concatenate([ss_spill, ss_main + n_spill])
        new_se = np.concatenate([se_spill, se_main + n_spill])
        new_fa = np.concatenate([fa_spill, fa_main])
        new_sm = np.concatenate([sm_spill, sm_main], axis=0)
        new_sn = np.arange(n_sweeps_new, dtype=ss_main.dtype)

        print(f"  Spillover rays: {n_spill}, main rays: {n_main}, total: {n_spill + n_main}")
        print(f"  time_coverage_start: {tcs_str}")
        print(f"  time_coverage_end:   {tce_str}")
        print(f"  Output: {main_vpt_file}")

        if dry_run:
            return

        if not no_backup:
            bak = main_vpt_file.replace('.nc', '_original.nc')
            if os.path.exists(bak) and not force:
                print(f"  Backup already exists (use --force to overwrite): {bak}")
            else:
                shutil.copy2(main_vpt_file, bak)
                print(f"  Backup: {bak}")

        tmp = main_vpt_file + '.mergetmp'
        with nc4.Dataset(tmp, 'w', format='NETCDF4') as out:
            out.createDimension('time', None)
            out.createDimension('range', ds_main.dimensions['range'].size)
            out.createDimension('sweep', n_sweeps_new)
            out.createDimension('string_length', sl)
            for dim in ('frequency',):
                if dim in ds_main.dimensions:
                    out.createDimension(dim, ds_main.dimensions[dim].size)

            _copy_global_attrs(ds_main, out, overrides={
                'time_coverage_start': tcs_str,
                'time_coverage_end':   tce_str,
            })

            # Time
            src_t = ds_main.variables['time']
            dst_t = out.createVariable('time', src_t.dtype, ('time',))
            _copy_var_attrs(src_t, dst_t)
            dst_t.units = new_time_units
            dst_t[:n_spill] = spill_seconds
            dst_t[n_spill:] = main_seconds

            # Ray-indexed variables
            ray_vars = [v for v in ds_main.variables
                        if ds_main.variables[v].dimensions
                        and ds_main.variables[v].dimensions[0] == 'time'
                        and v != 'time']
            for vname in ray_vars:
                src_m = ds_main.variables[vname]
                src_s = ds_spill.variables.get(vname)
                fv = src_m._FillValue if hasattr(src_m, '_FillValue') else None
                dst = out.createVariable(vname, src_m.dtype, src_m.dimensions,
                                         fill_value=fv, zlib=True, complevel=4)
                _copy_var_attrs(src_m, dst)
                fill = fv if fv is not None else 0
                if src_m.ndim == 1:
                    dst[:n_spill] = src_s[:] if src_s is not None else fill
                    dst[n_spill:] = src_m[:]
                else:
                    dst[:n_spill] = src_s[:, :] if src_s is not None else fill
                    dst[n_spill:] = src_m[:, :]

            # Range
            src = ds_main.variables['range']
            dst = out.createVariable('range', src.dtype, src.dimensions,
                                     zlib=True, complevel=4)
            _copy_var_attrs(src, dst)
            dst[:] = src[:]

            # Sweep variables
            for vname, data in [
                ('sweep_start_ray_index', new_ss),
                ('sweep_end_ray_index',   new_se),
                ('sweep_number',          new_sn),
                ('fixed_angle',           new_fa),
                ('sweep_mode',            new_sm),
            ]:
                src = ds_main.variables[vname]
                dst = out.createVariable(vname, src.dtype, src.dimensions,
                                         zlib=True, complevel=4)
                _copy_var_attrs(src, dst)
                dst[:] = data

            for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                          'frequency', 'azimuth_correction'):
                if vname not in ds_main.variables:
                    continue
                src = ds_main.variables[vname]
                dst = out.createVariable(vname, src.dtype, src.dimensions)
                _copy_var_attrs(src, dst)
                dst[:] = src[:]

            for vname in ('time_reference',):
                if vname not in ds_main.variables:
                    continue
                src = ds_main.variables[vname]
                dst = out.createVariable(vname, src.dtype, src.dimensions)
                _copy_var_attrs(src, dst)
                dst[:] = src[:]

            _write_char_var(out, 'time_coverage_start', tcs_str, sl,
                            src_var=ds_main.variables.get('time_coverage_start'))
            _write_char_var(out, 'time_coverage_end', tce_str, sl,
                            src_var=ds_main.variables.get('time_coverage_end'))

            existing_hist = getattr(out, 'history', '')
            new_entry = _history_entry(
                f'Prepended {n_spill} post-midnight spillover ray(s) from standalone '
                f'{os.path.basename(spillover_file)}; existing {n_sweeps_main} sweep(s) '
                f'renumbered {n_sweeps_spill}-{n_sweeps_new - 1}'
            )
            out.history = (new_entry + '\n' + existing_hist).strip()
            out.last_revised_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

    os.replace(tmp, main_vpt_file)
    print(f"  Written (merged): {main_vpt_file}")
    os.unlink(spillover_file)
    print(f"  Removed standalone spillover: {os.path.basename(spillover_file)}")

    # Rename the merged file so its timestamp token reflects the new (earlier) start time.
    correct_time_token = times_spill[0].strftime('%Y%m%d-%H%M%S')
    correct_base = re.sub(r'\d{8}-\d{6}|\d{8}', correct_time_token,
                          os.path.basename(main_vpt_file), count=1)
    correct_path = os.path.join(os.path.dirname(main_vpt_file), correct_base)
    if correct_path != main_vpt_file:
        os.rename(main_vpt_file, correct_path)
        print(f"  Renamed to: {correct_base}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Fix a VPT CF-Radial file whose rays spill past midnight.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('source', help='Source VPT .nc file that crosses midnight')
    parser.add_argument('next_day', nargs='?', default=None,
                        help='Existing VPT .nc file for the following day (optional)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print what would be done without modifying files')
    parser.add_argument('--no-backup', action='store_true',
                        help='Skip creating *_original.nc backups')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing backup files')
    args = parser.parse_args()

    src = os.path.abspath(args.source)
    if not os.path.isfile(src):
        sys.exit(f'Source file not found: {src}')

    # Find split point
    print(f'Source: {src}')
    with nc4.Dataset(src) as ds:
        try:
            split = _find_split(ds.variables['time'])
        except ValueError as e:
            sys.exit(str(e))
        t = ds.variables['time']
        times = nc4.num2date(t[:], t.units, only_use_cftime_datetimes=False)

    print(f'  Midnight crossing at ray index {split}')
    print(f'    Last ray before midnight: [{split-1}] {times[split-1]}')
    print(f'    First spillover ray:      [{split}]   {times[split]}')
    print(f'  Total rays: {len(times)}, spillover: {len(times) - split}')

    # Locate or validate next-day file
    src_date = datetime.date(times[0].year, times[0].month, times[0].day)
    if args.next_day:
        next_day_file = os.path.abspath(args.next_day)
        if not os.path.isfile(next_day_file):
            sys.exit(f'Next-day file not found: {next_day_file}')
    else:
        next_day_file = _detect_next_day_file(src, src_date)
        if next_day_file:
            print(f'  Auto-detected next-day file: {next_day_file}')
        else:
            print('  No next-day file found — spillover will be written standalone.')

    first_spill = datetime.datetime(times[split].year, times[split].month,
                                    times[split].day, times[split].hour,
                                    times[split].minute, times[split].second)
    out_path = _output_path_for_spillover(src, first_spill, next_day_file)

    print()
    print('=== Step 1: Trim source file ===')
    trim_source(src, split, dry_run=args.dry_run, no_backup=args.no_backup,
                force=args.force)

    print()
    print('=== Step 2: Write next-day file with prepended spillover ===')
    # Use the backup as the source for reading spillover (trim_source may have
    # overwritten the original).
    src_bak = src.replace('.nc', '_original.nc') if not args.no_backup else src
    if args.no_backup:
        # Re-read split from original (src still has trimmed data if dry_run=False)
        # We need the original — can't proceed without backup in non-dry-run mode.
        if not args.dry_run:
            sys.exit('--no-backup is incompatible with a non-dry-run after the source '
                     'has been trimmed. Re-run with backup enabled.')
        src_bak = src
    prepend_spillover(src_bak, split, next_day_file, out_path,
                      dry_run=args.dry_run, no_backup=args.no_backup, force=args.force)

    print()
    if args.dry_run:
        print('Dry run complete — no files modified.')
    else:
        print('Done.')


if __name__ == '__main__':
    main()
