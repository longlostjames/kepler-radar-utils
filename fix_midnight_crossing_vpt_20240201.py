"""
fix_midnight_crossing_vpt_20240201.py
---------------------------------------
One-off fix for ncas-mobile-ka-band-radar-1_cao_20240201-171258_vpt_l1_v1.0.1.nc
which runs past midnight into 2024-02-02.

Operations
----------
1. Trim the 20240201 file to rays 0-9535 (17:12:58 - 23:59:59Z on 2024-02-01).
   The original file is backed up as *_original.nc before overwriting.
   time_coverage_end, history, and last_revised_date are updated.

2. Create a new standalone 20240202-000000_vpt file from the spillover rays
   9536-26602 (2024-02-02T00:00:00Z - 2024-02-02T12:08:13Z).
   This new single-sweep file sits alongside the already-existing
   20240202-143823_vpt file in the 20240202 L1 directory (there is a ~2.5 h
   gap between the two files for which no raw data are available).

File structure
--------------
20240201-171258_vpt  (original, 26603 rays, 1 sweep)
  time units: seconds since 2024-02-01T00:00:00Z
  ray 0:     17:12:58 UTC 2024-02-01
  ray 9535:  23:59:59 UTC 2024-02-01   <- trim here
  ray 9536:  00:00:00 UTC 2024-02-02   <- first spillover ray

New 20240202-000000_vpt  (17067 rays, 1 sweep)
  time units: seconds since 2024-02-02T00:00:00Z
  00:00:00Z - 12:08:13Z on 2024-02-02
"""

import datetime
import os
import shutil
import socket

import netCDF4 as nc4
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = (
    '/gws/pw/j07/ncas_obs_vol2/cao/processing/'
    'ncas-mobile-ka-band-radar-1/ccrest-m/L1_v1.0.1'
)

SRC01 = os.path.join(
    BASE, '20240201',
    'ncas-mobile-ka-band-radar-1_cao_20240201-171258_vpt_l1_v1.0.1.nc'
)
OUT02 = os.path.join(
    BASE, '20240202',
    'ncas-mobile-ka-band-radar-1_cao_20240202-000000_vpt_l1_v1.0.1.nc'
)

# Index of the first ray on 2024-02-02 inside the 20240201 file
SPLIT = 9536


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _history_entry(msg):
    now = datetime.datetime.utcnow()
    user = os.environ.get('USER', 'unknown')
    host = socket.gethostname()
    timestamp = now.strftime('%a %b %d %H:%M:%S %Y')
    return f'{timestamp} - user:{user} machine:{host} {msg}'


def copy_var_attrs(src_var, dst_var):
    for attr in src_var.ncattrs():
        if attr != '_FillValue':
            setattr(dst_var, attr, getattr(src_var, attr))


def copy_global_attrs(src_ds, dst_ds, overrides=None):
    overrides = overrides or {}
    for attr in src_ds.ncattrs():
        val = overrides.get(attr, getattr(src_ds, attr))
        setattr(dst_ds, attr, val)
    for k, v in overrides.items():
        if k not in src_ds.ncattrs():
            setattr(dst_ds, k, v)


def str_to_char_arr(s, length):
    arr = np.zeros(length, dtype='S1')
    for i, b in enumerate(s.encode('ascii')):
        arr[i] = bytes([b])
    return arr


def write_char_var(out_ds, vname, s, length, src_var=None):
    dst = out_ds.createVariable(vname, 'S1', ('string_length',))
    if src_var is not None:
        copy_var_attrs(src_var, dst)
    dst[:] = str_to_char_arr(s, length)
    return dst


# ---------------------------------------------------------------------------
# 1.  Trim the 20240201 file (rays 0 .. SPLIT-1)
# ---------------------------------------------------------------------------
print("=== Step 1: Trimming 20240201 file ===")

bak01 = SRC01.replace('.nc', '_original.nc')
if not os.path.exists(bak01):
    shutil.copy2(SRC01, bak01)
    print(f"  Backup: {bak01}")
else:
    print(f"  Backup already exists, skipping copy: {bak01}")

ds01 = nc4.Dataset(SRC01, 'r')
t01 = ds01.variables['time']
times01 = nc4.num2date(t01[:], t01.units, only_use_cftime_datetimes=False)

# Sanity checks
assert times01[SPLIT - 1].day == 1, f"Ray {SPLIT-1} should be on 2024-02-01"
assert times01[SPLIT].day == 2, f"Ray {SPLIT} should be on 2024-02-02"
print(f"  Last ray before midnight: [{SPLIT-1}] {times01[SPLIT-1]}")
print(f"  First ray after midnight: [{SPLIT}]  {times01[SPLIT]}")

sl = ds01.dimensions['string_length'].size
tce01_str = times01[SPLIT - 1].strftime('%Y-%m-%dT%H:%M:%SZ')
print(f"  Updated time_coverage_end: {tce01_str}")
print(f"  Keeping {SPLIT} rays in 20240201 file")

tmp01 = SRC01 + '.tmp'
with nc4.Dataset(tmp01, 'w', format='NETCDF4') as out:
    out.createDimension('time', None)
    out.createDimension('range', ds01.dimensions['range'].size)
    out.createDimension('sweep', 1)
    out.createDimension('string_length', sl)
    out.createDimension('frequency', ds01.dimensions['frequency'].size)

    copy_global_attrs(ds01, out, overrides={'time_coverage_end': tce01_str})

    # Ray-indexed variables
    ray_vars = [
        v for v in ds01.variables
        if ds01.variables[v].dimensions
        and ds01.variables[v].dimensions[0] == 'time'
    ]
    for vname in ray_vars:
        src = ds01.variables[vname]
        fv = src._FillValue if hasattr(src, '_FillValue') else None
        zlib_kw = {} if vname == 'time' else dict(zlib=True, complevel=4)
        dst = out.createVariable(vname, src.dtype, src.dimensions,
                                 fill_value=fv, **zlib_kw)
        copy_var_attrs(src, dst)
        if src.ndim == 1:
            dst[:] = src[:SPLIT]
        else:
            dst[:] = src[:SPLIT, :]

    # Range
    src = ds01.variables['range']
    dst = out.createVariable('range', src.dtype, src.dimensions, zlib=True, complevel=4)
    copy_var_attrs(src, dst)
    dst[:] = src[:]

    # Sweep variables: 1 sweep, end ray trimmed
    for vname, data in [
        ('sweep_start_ray_index', np.array([0], dtype=ds01.variables['sweep_start_ray_index'].dtype)),
        ('sweep_end_ray_index',   np.array([SPLIT - 1], dtype=ds01.variables['sweep_end_ray_index'].dtype)),
    ]:
        src = ds01.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = data

    for vname in ('sweep_number', 'fixed_angle', 'sweep_mode'):
        src = ds01.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # Scalar / frequency variables
    for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                  'frequency', 'azimuth_correction'):
        if vname not in ds01.variables:
            continue
        src = ds01.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # String char variables
    for vname in ('time_coverage_start', 'time_reference'):
        if vname not in ds01.variables:
            continue
        src = ds01.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    write_char_var(out, 'time_coverage_end', tce01_str, sl,
                   src_var=ds01.variables.get('time_coverage_end'))

    out.featureType = 'timeSeriesProfile'

    existing_hist = getattr(out, 'history', '')
    new_entry = _history_entry(
        f'Trimmed to {SPLIT} rays (removed post-midnight data; '
        f'continuation in ncas-mobile-ka-band-radar-1_cao_20240202-000000_vpt_l1_v1.0.1.nc)'
    )
    out.history = (new_entry + '\n' + existing_hist).strip()
    out.last_revised_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

ds01.close()
os.replace(tmp01, SRC01)
print(f"  Written: {SRC01}")


# ---------------------------------------------------------------------------
# 2.  Create new 20240202-000000_vpt from spillover rays
# ---------------------------------------------------------------------------
print("\n=== Step 2: Creating 20240202-000000_vpt file ===")

# Ensure output directory exists
os.makedirs(os.path.dirname(OUT02), exist_ok=True)

if os.path.exists(OUT02):
    print(f"  WARNING: {OUT02} already exists — overwriting")

ds01 = nc4.Dataset(bak01, 'r')   # read original (pre-trim) file
t01 = ds01.variables['time']
times01 = nc4.num2date(t01[:], t01.units, only_use_cftime_datetimes=False)

n_spill = int(ds01.dimensions['time'].size) - SPLIT   # 17067
sl = ds01.dimensions['string_length'].size

# Time values: convert from "seconds since 2024-02-01T00:00:00Z"
# to "seconds since 2024-02-02T00:00:00Z" by subtracting 86400.
spill_times = t01[SPLIT:] - 86400.0
new_time_units = 'seconds since 2024-02-02T00:00:00Z'

tcs02_str = times01[SPLIT].strftime('%Y-%m-%dT%H:%M:%SZ')      # 2024-02-02T00:00:00Z
tce02_str = times01[-1].strftime('%Y-%m-%dT%H:%M:%SZ')         # 2024-02-02T12:08:13Z

print(f"  Spillover rays: {n_spill} (indices {SPLIT}-{ds01.dimensions['time'].size - 1})")
print(f"  time_coverage_start: {tcs02_str}")
print(f"  time_coverage_end:   {tce02_str}")

with nc4.Dataset(OUT02, 'w', format='NETCDF4') as out:
    out.createDimension('time', None)
    out.createDimension('range', ds01.dimensions['range'].size)
    out.createDimension('sweep', 1)
    out.createDimension('string_length', sl)
    out.createDimension('frequency', ds01.dimensions['frequency'].size)

    copy_global_attrs(ds01, out, overrides={
        'time_coverage_start': tcs02_str,
        'time_coverage_end':   tce02_str,
    })

    # Time variable with new units and spillover values
    src_t = ds01.variables['time']
    dst_t = out.createVariable('time', src_t.dtype, ('time',))
    copy_var_attrs(src_t, dst_t)
    dst_t.units = new_time_units
    dst_t[:] = spill_times

    # Other ray-indexed variables
    ray_vars = [
        v for v in ds01.variables
        if ds01.variables[v].dimensions
        and ds01.variables[v].dimensions[0] == 'time'
        and v != 'time'
    ]
    for vname in ray_vars:
        src = ds01.variables[vname]
        fv = src._FillValue if hasattr(src, '_FillValue') else None
        dst = out.createVariable(vname, src.dtype, src.dimensions,
                                 fill_value=fv, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        if src.ndim == 1:
            dst[:] = src[SPLIT:]
        else:
            dst[:] = src[SPLIT:, :]

    # Range
    src = ds01.variables['range']
    dst = out.createVariable('range', src.dtype, src.dimensions, zlib=True, complevel=4)
    copy_var_attrs(src, dst)
    dst[:] = src[:]

    # Sweep variables: single sweep covering all spillover rays
    for vname, data in [
        ('sweep_start_ray_index', np.array([0], dtype=ds01.variables['sweep_start_ray_index'].dtype)),
        ('sweep_end_ray_index',   np.array([n_spill - 1], dtype=ds01.variables['sweep_end_ray_index'].dtype)),
        ('sweep_number',          np.array([0], dtype=ds01.variables['sweep_number'].dtype)),
    ]:
        src = ds01.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = data

    for vname in ('fixed_angle', 'sweep_mode'):
        src = ds01.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = src[:]   # single sweep in original, keep as-is

    # Scalar / frequency variables
    for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                  'frequency', 'azimuth_correction'):
        if vname not in ds01.variables:
            continue
        src = ds01.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # String char variables
    write_char_var(out, 'time_coverage_start', tcs02_str, sl,
                   src_var=ds01.variables.get('time_coverage_start'))
    write_char_var(out, 'time_coverage_end', tce02_str, sl,
                   src_var=ds01.variables.get('time_coverage_end'))
    for vname in ('time_reference',):
        if vname not in ds01.variables:
            continue
        src = ds01.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    out.featureType = 'timeSeriesProfile'

    existing_hist = getattr(out, 'history', '')
    new_entry = _history_entry(
        f'Created from spillover rays {SPLIT}-{ds01.dimensions["time"].size - 1} of '
        f'ncas-mobile-ka-band-radar-1_cao_20240201-171258_vpt_l1_v1.0.1.nc '
        f'(post-midnight data split at day boundary)'
    )
    out.history = (new_entry + '\n' + existing_hist).strip()
    out.last_revised_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

ds01.close()
print(f"  Written: {OUT02}")
print("\nDone.")
print(f"  20240201 trimmed file:  {SRC01}")
print(f"  New 20240202 early VPT: {OUT02}")
print(f"  Note: ~2.5 h gap between {tce02_str} and the existing 20240202-143823_vpt file is expected (no raw data available).")
