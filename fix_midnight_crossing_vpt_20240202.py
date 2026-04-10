"""
fix_midnight_crossing_vpt_20240202.py
---------------------------------------
One-off fix for ncas-mobile-ka-band-radar-1_cao_20240202-143823_vpt_l1_v1.0.1.nc
which runs past midnight into 2024-02-03.

NOTE on "earlier data"
----------------------
The only raw mmclx file for 2024-02-02 is 20240202_143853.mmclx.gz, whose first
ray is at 14:38:23 UTC.  No raw files from earlier in the day are available in
the CCREST-M mom directory.  If earlier VPT data exist in another format (e.g.
MIRA native .nc.gz from a different recording session), they would need to be
ingested separately before being prepended to the L1 file.  This script handles
only the midnight-crossing problem.

Operations
----------
1. Trim the 20240202 file to rays 0-6566 (up to 2024-02-02T23:59:59Z).
   The single sweep is truncated at its new end ray.
   The original file is backed up as *_original.nc before overwriting.
   time_coverage_end and history/last_revised_date are updated.

2. Create a merged 20240203 VPT file by prepending the spillover rays
   6567-13282 (2024-02-03T00:00:03Z - 2024-02-03T09:34:03Z) as a new sweep 0
   to the existing 20240203-123024 file (which covers 12:30:24 - 15:08:47 UTC).
   The merged file is written as …20240203-000003_vpt_l1_v1.0.1.nc.
   The old 123024 file is backed up as *_original.nc.
   time_coverage_start, history, and last_revised_date are updated.

File structure
--------------
20240202-143823_vpt  (original, 13283 rays, 1 sweep)
  time units: seconds since 2024-02-02T00:00:00Z
  ray 0:    14:38:23 UTC  (t=52703.66 s)
  ray 6566: 23:59:59 UTC  (t=86399.x  s)  <- trim here
  ray 6567: 00:00:03 UTC 20240203           <- first spillover ray

20240203-123024_vpt  (original, 824 rays, 8 sweeps)
  time units: seconds since 2024-02-03T00:00:00Z
  ray 0:  12:30:24 UTC
  ray 823: 15:08:47 UTC

Merged 20240203-000003_vpt  (6716 + 824 = 7540 rays, 9 sweeps)
  sweep 0: spillover (rays 0-6715), 00:00:03 - 09:34:03 UTC
  sweeps 1-8: existing 20240203 sweeps, ray indices offset by 6716
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

SRC02 = os.path.join(
    BASE, '20240202',
    'ncas-mobile-ka-band-radar-1_cao_20240202-143823_vpt_l1_v1.0.1.nc'
)
SRC03 = os.path.join(
    BASE, '20240203',
    'ncas-mobile-ka-band-radar-1_cao_20240203-123024_vpt_l1_v1.0.1.nc'
)
OUT03 = os.path.join(
    BASE, '20240203',
    'ncas-mobile-ka-band-radar-1_cao_20240203-000003_vpt_l1_v1.0.1.nc'
)

# Index of the first ray on 2024-02-03 inside the 20240202 file
SPLIT = 6567


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
    """Create/write a string_length character variable."""
    dst = out_ds.createVariable(vname, 'S1', ('string_length',))
    if src_var is not None:
        copy_var_attrs(src_var, dst)
    dst[:] = str_to_char_arr(s, length)
    return dst


# ---------------------------------------------------------------------------
# 1.  Trim the 20240202 file (rays 0 .. SPLIT-1)
# ---------------------------------------------------------------------------
print("=== Step 1: Trimming 20240202 file ===")

bak02 = SRC02.replace('.nc', '_original.nc')
if not os.path.exists(bak02):
    shutil.copy2(SRC02, bak02)
    print(f"  Backup: {bak02}")
else:
    print(f"  Backup already exists, skipping copy: {bak02}")

ds02 = nc4.Dataset(SRC02, 'r')
t02 = ds02.variables['time']
times02 = nc4.num2date(t02[:], t02.units, only_use_cftime_datetimes=False)

# Verify SPLIT makes sense
assert times02[SPLIT - 1].day == 2, f"Ray {SPLIT-1} should be on 2024-02-02"
assert times02[SPLIT].day == 3, f"Ray {SPLIT} should be on 2024-02-03"
print(f"  Last ray before midnight: [{SPLIT-1}] {times02[SPLIT-1]}")
print(f"  First ray after midnight: [{SPLIT}]  {times02[SPLIT]}")

sl = ds02.dimensions['string_length'].size
tce02_str = times02[SPLIT - 1].strftime('%Y-%m-%dT%H:%M:%SZ')

n_rays_trimmed = SPLIT  # 6567

print(f"  Keeping {n_rays_trimmed} rays in 20240202 file")
print(f"  Updated time_coverage_end: {tce02_str}")

tmp02 = SRC02 + '.tmp'
with nc4.Dataset(tmp02, 'w', format='NETCDF4') as out:
    # Dimensions: 1 sweep, SPLIT rays
    out.createDimension('time', None)
    out.createDimension('range', ds02.dimensions['range'].size)
    out.createDimension('sweep', 1)
    out.createDimension('string_length', sl)
    out.createDimension('frequency', ds02.dimensions['frequency'].size)

    copy_global_attrs(ds02, out, overrides={'time_coverage_end': tce02_str})

    # Ray-indexed variables
    ray_vars = [
        v for v in ds02.variables
        if ds02.variables[v].dimensions
        and ds02.variables[v].dimensions[0] == 'time'
    ]
    for vname in ray_vars:
        src = ds02.variables[vname]
        fv = src._FillValue if hasattr(src, '_FillValue') else None
        zlib_kw = {} if vname == 'time' else dict(zlib=True, complevel=4)
        dst = out.createVariable(vname, src.dtype, src.dimensions,
                                 fill_value=fv, **zlib_kw)
        copy_var_attrs(src, dst)
        if src.ndim == 1:
            dst[:] = src[:SPLIT]
        else:
            dst[:] = src[:SPLIT, :]

    # Range variable
    src = ds02.variables['range']
    dst = out.createVariable('range', src.dtype, src.dimensions, zlib=True, complevel=4)
    copy_var_attrs(src, dst)
    dst[:] = src[:]

    # Sweep variables: 1 sweep, trimmed end ray
    src_ss = ds02.variables['sweep_start_ray_index']
    src_se = ds02.variables['sweep_end_ray_index']
    for vname, data in [
        ('sweep_start_ray_index', np.array([0], dtype=src_ss.dtype)),
        ('sweep_end_ray_index',   np.array([SPLIT - 1], dtype=src_se.dtype)),
    ]:
        src = ds02.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = data

    for vname in ('sweep_number', 'fixed_angle', 'sweep_mode'):
        src = ds02.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = src[:]   # only 1 sweep in original, keep as-is

    # Scalar / frequency variables
    for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                  'frequency', 'azimuth_correction'):
        if vname not in ds02.variables:
            continue
        src = ds02.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # String character variables
    for vname in ('time_coverage_start', 'time_reference'):
        if vname not in ds02.variables:
            continue
        src = ds02.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # Override time_coverage_end
    write_char_var(out, 'time_coverage_end', tce02_str, sl,
                   src_var=ds02.variables.get('time_coverage_end'))

    # featureType
    out.featureType = 'timeSeriesProfile'

    # History
    existing_hist = getattr(out, 'history', '')
    new_entry = _history_entry(
        f'Trimmed to {n_rays_trimmed} rays (removed post-midnight data; '
        f'continuation in ncas-mobile-ka-band-radar-1_cao_20240203-000003_vpt_l1_v1.0.1.nc)'
    )
    out.history = (new_entry + '\n' + existing_hist).strip()
    out.last_revised_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

ds02.close()
os.replace(tmp02, SRC02)
print(f"  Written: {SRC02}")


# ---------------------------------------------------------------------------
# 2.  Create merged 20240203 file
# ---------------------------------------------------------------------------
print("\n=== Step 2: Creating merged 20240203 file ===")

bak03 = SRC03.replace('.nc', '_original.nc')
if not os.path.exists(bak03):
    shutil.copy2(SRC03, bak03)
    print(f"  Backup: {bak03}")
else:
    print(f"  Backup already exists, skipping copy: {bak03}")

# Re-open original (pre-trim) 20240202 file to read spillover
ds02 = nc4.Dataset(bak02, 'r')
ds03 = nc4.Dataset(SRC03, 'r')

t02 = ds02.variables['time']
t03 = ds03.variables['time']

times02 = nc4.num2date(t02[:], t02.units, only_use_cftime_datetimes=False)
times03 = nc4.num2date(t03[:], t03.units, only_use_cftime_datetimes=False)

n_spill = int(ds02.dimensions['time'].size) - SPLIT   # 6716
n_03    = int(ds03.dimensions['time'].size)            # 824
n_total = n_spill + n_03                               # 7540

# Spillover time values: convert from "seconds since 2024-02-02T00:00:00Z"
# to "seconds since 2024-02-03T00:00:00Z" by subtracting 86400.
spill_times = t02[SPLIT:] - 86400.0
exist_times = t03[:]

new_time_units = 'seconds since 2024-02-03T00:00:00Z'

# Sweep structure
ss02 = ds02.variables['sweep_start_ray_index'][:]
se02 = ds02.variables['sweep_end_ray_index'][:]
fa02 = ds02.variables['fixed_angle'][:]
sm02 = ds02.variables['sweep_mode'][:]    # (1, string_length)
sn02 = ds02.variables['sweep_number'][:]

ss03 = ds03.variables['sweep_start_ray_index'][:]
se03 = ds03.variables['sweep_end_ray_index'][:]
fa03 = ds03.variables['fixed_angle'][:]
sm03 = ds03.variables['sweep_mode'][:]    # (8, string_length)
sn03 = ds03.variables['sweep_number'][:]

n_sweeps03 = int(ds03.dimensions['sweep'].size)   # 8
n_new_sweeps = 1 + n_sweeps03                      # 9
sl = ds03.dimensions['string_length'].size

new_ss = np.zeros(n_new_sweeps, dtype=ss03.dtype)
new_se = np.zeros(n_new_sweeps, dtype=se03.dtype)
new_fa = np.zeros(n_new_sweeps, dtype=fa03.dtype)
new_sm = np.zeros((n_new_sweeps, sl), dtype=sm03.dtype)
new_sn = np.arange(n_new_sweeps, dtype=sn03.dtype)

# Sweep 0: spillover from 20240202 (same fixed_angle/mode as original single sweep)
new_ss[0] = 0
new_se[0] = n_spill - 1
new_fa[0] = fa02[0]
new_sm[0, :] = sm02[0, :]

# Sweeps 1..8: existing 20240203 sweeps, offset by n_spill
for i in range(n_sweeps03):
    new_ss[1 + i] = int(ss03[i]) + n_spill
    new_se[1 + i] = int(se03[i]) + n_spill
    new_fa[1 + i] = fa03[i]
    new_sm[1 + i, :] = sm03[i, :]

# Time coverage strings
tcs03_str = times02[SPLIT].strftime('%Y-%m-%dT%H:%M:%SZ')    # 2024-02-03T00:00:03Z
tce03_str = times03[-1].strftime('%Y-%m-%dT%H:%M:%SZ')       # 2024-02-03T15:08:47Z

print(f"  Spillover rays: {n_spill} (indices {SPLIT}-{ds02.dimensions['time'].size - 1} in original 20240202 file)")
print(f"  Existing 20240203 rays: {n_03}")
print(f"  Merged total: {n_total} rays, {n_new_sweeps} sweeps")
print(f"  time_coverage_start: {tcs03_str}")
print(f"  time_coverage_end:   {tce03_str}")

with nc4.Dataset(OUT03, 'w', format='NETCDF4') as out:
    out.createDimension('time', None)
    out.createDimension('range', ds03.dimensions['range'].size)
    out.createDimension('sweep', n_new_sweeps)
    out.createDimension('string_length', sl)
    out.createDimension('frequency', ds03.dimensions['frequency'].size)

    copy_global_attrs(ds03, out, overrides={
        'time_coverage_start': tcs03_str,
        'time_coverage_end':   tce03_str,
    })

    # --- Time variable: concatenate spillover + existing, new units ---
    src_t02 = ds02.variables['time']
    src_t03 = ds03.variables['time']
    dst_t = out.createVariable('time', src_t03.dtype, ('time',))
    copy_var_attrs(src_t03, dst_t)
    dst_t.units = new_time_units
    dst_t[:n_spill] = spill_times
    dst_t[n_spill:] = exist_times

    # --- Other ray-indexed variables ---
    ray_vars = [
        v for v in ds03.variables
        if ds03.variables[v].dimensions
        and ds03.variables[v].dimensions[0] == 'time'
        and v != 'time'
    ]
    for vname in ray_vars:
        src2 = ds02.variables.get(vname)
        src3 = ds03.variables[vname]
        fv = src3._FillValue if hasattr(src3, '_FillValue') else None
        dst = out.createVariable(vname, src3.dtype, src3.dimensions,
                                 fill_value=fv, zlib=True, complevel=4)
        copy_var_attrs(src3, dst)
        if src2 is not None:
            if src3.ndim == 1:
                dst[:n_spill] = src2[SPLIT:]
                dst[n_spill:] = src3[:]
            else:
                dst[:n_spill, :] = src2[SPLIT:, :]
                dst[n_spill:, :] = src3[:, :]
        else:
            # Variable not in 20240202 file — fill spillover with fill value
            if src3.ndim == 1:
                dst[n_spill:] = src3[:]
            else:
                dst[n_spill:, :] = src3[:, :]

    # --- Range variable ---
    src = ds03.variables['range']
    dst = out.createVariable('range', src.dtype, src.dimensions, zlib=True, complevel=4)
    copy_var_attrs(src, dst)
    dst[:] = src[:]

    # --- Sweep variables ---
    for vname in ('sweep_start_ray_index', 'sweep_end_ray_index',
                  'sweep_number', 'fixed_angle', 'sweep_mode'):
        src = ds03.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        if vname == 'sweep_start_ray_index':
            dst[:] = new_ss
        elif vname == 'sweep_end_ray_index':
            dst[:] = new_se
        elif vname == 'sweep_number':
            dst[:] = new_sn
        elif vname == 'fixed_angle':
            dst[:] = new_fa
        elif vname == 'sweep_mode':
            dst[:] = new_sm

    # --- Scalar / frequency variables ---
    for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                  'frequency', 'azimuth_correction'):
        if vname not in ds03.variables:
            continue
        src = ds03.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # --- String character variables ---
    write_char_var(out, 'time_coverage_start', tcs03_str, sl,
                   src_var=ds03.variables.get('time_coverage_start'))
    write_char_var(out, 'time_coverage_end', tce03_str, sl,
                   src_var=ds03.variables.get('time_coverage_end'))
    for vname in ('time_reference',):
        if vname not in ds03.variables:
            continue
        src = ds03.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # featureType
    out.featureType = 'timeSeriesProfile'

    # --- History and last_revised_date ---
    existing_hist = getattr(out, 'history', '')
    new_entry = _history_entry(
        f'Prepended {n_spill} spillover rays (2024-02-03T00:00:03Z - '
        f'2024-02-03T09:34:03Z) from ncas-mobile-ka-band-radar-1_cao_'
        f'20240202-143823_vpt_l1_v1.0.1.nc as new sweep 0; '
        f'existing sweeps renumbered 1-{n_sweeps03}'
    )
    out.history = (new_entry + '\n' + existing_hist).strip()
    out.last_revised_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

ds02.close()
ds03.close()

print(f"  Written: {OUT03}")
print("\nDone.  Next step: remove or archive the old 20240203-123024_vpt file if no longer needed.")
print(f"  Old file backup: {bak03}")
