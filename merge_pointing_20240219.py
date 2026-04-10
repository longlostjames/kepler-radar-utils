"""
merge_pointing_20240219.py
--------------------------
One-off fix for 20240219: campaign processing created 42 individual single-sweep
pointing files instead of one multi-sweep file.  This script merges them into a
single CfRadial file with 42 sweeps.

Operations
----------
1. Read all *20240219*pointing*.nc files from the L1 directory, sorted by name
   (= time order).
2. Write a merged file named after the first file's timestamp:
     ncas-mobile-ka-band-radar-1_cao_20240219-083218_pointing_l1_v1.0.1.nc
3. Move the 42 original single-sweep files into a backup subdirectory
     <L1_DATE_DIR>/pointing_originals_20240219/
   so they are preserved but no longer interfere with processing.

File details
------------
42 input files covering 08:32:18 – 08:36:10 UTC
Each file: 1 sweep, 15-17 rays, 733 gates
Merged output: 42 sweeps, variable rays per sweep, 733 gates
"""

import datetime
import glob
import os
import shutil
import socket

import cftime
import netCDF4 as nc4
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = (
    '/gws/pw/j07/ncas_obs_vol2/cao/processing/'
    'ncas-mobile-ka-band-radar-1/ccrest-m/L1_v1.0.1'
)
DATE_DIR = os.path.join(BASE, '20240219')
BACKUP_DIR = os.path.join(DATE_DIR, 'pointing_originals_20240219')

SRC_FILES = sorted(glob.glob(os.path.join(DATE_DIR, '*20240219*pointing*.nc')))
assert len(SRC_FILES) > 0, "No pointing files found for 20240219"
print(f"Found {len(SRC_FILES)} pointing files to merge")

# Output uses timestamp of the first file
first_basename = os.path.basename(SRC_FILES[0])
OUT_FILE = os.path.join(DATE_DIR, first_basename)   # same name as first file


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


def write_char_var(out_ds, vname, s, sl, src_var=None):
    dst = out_ds.createVariable(vname, 'S1', ('string_length',))
    if src_var is not None:
        copy_var_attrs(src_var, dst)
    dst[:] = str_to_char_arr(s, sl)
    return dst


# ---------------------------------------------------------------------------
# Read all source datasets
# ---------------------------------------------------------------------------
datasets = [nc4.Dataset(f, 'r') for f in SRC_FILES]
ds0 = datasets[0]

n_sweeps = len(datasets)
sl = ds0.dimensions['string_length'].size
n_range = ds0.dimensions['range'].size
freq_size = ds0.dimensions['frequency'].size

# Compute per-sweep ray counts and cumulative ray offsets
ray_counts = [int(ds.dimensions['time'].size) for ds in datasets]
n_rays_total = sum(ray_counts)
sweep_start = []
sweep_end = []
offset = 0
for rc in ray_counts:
    sweep_start.append(offset)
    sweep_end.append(offset + rc - 1)
    offset += rc

print(f"Total sweeps : {n_sweeps}")
print(f"Ray counts   : min={min(ray_counts)} max={max(ray_counts)}")
print(f"Total rays   : {n_rays_total}")
print(f"Gates        : {n_range}")

# Determine time_coverage_start / end from first and last files
t0 = datasets[0].variables['time']
t1 = datasets[-1].variables['time']
tcs = cftime.num2pydate(t0[0], t0.units).strftime('%Y-%m-%dT%H:%M:%SZ')
tce = cftime.num2pydate(t1[-1], t1.units).strftime('%Y-%m-%dT%H:%M:%SZ')
print(f"time_coverage_start: {tcs}")
print(f"time_coverage_end  : {tce}")

# ---------------------------------------------------------------------------
# Write merged file (to a temp path, then rename)
# ---------------------------------------------------------------------------
tmp_out = OUT_FILE + '.tmp'

print(f"\nWriting merged file: {OUT_FILE}")

with nc4.Dataset(tmp_out, 'w', format='NETCDF4') as out:
    # --- Dimensions ---
    out.createDimension('time', None)          # unlimited
    out.createDimension('range', n_range)
    out.createDimension('sweep', n_sweeps)
    out.createDimension('string_length', sl)
    out.createDimension('frequency', freq_size)

    # --- Global attributes ---
    now_str = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
    copy_global_attrs(ds0, out, overrides={
        'time_coverage_start': tcs,
        'time_coverage_end': tce,
        'history': _history_entry(
            f'Merged {n_sweeps} single-sweep pointing files into one multi-sweep file'
        ),
        'last_revised_date': now_str,
        'date_created': now_str,
    })

    # --- Range (copy from first file) ---
    src = ds0.variables['range']
    dst = out.createVariable('range', src.dtype, ('range',), zlib=True, complevel=4)
    copy_var_attrs(src, dst)
    dst[:] = src[:]

    # --- Scalar variables (copy from first file) ---
    for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                  'azimuth_correction'):
        if vname in ds0.variables:
            src = ds0.variables[vname]
            dst = out.createVariable(vname, src.dtype, src.dimensions)
            copy_var_attrs(src, dst)
            dst.assignValue(src[...])

    # --- Frequency variable ---
    if 'frequency' in ds0.variables:
        src = ds0.variables['frequency']
        dst = out.createVariable('frequency', src.dtype, ('frequency',),
                                 zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # --- Ray-indexed variables (concatenate across files) ---
    ray_dim_vars = [
        v for v in ds0.variables
        if ds0.variables[v].dimensions and ds0.variables[v].dimensions[0] == 'time'
    ]
    for vname in ray_dim_vars:
        src0 = ds0.variables[vname]
        fv = src0._FillValue if hasattr(src0, '_FillValue') else None
        if vname == 'time':
            dst = out.createVariable(vname, src0.dtype, src0.dimensions,
                                     fill_value=fv)
        else:
            dst = out.createVariable(vname, src0.dtype, src0.dimensions,
                                     fill_value=fv, zlib=True, complevel=4)
        copy_var_attrs(src0, dst)

        if src0.ndim == 1:
            data = np.concatenate([ds.variables[vname][:] for ds in datasets])
            dst[:] = data
        else:
            data = np.concatenate([ds.variables[vname][:, :] for ds in datasets], axis=0)
            dst[:, :] = data

    # --- Sweep-indexed variables ---
    # sweep_number: 0 .. n_sweeps-1
    src = ds0.variables['sweep_number']
    dst = out.createVariable('sweep_number', src.dtype, ('sweep',),
                             zlib=True, complevel=4)
    copy_var_attrs(src, dst)
    dst[:] = np.arange(n_sweeps, dtype=src.dtype)

    # sweep_start_ray_index / sweep_end_ray_index
    for vname, data in [
        ('sweep_start_ray_index', np.array(sweep_start, dtype=np.int32)),
        ('sweep_end_ray_index',   np.array(sweep_end,   dtype=np.int32)),
    ]:
        src = ds0.variables[vname]
        dst = out.createVariable(vname, src.dtype, ('sweep',),
                                 zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = data

    # fixed_angle and target_scan_rate: concatenate from each file's sweep dim
    for vname in ('fixed_angle', 'target_scan_rate'):
        if vname in ds0.variables:
            src0 = ds0.variables[vname]
            dst = out.createVariable(vname, src0.dtype, ('sweep',),
                                     zlib=True, complevel=4)
            copy_var_attrs(src0, dst)
            dst[:] = np.concatenate([ds.variables[vname][:] for ds in datasets])

    # sweep_mode: (sweep, string_length) — take from each file's single sweep
    src0 = ds0.variables['sweep_mode']
    dst = out.createVariable('sweep_mode', 'S1', ('sweep', 'string_length'))
    copy_var_attrs(src0, dst)
    mode_arr = np.concatenate([ds.variables['sweep_mode'][:] for ds in datasets], axis=0)
    dst[:] = mode_arr

    # --- String time variables (from first / last file) ---
    write_char_var(out, 'time_coverage_start', tcs, sl,
                   src_var=ds0.variables['time_coverage_start'])
    write_char_var(out, 'time_coverage_end', tce, sl,
                   src_var=ds0.variables['time_coverage_end'])
    write_char_var(out, 'time_reference', tcs, sl,
                   src_var=ds0.variables['time_reference'])

print("  Written OK")

# ---------------------------------------------------------------------------
# Close all source datasets before moving files
# ---------------------------------------------------------------------------
for ds in datasets:
    ds.close()

# ---------------------------------------------------------------------------
# Validate: open the merged file with pyart
# ---------------------------------------------------------------------------
print("\nValidating merged file with pyart...")
import pyart
radar = pyart.io.read_cfradial(tmp_out)
print(f"  nsweeps = {radar.nsweeps}")
print(f"  nrays   = {radar.nrays}")
print(f"  ngates  = {radar.ngates}")
assert radar.nsweeps == n_sweeps, "Sweep count mismatch!"
assert radar.nrays == n_rays_total, "Ray count mismatch!"
print("  Validation passed")

# ---------------------------------------------------------------------------
# Back up and install merged file
# ---------------------------------------------------------------------------
os.makedirs(BACKUP_DIR, exist_ok=True)
print(f"\nMoving {n_sweeps} original files to: {BACKUP_DIR}")
for f in SRC_FILES:
    shutil.move(f, os.path.join(BACKUP_DIR, os.path.basename(f)))
    print(f"  Moved: {os.path.basename(f)}")

os.rename(tmp_out, OUT_FILE)
print(f"\nMerged file installed: {OUT_FILE}")
print("Done.")
