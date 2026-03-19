"""
fix_midnight_crossing_vpt_20171212.py
--------------------------------------
One-off fix for ncas-mobile-ka-band-radar-1_cao_20171212-002916_vpt_l1_v1.0.0.nc
which runs past midnight into 2017-12-13.

Operations
----------
1. Trim the 20171212 file to rays 0-10614 (up to midnight).
   Sweep 3 is truncated at its new end ray; all other sweeps are unchanged.
   The original file is backed up as *_original.nc before overwriting.
   History attribute is updated to record the trim.

2. Create a merged 20171213 file by prepending the spillover rays 10615-11361
   (00:00:07 - 02:32:54) to the existing 20171213-023315 file.
   The merged file is written as ...20171213-000007_vpt_l1_v1.0.0.nc.
   The old 023315 file is backed up as *_original.nc.
   History attribute is updated to record the merge.

NOTE: The *_original.nc backups must exist and contain the original (pre-trim)
data for this script to run correctly.  If they have been removed, you cannot
re-run this script.
"""

import datetime
import os
import shutil
import socket

import netCDF4 as nc4
import numpy as np


def _history_entry(msg):
    """Return a new history line matching the existing format: 'DDD Mon DD HH:MM:SS YYYY - <msg>'"""
    now = datetime.datetime.utcnow()
    user = os.environ.get('USER', 'unknown')
    host = socket.gethostname()
    timestamp = now.strftime('%a %b %d %H:%M:%S %Y')
    return f'{timestamp} - user:{user} machine:{host} {msg}'

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/picasso/L1_v1.0.0'

SRC12  = os.path.join(BASE, '20171212', 'ncas-mobile-ka-band-radar-1_cao_20171212-002916_vpt_l1_v1.0.0.nc')
SRC13  = os.path.join(BASE, '20171213', 'ncas-mobile-ka-band-radar-1_cao_20171213-023315_vpt_l1_v1.0.0.nc')
OUT13  = os.path.join(BASE, '20171213', 'ncas-mobile-ka-band-radar-1_cao_20171213-000007_vpt_l1_v1.0.0.nc')

# Ray index of the first ray on 2017-12-13 inside the 20171212 file
SPLIT = 10615

# ---------------------------------------------------------------------------
# Helper: copy variable attributes
# ---------------------------------------------------------------------------
def copy_var_attrs(src_var, dst_var):
    for attr in src_var.ncattrs():
        if attr != '_FillValue':
            setattr(dst_var, attr, getattr(src_var, attr))


def copy_global_attrs(src_ds, dst_ds, overrides=None):
    for attr in src_ds.ncattrs():
        val = overrides.get(attr, getattr(src_ds, attr)) if overrides else getattr(src_ds, attr)
        setattr(dst_ds, attr, val)
    if overrides:
        for k, v in overrides.items():
            if k not in src_ds.ncattrs():
                setattr(dst_ds, k, v)


# ---------------------------------------------------------------------------
# 1.  Trim the 20171212 file (rays 0 .. SPLIT-1)
# ---------------------------------------------------------------------------
print("=== Step 1: Trimming 20171212 file ===")

bak12 = SRC12.replace('.nc', '_original.nc')
if not os.path.exists(bak12):
    shutil.copy2(SRC12, bak12)
    print(f"  Backup: {bak12}")
else:
    print(f"  Backup already exists, skipping: {bak12}")

ds12 = nc4.Dataset(SRC12, 'r')
t12  = ds12.variables['time']
times12 = nc4.num2date(t12[:], t12.units, only_use_cftime_datetimes=False)

# Compute trimmed time_coverage_end
last_ray_time = times12[SPLIT - 1]
tce12 = last_ray_time.strftime('%Y-%m-%dT%H:%M:%SZ')

# Sweep arrays (4 sweeps in original)
ss12 = ds12.variables['sweep_start_ray_index'][:]
se12 = ds12.variables['sweep_end_ray_index'][:]
# Find sweeps that have any rays < SPLIT
n_sweeps12 = int(ds12.dimensions['sweep'].size)
new_ss = []
new_se = []
keep_sweeps = []
for i in range(n_sweeps12):
    s, e = int(ss12[i]), int(se12[i])
    if s >= SPLIT:
        break  # sweep entirely after midnight - drop
    new_ss.append(s)
    new_se.append(min(e, SPLIT - 1))  # truncate if straddles midnight
    keep_sweeps.append(i)

n_new_sweeps12 = len(keep_sweeps)
print(f"  Keeping {n_new_sweeps12} sweeps, {SPLIT} rays")

tmp12 = SRC12 + '.tmp'
with nc4.Dataset(tmp12, 'w', format='NETCDF4') as out:
    # Dimensions
    out.createDimension('time', None)          # unlimited
    out.createDimension('range', ds12.dimensions['range'].size)
    out.createDimension('sweep', n_new_sweeps12)
    out.createDimension('string_length', ds12.dimensions['string_length'].size)
    out.createDimension('frequency', ds12.dimensions['frequency'].size)

    # Global attrs
    copy_global_attrs(ds12, out, overrides={'time_coverage_end': tce12})

    # Ray-indexed variables (time dim)
    ray_vars = [v for v in ds12.variables
                if ds12.variables[v].dimensions and ds12.variables[v].dimensions[0] == 'time']
    for vname in ray_vars:
        src = ds12.variables[vname]
        fv = src._FillValue if hasattr(src, '_FillValue') else None
        if vname == 'time':
            dst = out.createVariable(vname, src.dtype, src.dimensions, fill_value=fv)
        else:
            dst = out.createVariable(vname, src.dtype, src.dimensions,
                                     fill_value=fv, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        if src.ndim == 1:
            dst[:] = src[:SPLIT]
        else:
            dst[:] = src[:SPLIT, :]

    # Range variable
    for vname in ('range',):
        src = ds12.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # Sweep variables
    for vname in ('sweep_number', 'fixed_angle', 'sweep_start_ray_index',
                  'sweep_end_ray_index', 'sweep_mode'):
        src = ds12.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions, zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        if vname == 'sweep_start_ray_index':
            dst[:] = np.array(new_ss, dtype=src.dtype)
        elif vname == 'sweep_end_ray_index':
            dst[:] = np.array(new_se, dtype=src.dtype)
        elif src.ndim == 1:
            dst[:] = src[keep_sweeps]
        else:
            dst[:] = src[keep_sweeps, :]

    # Scalar / frequency variables
    for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                  'frequency', 'azimuth_correction'):
        src = ds12.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # String char variables (time_coverage_start/end, time_reference)
    sl = ds12.dimensions['string_length'].size
    for vname in ('time_coverage_start', 'time_reference'):
        src = ds12.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]
    # Override time_coverage_end
    src_tce = ds12.variables['time_coverage_end']
    dst_tce = out.createVariable('time_coverage_end', src_tce.dtype, src_tce.dimensions)
    copy_var_attrs(src_tce, dst_tce)
    tce_arr = np.zeros(sl, dtype='S1')
    tce_bytes = tce12.encode('ascii')
    for idx, b in enumerate(tce_bytes):
        tce_arr[idx] = bytes([b])
    dst_tce[:] = tce_arr

    # Update history and last_revised_date in the trimmed file
    existing_hist = getattr(out, 'history', '')
    new_entry = _history_entry(
        f'Trimmed to {SPLIT} rays (removed post-midnight data; '
        f'see ncas-mobile-ka-band-radar-1_cao_20171213-000007_vpt_l1_v1.0.0.nc for continuation)'
    )
    out.history = (new_entry + '\n' + existing_hist).strip()
    out.last_revised_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

ds12.close()

os.replace(tmp12, SRC12)
print(f"  Written: {SRC12}")


# ---------------------------------------------------------------------------
# 2.  Create merged 20171213 file
# ---------------------------------------------------------------------------
print("\n=== Step 2: Creating merged 20171213 file ===")

bak13 = SRC13.replace('.nc', '_original.nc')
if not os.path.exists(bak13):
    shutil.copy2(SRC13, bak13)
    print(f"  Backup: {bak13}")
else:
    print(f"  Backup already exists, skipping: {bak13}")

ds12 = nc4.Dataset(bak12, 'r')   # read original (pre-trim) 20171212 file
ds13 = nc4.Dataset(SRC13, 'r')

t12   = ds12.variables['time']
t13   = ds13.variables['time']
times12 = nc4.num2date(t12[:], t12.units, only_use_cftime_datetimes=False)
times13 = nc4.num2date(t13[:], t13.units, only_use_cftime_datetimes=False)

# Spillover rays in 20171212 file: SPLIT .. end
n_spill = int(ds12.dimensions['time'].size) - SPLIT   # 747
n_13    = int(ds13.dimensions['time'].size)             # 4850
n_total = n_spill + n_13

# New time values: units = seconds since 2017-12-13T00:00:00Z
# 20171212 time unit: seconds since 2017-12-12T00:00:00Z  -> subtract 86400
spill_times = t12[SPLIT:] - 86400.0
exist_times = t13[:]

new_time_units = 'seconds since 2017-12-13T00:00:00Z'

# Sweep structure for merged file
# Sweep 0: spillover (from 20171212 sweep 3, rays SPLIT..end)
#   -> ray indices 0 .. n_spill-1 in new file
# Sweeps 1..5: existing 20171213 sweeps with ray indices offset by n_spill
ss12 = ds12.variables['sweep_start_ray_index'][:]
se12 = ds12.variables['sweep_end_ray_index'][:]
fa12 = ds12.variables['fixed_angle'][:]
sm12 = ds12.variables['sweep_mode'][:]
sn12 = ds12.variables['sweep_number'][:]

ss13 = ds13.variables['sweep_start_ray_index'][:]
se13 = ds13.variables['sweep_end_ray_index'][:]
fa13 = ds13.variables['fixed_angle'][:]
sm13 = ds13.variables['sweep_mode'][:]
sn13 = ds13.variables['sweep_number'][:]

n_sweeps13 = int(ds13.dimensions['sweep'].size)
n_new_sweeps = 1 + n_sweeps13
sl = ds12.dimensions['string_length'].size

new_ss = np.zeros(n_new_sweeps, dtype=ss13.dtype)
new_se = np.zeros(n_new_sweeps, dtype=se13.dtype)
new_fa = np.zeros(n_new_sweeps, dtype=fa13.dtype)
new_sm = np.zeros((n_new_sweeps, sl), dtype=sm13.dtype)
new_sn = np.arange(n_new_sweeps, dtype=sn13.dtype)

# Sweep 0: spillover
new_ss[0] = 0
new_se[0] = n_spill - 1
new_fa[0] = fa12[-1]   # fixed angle from last sweep of 20171212
new_sm[0, :] = sm12[-1, :]

# Sweeps 1..n_sweeps13
for i in range(n_sweeps13):
    new_ss[1 + i] = int(ss13[i]) + n_spill
    new_se[1 + i] = int(se13[i]) + n_spill
    new_fa[1 + i] = fa13[i]
    new_sm[1 + i, :] = sm13[i, :]

# Time coverage strings for merged file
tcs13_str = times12[SPLIT].strftime('%Y-%m-%dT%H:%M:%SZ')
tce13_str = times13[-1].strftime('%Y-%m-%dT%H:%M:%SZ')

def str_to_char_arr(s, length):
    arr = np.zeros(length, dtype='S1')
    for i, b in enumerate(s.encode('ascii')):
        arr[i] = bytes([b])
    return arr

print(f"  Merged file: {n_total} rays, {n_new_sweeps} sweeps")
print(f"  time_coverage_start: {tcs13_str}")
print(f"  time_coverage_end:   {tce13_str}")

with nc4.Dataset(OUT13, 'w', format='NETCDF4') as out:
    # Dimensions
    out.createDimension('time', None)
    out.createDimension('range', ds13.dimensions['range'].size)
    out.createDimension('sweep', n_new_sweeps)
    out.createDimension('string_length', sl)
    out.createDimension('frequency', ds13.dimensions['frequency'].size)

    # Global attrs from 20171213 file, with updated coverage bounds + filename
    overrides = {
        'time_coverage_start': tcs13_str,
        'time_coverage_end':   tce13_str,
    }
    copy_global_attrs(ds13, out, overrides=overrides)

    # --- time variable ---
    t_out = out.createVariable('time', t13.dtype, ('time',))
    copy_var_attrs(t13, t_out)
    t_out.units = new_time_units
    t_out[:] = np.concatenate([spill_times, exist_times])

    # --- range ---
    r_src = ds13.variables['range']
    r_out = out.createVariable('range', r_src.dtype, ('range',), zlib=True, complevel=4)
    copy_var_attrs(r_src, r_out)
    r_out[:] = r_src[:]

    # --- ray-indexed vars (excluding time) ---
    ray_vars_1d = ['azimuth', 'elevation', 'prt', 'pulse_width', 'radar_measured_transmit_power_h']
    ray_vars_2d = ['DBZ', 'VEL', 'WIDTH', 'LDR', 'SNR']

    for vname in ray_vars_1d:
        src12v = ds12.variables[vname]
        src13v = ds13.variables[vname]
        fv = src13v._FillValue if hasattr(src13v, '_FillValue') else None
        dst = out.createVariable(vname, src13v.dtype, src13v.dimensions,
                                 fill_value=fv, zlib=True, complevel=4)
        copy_var_attrs(src13v, dst)
        dst[:] = np.concatenate([src12v[SPLIT:], src13v[:]])

    for vname in ray_vars_2d:
        src12v = ds12.variables[vname]
        src13v = ds13.variables[vname]
        fv = src13v._FillValue if hasattr(src13v, '_FillValue') else None
        dst = out.createVariable(vname, src13v.dtype, src13v.dimensions,
                                 fill_value=fv, zlib=True, complevel=4)
        copy_var_attrs(src13v, dst)
        dst[:] = np.concatenate([src12v[SPLIT:, :], src13v[:, :]], axis=0)

    # --- sweep variables ---
    for vname, arr in [('sweep_start_ray_index', new_ss),
                       ('sweep_end_ray_index',   new_se),
                       ('fixed_angle',           new_fa),
                       ('sweep_number',          new_sn)]:
        src = ds13.variables[vname]
        dst = out.createVariable(vname, src.dtype, ('sweep',), zlib=True, complevel=4)
        copy_var_attrs(src, dst)
        dst[:] = arr

    sm_out = out.createVariable('sweep_mode', sm13.dtype, ('sweep', 'string_length'), zlib=True, complevel=4)
    copy_var_attrs(ds13.variables['sweep_mode'], sm_out)
    sm_out[:] = new_sm

    # --- scalar variables ---
    for vname in ('latitude', 'longitude', 'altitude', 'volume_number',
                  'frequency', 'azimuth_correction'):
        src = ds13.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # --- string char variables ---
    for vname, val_str in [('time_coverage_start', tcs13_str),
                           ('time_coverage_end',   tce13_str)]:
        src = ds13.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = str_to_char_arr(val_str, sl)

    for vname in ('time_reference',):
        src = ds13.variables[vname]
        dst = out.createVariable(vname, src.dtype, src.dimensions)
        copy_var_attrs(src, dst)
        dst[:] = src[:]

    # Update history and last_revised_date in the merged file
    existing_hist = getattr(ds13, 'history', '')
    new_entry = _history_entry(
        f'Merged: prepended {n_spill} spillover rays (00:00:07-02:32:54 from '
        f'ncas-mobile-ka-band-radar-1_cao_20171212-002916_vpt_l1_v1.0.0.nc) '
        f'to former ncas-mobile-ka-band-radar-1_cao_20171213-023315_vpt_l1_v1.0.0.nc; '
        f'{n_new_sweeps} sweeps total'
    )
    out.history = (new_entry + '\n' + existing_hist).strip()
    out.last_revised_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

ds12.close()
ds13.close()

print(f"  Written: {OUT13}")
print("\nDone. Old files backed up as *_original.nc")
print(f"  {bak12}")
print(f"  {bak13}")
print("\nVerify outputs, then you can delete the *_original.nc backups.")
