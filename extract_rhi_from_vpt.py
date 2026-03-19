"""
extract_rhi_from_vpt.py
-----------------------
Extract rays flagged antenna_transition=1 from a VPT CF-Radial L1 file and
write them to a new CF-Radial file with the correct scan type.  Also writes
an amended VPT file with those rays removed.

In VPT files, antenna_transition=1 is set for any ray with elevation < 85°.
This catches non-vertical scan segments (e.g. RHI or manual scans) that were
embedded in a VPT file and incorrectly flagged as transitions.

Output files
------------
1. Extracted scan file  – new filename: 'vpt' replaced by scan type, timestamp
   updated to the start of the first extracted ray.  Written to the same
   directory as the input.

2. Amended VPT file     – original file backed up as *_original.nc, then the
   original path is overwritten with the transition-free version.  Original
   sweep boundaries are preserved; sweeps that become empty are dropped.

Usage
-----
    python extract_rhi_from_vpt.py <vpt_cfradial_file> [--scan-type rhi]
"""

import argparse
import datetime
import os
import re
import shutil
import socket

import netCDF4 as nc4
import numpy as np


# ---------------------------------------------------------------------------
# Helpers  (same style as fix_midnight_crossing_vpt_20171212.py)
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
    for attr in src_ds.ncattrs():
        val = overrides.get(attr, getattr(src_ds, attr)) if overrides else getattr(src_ds, attr)
        setattr(dst_ds, attr, val)
    if overrides:
        for k, v in overrides.items():
            if k not in src_ds.ncattrs():
                setattr(dst_ds, k, v)


def str_to_char_arr(s, length):
    arr = np.zeros(length, dtype='S1')
    for i, b in enumerate(s.encode('ascii')):
        arr[i] = bytes([b])
    return arr


def find_contiguous_groups(mask):
    """Return list of (start, end) index pairs for contiguous True runs."""
    groups = []
    in_group = False
    for i, val in enumerate(mask):
        if val and not in_group:
            start = i
            in_group = True
        elif not val and in_group:
            groups.append((start, i - 1))
            in_group = False
    if in_group:
        groups.append((start, len(mask) - 1))
    return groups


# ---------------------------------------------------------------------------
# Shared CF-Radial writer
# ---------------------------------------------------------------------------

def _write_cfradial(ds, out_path,
                    all_idx, new_ss, new_se, fixed_angles,
                    sweep_mode_rows,
                    scan_type_attr, tcs_str, tce_str,
                    ray_vars, range_vars, scalar_vars,
                    history_msg, phase_sequence=None):
    """
    Write a CF-Radial file selecting a subset of rays from *ds*.

    Parameters
    ----------
    ds            : open nc4.Dataset (source, read-only)
    out_path      : destination path (will be created / overwritten)
    all_idx       : 1-D array of ray indices to copy from ds (in order)
    new_ss/se     : lists of new sweep start/end indices (re-indexed from 0)
    fixed_angles  : list of fixed_angle values per sweep
    sweep_mode_rows : (n_sweeps, string_length) uint8/S1 array for sweep_mode
    scan_type_attr : string for the global 'scan_type' attribute
    tcs_str/tce_str : ISO-8601 time coverage strings
    ray_vars      : variable names whose first dim is 'time'
    range_vars    : variable names whose first dim is 'range'
    scalar_vars   : all other variable names (scalars / frequency)
    history_msg   : string to prepend to history
    """
    n_range = int(ds.dimensions['range'].size)
    sl = int(ds.dimensions['string_length'].size)
    has_freq_dim = 'frequency' in ds.dimensions
    n_sweeps_out = len(new_ss)

    sweep_handled = {'sweep_number', 'fixed_angle',
                     'sweep_start_ray_index', 'sweep_end_ray_index', 'sweep_mode'}
    str_handled = {'time_coverage_start', 'time_coverage_end', 'time_reference'}

    ss_src = ds.variables.get('sweep_start_ray_index')
    se_src = ds.variables.get('sweep_end_ray_index')
    fa_src = ds.variables.get('fixed_angle')
    sn_src = ds.variables.get('sweep_number')
    sm_src = ds.variables.get('sweep_mode')

    ss_dtype = ss_src.dtype if ss_src is not None else np.int32
    se_dtype = se_src.dtype if se_src is not None else np.int32
    fa_dtype = fa_src.dtype if fa_src is not None else np.float32
    sn_dtype = sn_src.dtype if sn_src is not None else np.int32
    sm_dtype = sm_src.dtype if sm_src is not None else 'S1'

    with nc4.Dataset(out_path, 'w', format='NETCDF4') as out:

        # Dimensions
        out.createDimension('time', None)
        out.createDimension('range', n_range)
        out.createDimension('sweep', n_sweeps_out)
        out.createDimension('string_length', sl)
        if has_freq_dim:
            out.createDimension('frequency',
                                int(ds.dimensions['frequency'].size))

        # Global attributes
        overrides = {
            'time_coverage_start': tcs_str,
            'time_coverage_end': tce_str,
            'scan_type': scan_type_attr,
        }
        if phase_sequence is not None:
            overrides['phase_sequence'] = phase_sequence
        copy_global_attrs(ds, out, overrides=overrides)

        # ---- Ray-indexed variables ----
        for vname in ray_vars:
            src = ds.variables[vname]
            fv = src._FillValue if hasattr(src, '_FillValue') else None
            if vname == 'time':
                dst = out.createVariable(vname, src.dtype, src.dimensions,
                                         fill_value=fv)
            else:
                dst = out.createVariable(vname, src.dtype, src.dimensions,
                                         fill_value=fv, zlib=True, complevel=4)
            copy_var_attrs(src, dst)
            if src.ndim == 1:
                dst[:] = src[all_idx]
            else:
                dst[:] = src[all_idx, :]

        # ---- Range variables ----
        for vname in range_vars:
            src = ds.variables[vname]
            fv = src._FillValue if hasattr(src, '_FillValue') else None
            dst = out.createVariable(vname, src.dtype, src.dimensions,
                                     fill_value=fv, zlib=True, complevel=4)
            copy_var_attrs(src, dst)
            dst[:] = src[:]

        # ---- Sweep variables ----
        dst_ss = out.createVariable('sweep_start_ray_index', ss_dtype,
                                    ('sweep',), zlib=True, complevel=4)
        if ss_src is not None:
            copy_var_attrs(ss_src, dst_ss)
        dst_ss[:] = np.array(new_ss, dtype=ss_dtype)

        dst_se = out.createVariable('sweep_end_ray_index', se_dtype,
                                    ('sweep',), zlib=True, complevel=4)
        if se_src is not None:
            copy_var_attrs(se_src, dst_se)
        dst_se[:] = np.array(new_se, dtype=se_dtype)

        dst_fa = out.createVariable('fixed_angle', fa_dtype,
                                    ('sweep',), zlib=True, complevel=4)
        if fa_src is not None:
            copy_var_attrs(fa_src, dst_fa)
        dst_fa[:] = np.array(fixed_angles, dtype=fa_dtype)

        dst_sn = out.createVariable('sweep_number', sn_dtype,
                                    ('sweep',), zlib=True, complevel=4)
        if sn_src is not None:
            copy_var_attrs(sn_src, dst_sn)
        dst_sn[:] = np.arange(n_sweeps_out, dtype=sn_dtype)

        dst_sm = out.createVariable('sweep_mode', sm_dtype,
                                    ('sweep', 'string_length'),
                                    zlib=True, complevel=4)
        if sm_src is not None:
            copy_var_attrs(sm_src, dst_sm)
        dst_sm[:] = sweep_mode_rows

        # ---- String char variables ----
        for vname, val_str in [('time_coverage_start', tcs_str),
                                ('time_coverage_end', tce_str)]:
            if vname in ds.variables:
                src = ds.variables[vname]
                dst = out.createVariable(vname, src.dtype, src.dimensions)
                copy_var_attrs(src, dst)
                dst[:] = str_to_char_arr(val_str, sl)

        if 'time_reference' in ds.variables:
            src = ds.variables['time_reference']
            dst = out.createVariable('time_reference', src.dtype,
                                     src.dimensions)
            copy_var_attrs(src, dst)
            dst[:] = src[:]

        # ---- Scalar / frequency variables ----
        for vname in scalar_vars:
            src = ds.variables[vname]
            dst = out.createVariable(vname, src.dtype, src.dimensions)
            copy_var_attrs(src, dst)
            dst[:] = src[:]

        # ---- History ----
        existing_hist = getattr(ds, 'history', '')
        new_entry = _history_entry(history_msg)
        out.history = (new_entry + '\n' + existing_hist).strip()
        out.last_revised_date = (
            datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Extract antenna_transition=1 rays from a VPT CF-Radial '
                    'file into a new scan file, and write an amended VPT file '
                    'with those rays removed.'
    )
    parser.add_argument('vpt_file', help='Input VPT CF-Radial L1 NetCDF file')
    parser.add_argument(
        '--scan-type', default='rhi',
        help='Scan type label for the extracted file and its sweep_mode '
             '(default: rhi).  Use "manual_rhi" for manually-pointed RHI scans.'
    )
    args = parser.parse_args()

    src_path = args.vpt_file
    scan_type = args.scan_type.lower()

    if not os.path.exists(src_path):
        raise FileNotFoundError(f'Input file not found: {src_path}')

    ds = nc4.Dataset(src_path, 'r')

    n_rays = int(ds.dimensions['time'].size)
    sl = int(ds.dimensions['string_length'].size)

    # ------------------------------------------------------------------
    # Find transition rays
    # ------------------------------------------------------------------
    if 'antenna_transition' in ds.variables:
        at = ds.variables['antenna_transition'][:]
    else:
        print('  WARNING: antenna_transition variable not found; '
              'deriving from elevation < 85°')
        at = (ds.variables['elevation'][:] < 85.0).astype(np.int8)

    trans_mask = at.astype(bool)
    n_trans = int(np.sum(trans_mask))

    if n_trans == 0:
        print('No antenna_transition=1 rays found. Nothing to do.')
        ds.close()
        return

    print(f'Found {n_trans} antenna_transition=1 rays out of {n_rays} total.')

    elev = ds.variables['elevation'][:]
    azim = ds.variables['azimuth'][:]
    time_var = ds.variables['time']
    times = nc4.num2date(time_var[:], time_var.units,
                         only_use_cftime_datetimes=False)

    # ------------------------------------------------------------------
    # Classify source variables (shared between both output files)
    # ------------------------------------------------------------------
    all_var_names = list(ds.variables.keys())

    ray_vars = [v for v in all_var_names
                if ds.variables[v].dimensions
                and ds.variables[v].dimensions[0] == 'time']
    range_vars = [v for v in all_var_names
                  if ds.variables[v].dimensions
                  and ds.variables[v].dimensions[0] == 'range']
    sweep_handled = {'sweep_number', 'fixed_angle',
                     'sweep_start_ray_index', 'sweep_end_ray_index', 'sweep_mode'}
    str_handled = {'time_coverage_start', 'time_coverage_end', 'time_reference'}
    already_handled = set(ray_vars) | set(range_vars) | sweep_handled | str_handled
    scalar_vars = [v for v in all_var_names if v not in already_handled]

    src_dir = os.path.dirname(os.path.abspath(src_path))
    src_base = os.path.basename(src_path)

    # original scan_type global attribute (for the amended VPT output)
    orig_scan_type = getattr(ds, 'scan_type', 'vertical_pointing')
    n_sweeps_orig = int(ds.dimensions['sweep'].size)

    # original phase_sequence (may not exist)
    orig_phase_seq_str = getattr(ds, 'phase_sequence', None)
    orig_phases = None
    if orig_phase_seq_str:
        parts = [p.strip() for p in orig_phase_seq_str.split(',')]
        if len(parts) == n_sweeps_orig:
            orig_phases = parts
        else:
            print(f'  INFO: phase_sequence has {len(parts)} entries but file has '
                  f'{n_sweeps_orig} sweeps; ignoring for amended VPT.')

    # ================================================================
    # 1. Extracted scan files  (one file per contiguous group)
    # ================================================================
    print(f'\n=== Step 1: Writing extracted {scan_type} files ===')

    rhi_groups = find_contiguous_groups(trans_mask)
    n_rhi_sweeps = len(rhi_groups)

    print(f'Contiguous group(s) → {n_rhi_sweeps} sweep(s):')
    for i, (s, e) in enumerate(rhi_groups):
        print(f'  Sweep {i}: rays {s}–{e}  ({e - s + 1} rays)  '
              f'el={elev[s:e+1].min():.1f}°–{elev[s:e+1].max():.1f}°  '
              f'az={azim[s:e+1].min():.1f}°–{azim[s:e+1].max():.1f}°')

    if scan_type in ('rhi', 'manual_rhi'):
        rhi_fixed = [float(np.mean(azim[s:e+1])) for s, e in rhi_groups]
    elif scan_type == 'ppi':
        rhi_fixed = [float(np.mean(elev[s:e+1])) for s, e in rhi_groups]
    else:
        rhi_fixed = [float(np.mean(elev[s:e+1])) for s, e in rhi_groups]

    # sweep_mode row (single sweep per file)
    sm_row = np.zeros((1, sl), dtype='S1')
    for j, b in enumerate(scan_type.encode('ascii')[:sl]):
        sm_row[0, j] = bytes([b])

    rhi_basenames = []   # collected for the VPT history message

    for i, (s, e) in enumerate(rhi_groups):
        sweep_idx = np.arange(s, e + 1)
        tcs = times[s].strftime('%Y-%m-%dT%H:%M:%SZ')
        tce = times[e].strftime('%Y-%m-%dT%H:%M:%SZ')

        new_tstamp = times[s].strftime('%Y%m%d-%H%M%S')
        out_base = re.sub(r'\d{8}-\d{6}', new_tstamp, src_base)
        out_base = re.sub(r'_vpt_', f'_{scan_type}_', out_base)
        out_path = os.path.join(src_dir, out_base)

        _write_cfradial(
            ds, out_path,
            sweep_idx, [0], [e - s], [rhi_fixed[i]],
            sm_row, scan_type, tcs, tce,
            ray_vars, range_vars, scalar_vars,
            history_msg=(
                f'Extracted {e - s + 1} antenna_transition=1 ray(s) '
                f'(group {i + 1} of {n_rhi_sweeps}) from {src_base} as {scan_type} scan'
            ),
            phase_sequence=scan_type,
        )
        rhi_basenames.append(out_base)
        print(f'  Written: {out_path}')

    # Reference name used in the VPT history (first file, or summary)
    if n_rhi_sweeps == 1:
        rhi_base = rhi_basenames[0]
    else:
        rhi_base = f'{n_rhi_sweeps} {scan_type} files (first: {rhi_basenames[0]})'

    # ================================================================
    # 2. Amended VPT file  (non-transition rays only)
    # ================================================================
    print(f'\n=== Step 2: Writing amended VPT file ===')

    non_trans_mask = ~trans_mask
    ss_orig = ds.variables['sweep_start_ray_index'][:]
    se_orig = ds.variables['sweep_end_ray_index'][:]
    fa_orig = ds.variables['fixed_angle'][:]
    sm_orig = ds.variables['sweep_mode'][:] if 'sweep_mode' in ds.variables else None

    # For each original sweep, keep only its non-transition rays.
    # Sweeps that are entirely transition rays are dropped.
    vpt_idx_list = []
    vpt_ss, vpt_se, vpt_fa = [], [], []
    vpt_sm_rows = []
    vpt_keep_indices = []   # original sweep indices that survive
    offset = 0

    for i in range(n_sweeps_orig):
        s, e = int(ss_orig[i]), int(se_orig[i])
        keep = np.where(non_trans_mask[s:e + 1])[0] + s
        if len(keep) == 0:
            print(f'  Original sweep {i} (rays {s}–{e}): entirely transition, dropping.')
            continue
        dropped = (e - s + 1) - len(keep)
        if dropped:
            print(f'  Original sweep {i} (rays {s}–{e}): keeping {len(keep)}, '
                  f'dropping {dropped} transition ray(s).')
        else:
            print(f'  Original sweep {i} (rays {s}–{e}): all {len(keep)} rays kept.')
        vpt_keep_indices.append(i)
        vpt_idx_list.append(keep)
        vpt_ss.append(offset)
        vpt_se.append(offset + len(keep) - 1)
        vpt_fa.append(float(fa_orig[i]))
        vpt_sm_rows.append(sm_orig[i] if sm_orig is not None
                           else str_to_char_arr('vertical_pointing', sl))
        offset += len(keep)

    if not vpt_idx_list:
        print('  All rays are transition rays; amended VPT would be empty. Skipping.')
    else:
        vpt_all_idx = np.concatenate(vpt_idx_list)
        n_vpt_rays = len(vpt_all_idx)
        n_vpt_sweeps = len(vpt_ss)
        vpt_sm_arr = np.array(vpt_sm_rows)

        vpt_tcs = times[vpt_all_idx[0]].strftime('%Y-%m-%dT%H:%M:%SZ')
        vpt_tce = times[vpt_all_idx[-1]].strftime('%Y-%m-%dT%H:%M:%SZ')

        # phase_sequence: slice original to kept sweeps only
        vpt_phase_seq = None
        if orig_phases is not None:
            vpt_phase_seq = ', '.join(orig_phases[i] for i in vpt_keep_indices)

        # Back up original then overwrite it
        bak_path = src_path.replace('.nc', '_original.nc')
        if not os.path.exists(bak_path):
            shutil.copy2(src_path, bak_path)
            print(f'  Backup: {bak_path}')
        else:
            print(f'  Backup already exists, skipping: {bak_path}')

        tmp_path = src_path + '.tmp'
        _write_cfradial(
            ds, tmp_path,
            vpt_all_idx, vpt_ss, vpt_se, vpt_fa,
            vpt_sm_arr, orig_scan_type, vpt_tcs, vpt_tce,
            ray_vars, range_vars, scalar_vars,
            history_msg=(
                f'Removed {n_trans} antenna_transition=1 ray(s) '
                f'({n_rhi_sweeps} group(s)) extracted to {rhi_base}; '
                f'{n_vpt_rays} ray(s), {n_vpt_sweeps} sweep(s) remain'
            ),
            phase_sequence=vpt_phase_seq,
        )
        ds.close()
        os.replace(tmp_path, src_path)
        print(f'  Written: {src_path}  '
              f'({n_vpt_rays} rays, {n_vpt_sweeps} sweeps)')
        print(f'\nDone.')
        return

    ds.close()
    print(f'\nDone.')


if __name__ == '__main__':
    main()
