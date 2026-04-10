#!/usr/bin/env python3
"""
remove_antenna_transition_rays_vpt.py
--------------------------------------
Remove rays flagged as antenna_transition=1 from CF-Radial VPT L1 files.

For each input file:
  - All rays with antenna_transition==1 are removed from every ray-indexed variable.
  - Sweep start/end ray indices are recalculated for the retained rays.
  - Sweeps that become entirely empty (all their rays were in transition) are dropped.
  - time_coverage_start, time_coverage_end, history, and last_revised_date are
    updated to reflect the retained data.
  - The original file is backed up as <file>.bak before being replaced
    (unless --no-backup is given).

Usage:
    python remove_antenna_transition_rays_vpt.py [options] <file> [<file> ...]

Options:
    --dry-run    Print what would be done but do not modify any files.
    --no-backup  Do not create .bak backups before overwriting.

Examples:
    python remove_antenna_transition_rays_vpt.py \\
        /path/to/L1/20241216/ncas-mobile-ka-band-radar-1_cao_20241216-000007_vpt_l1_v1.0.3.nc

    python remove_antenna_transition_rays_vpt.py --dry-run \\
        /path/to/L1/20241216/*.nc
"""

import argparse
import datetime
import os
import shutil
import socket
import sys

import netCDF4 as nc4
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _history_entry(msg: str) -> str:
    """Return a history line: 'Www Mon DD HH:MM:SS YYYY - user:<u> machine:<h> <msg>'"""
    now = datetime.datetime.utcnow()
    user = os.environ.get('USER', 'unknown')
    host = socket.gethostname()
    ts = now.strftime('%a %b %d %H:%M:%S %Y')
    return f'{ts} - user:{user} machine:{host} {msg}'


def _copy_var_attrs(src_var, dst_var):
    for attr in src_var.ncattrs():
        if attr != '_FillValue':
            setattr(dst_var, attr, _nc_char(getattr(src_var, attr)))


def _nc_char(val):
    """Encode a str value as bytes so netCDF4 writes NC_CHAR not NC_STRING.
    Non-string values (int, float, arrays) are returned unchanged."""
    if isinstance(val, str):
        return val.encode('utf-8')
    return val


def _copy_global_attrs(src_ds, dst_ds, overrides=None):
    for attr in src_ds.ncattrs():
        val = overrides.get(attr, getattr(src_ds, attr)) if overrides else getattr(src_ds, attr)
        setattr(dst_ds, attr, _nc_char(val))
    if overrides:
        for k, v in overrides.items():
            if k not in src_ds.ncattrs():
                setattr(dst_ds, k, _nc_char(v))


def _str_to_char_array(s: str, length: int) -> np.ndarray:
    """Encode a string into a fixed-length char (S1) array."""
    arr = np.zeros(length, dtype='S1')
    b = s.encode('ascii')
    for i, c in enumerate(b[:length]):
        arr[i] = bytes([c])
    return arr


def _time_string(t_val: float, units: str) -> str:
    """Convert a scalar time value + units to an ISO 8601 string."""
    import cftime
    dt = cftime.num2pydate(np.array([t_val]), units)[0]
    return dt.strftime('%Y-%m-%dT%H:%M:%SZ')


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def remove_transition_rays(nc_file: str, backup: bool = True, dry_run: bool = False) -> bool:
    """
    Remove antenna_transition==1 rays from *nc_file* in-place.

    Returns True if the file was (or would be) modified, False if nothing
    needed to change (no transition rays found or file already clean).
    """
    print(f'\n{"="*70}')
    print(f'File: {os.path.basename(nc_file)}')

    with nc4.Dataset(nc_file, 'r') as ds:
        if 'antenna_transition' not in ds.variables:
            print('  No antenna_transition variable — skipping.')
            return False

        at = ds.variables['antenna_transition'][:]
        n_rays_orig = len(at)
        n_transition = int(np.sum(at == 1))

        if n_transition == 0:
            print('  No transition rays found — nothing to do.')
            return False

        print(f'  Original rays : {n_rays_orig}')
        print(f'  Transition    : {n_transition}  ({100*n_transition/n_rays_orig:.1f} %)')
        print(f'  Retained      : {n_rays_orig - n_transition}')

        # Build retained-ray index array
        keep_indices = np.where(at == 0)[0]   # old ray indices to keep (sorted)
        n_rays_new = len(keep_indices)

        # Build old→new ray index mapping (only for kept rays)
        old2new = np.full(n_rays_orig, -1, dtype=np.int64)
        old2new[keep_indices] = np.arange(n_rays_new, dtype=np.int64)

        # Rebuild sweep index arrays, dropping empty sweeps
        ss_orig = ds.variables['sweep_start_ray_index'][:]
        se_orig = ds.variables['sweep_end_ray_index'][:]
        nsweeps_orig = len(ss_orig)

        new_sweep_rows = []  # list of (new_ss, new_se, original_sweep_idx)
        for i in range(nsweeps_orig):
            s, e = int(ss_orig[i]), int(se_orig[i])
            # Which kept rays fall within this sweep?
            mask = (keep_indices >= s) & (keep_indices <= e)
            sweep_kept = keep_indices[mask]
            if len(sweep_kept) == 0:
                continue  # drop entirely-transitional sweep
            new_sweep_rows.append((
                int(old2new[sweep_kept[0]]),
                int(old2new[sweep_kept[-1]]),
                i
            ))

        nsweeps_new = len(new_sweep_rows)
        n_dropped_sweeps = nsweeps_orig - nsweeps_new
        print(f'  Original sweeps : {nsweeps_orig}')
        print(f'  Dropped sweeps  : {n_dropped_sweeps}')
        print(f'  Retained sweeps : {nsweeps_new}')

        if dry_run:
            print('  [DRY RUN] No changes written.')
            return True

        # Identify variable categories by their leading dimension
        ray_vars   = [v for v in ds.variables
                      if ds.variables[v].dimensions
                      and ds.variables[v].dimensions[0] == 'time']
        range_vars = [v for v in ds.variables
                      if ds.variables[v].dimensions
                      and ds.variables[v].dimensions[0] == 'range']
        sweep_vars = [v for v in ds.variables
                      if ds.variables[v].dimensions
                      and ds.variables[v].dimensions[0] == 'sweep'
                      and v not in ('sweep_start_ray_index', 'sweep_end_ray_index')]
        string_vars = [v for v in ds.variables
                       if ds.variables[v].dimensions
                       and ds.variables[v].dimensions[0] == 'string_length']
        scalar_vars = [v for v in ds.variables
                       if not ds.variables[v].dimensions]
        freq_vars   = [v for v in ds.variables
                       if ds.variables[v].dimensions
                       and ds.variables[v].dimensions[0] == 'frequency']

        # Derive updated time_coverage_start/end from retained rays
        time_var = ds.variables['time']
        time_units = time_var.units
        tcs = _time_string(float(time_var[keep_indices[0]]),  time_units)
        tce = _time_string(float(time_var[keep_indices[-1]]), time_units)

        string_length = ds.dimensions['string_length'].size

        # Collect all data to write while the source file is still open
        ray_data   = {}
        for v in ray_vars:
            src = ds.variables[v]
            ray_data[v] = src[keep_indices] if src.ndim == 1 else src[keep_indices, :]

        range_data = {}
        for v in range_vars:
            range_data[v] = ds.variables[v][:]

        sweep_kept_indices = [r[2] for r in new_sweep_rows]
        sweep_data = {}
        for v in sweep_vars:
            src = ds.variables[v]
            if src.ndim == 1:
                sweep_data[v] = src[sweep_kept_indices]
            else:
                sweep_data[v] = src[sweep_kept_indices, :]

        scalar_data = {v: ds.variables[v][:] for v in scalar_vars}
        freq_data   = {v: ds.variables[v][:] for v in freq_vars}
        string_data = {}
        for v in string_vars:
            if v == 'time_coverage_start':
                string_data[v] = _str_to_char_array(tcs, string_length)
            elif v == 'time_coverage_end':
                string_data[v] = _str_to_char_array(tce, string_length)
            else:
                string_data[v] = ds.variables[v][:]

        new_ss = np.array([r[0] for r in new_sweep_rows], dtype=ss_orig.dtype)
        new_se = np.array([r[1] for r in new_sweep_rows], dtype=se_orig.dtype)

        # Preserve variable metadata (dtype, fill_value, attributes)
        var_meta = {}
        for v in ds.variables:
            vv = ds.variables[v]
            fv = vv._FillValue if hasattr(vv, '_FillValue') else None
            attrs = {a: _nc_char(getattr(vv, a)) for a in vv.ncattrs() if a != '_FillValue'}
            var_meta[v] = {'dtype': vv.dtype, 'dims': vv.dimensions,
                           'fv': fv, 'attrs': attrs}

        # Global attributes (will override a few below)
        existing_history = getattr(ds, 'history', '')

    # -----------------------------------------------------------------------
    # Write temporary file, then atomically replace original
    # -----------------------------------------------------------------------
    tmp_file = nc_file + '.tmp'

    n_range = range_data[range_vars[0]].shape[0] if range_vars else 573

    with nc4.Dataset(tmp_file, 'w', format='NETCDF4') as out:
        # Dimensions
        out.createDimension('time',          None)          # unlimited
        out.createDimension('range',         n_range)
        out.createDimension('sweep',         nsweeps_new)
        out.createDimension('string_length', string_length)
        if freq_vars:
            out.createDimension('frequency', freq_data[freq_vars[0]].shape[0])

        # Global attributes
        history_entry = _history_entry(
            f'Removed {n_transition} antenna_transition=1 rays '
            f'({n_dropped_sweeps} empty sweeps also dropped). '
            f'Retained {n_rays_new} of {n_rays_orig} rays in {nsweeps_new} sweeps.'
        )
        new_history = (history_entry + '\n' + existing_history).strip()
        now_iso = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
        _copy_global_attrs(nc4.Dataset(nc_file, 'r'), out, overrides={
            'history': new_history,
            'last_revised_date': now_iso,
            'time_coverage_start': tcs,
            'time_coverage_end':   tce,
        })

        def _make_var(name, dims, compress=True):
            m = var_meta[name]
            kw = {'fill_value': m['fv']} if m['fv'] is not None else {}
            if compress and any(d in ('time', 'range', 'sweep') for d in dims):
                kw.update({'zlib': True, 'complevel': 4})
            v = out.createVariable(name, m['dtype'], dims, **kw)
            for a, val in m['attrs'].items():
                setattr(v, a, val)
            return v

        # Ray-indexed variables
        for vname in ray_vars:
            m = var_meta[vname]
            dst = _make_var(vname, m['dims'])
            dst[:] = ray_data[vname]

        # Range variable
        for vname in range_vars:
            m = var_meta[vname]
            dst = _make_var(vname, m['dims'])
            dst[:] = range_data[vname]

        # Sweep variables (other than start/end ray index)
        for vname in sweep_vars:
            m = var_meta[vname]
            dst = _make_var(vname, m['dims'])
            dst[:] = sweep_data[vname]

        # Sweep start/end ray index (rebuilt)
        for vname, arr in (('sweep_start_ray_index', new_ss),
                            ('sweep_end_ray_index',   new_se)):
            m = var_meta[vname]
            dst = _make_var(vname, m['dims'])
            dst[:] = arr

        # Scalar variables
        for vname in scalar_vars:
            m = var_meta[vname]
            dst = out.createVariable(vname, m['dtype'], m['dims'])
            for a, val in m['attrs'].items():
                setattr(dst, a, val)
            dst[:] = scalar_data[vname]

        # Frequency variables
        for vname in freq_vars:
            m = var_meta[vname]
            dst = _make_var(vname, m['dims'])
            dst[:] = freq_data[vname]

        # String (char-array) variables
        for vname in string_vars:
            m = var_meta[vname]
            dst = out.createVariable(vname, m['dtype'], m['dims'])
            for a, val in m['attrs'].items():
                setattr(dst, a, val)
            dst[:] = string_data[vname]

        # Ensure azimuth/elevation/range use proposed_standard_name, not standard_name.
        # These are NCAS proposed names, not CF standard names.
        for vname, psn in (('azimuth',   'ray_azimuth_angle'),
                           ('elevation', 'ray_elevation_angle'),
                           ('range',     'projection_range_coordinate')):
            if vname in out.variables:
                v = out.variables[vname]
                if hasattr(v, 'standard_name'):
                    v.delncattr('standard_name')
                v.setncattr('proposed_standard_name', _nc_char(psn))

    # Backup and replace
    if backup:
        bak = nc_file + '.bak'
        shutil.copy2(nc_file, bak)
        print(f'  Backup : {bak}')

    os.replace(tmp_file, nc_file)
    print(f'  Written: {nc_file}')
    print(f'  time_coverage: {tcs}  →  {tce}')
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Remove antenna_transition=1 rays from CF-Radial VPT L1 files.'
    )
    parser.add_argument('files', nargs='+', metavar='FILE',
                        help='One or more CF-Radial NetCDF files to process.')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print what would be done without modifying files.')
    parser.add_argument('--no-backup', action='store_true',
                        help='Do not create .bak backups.')
    return parser.parse_args()


def main():
    args = parse_args()
    n_modified = 0
    n_skipped  = 0

    for nc_file in args.files:
        if not os.path.isfile(nc_file):
            print(f'WARNING: File not found — {nc_file}')
            n_skipped += 1
            continue
        try:
            changed = remove_transition_rays(
                nc_file,
                backup=not args.no_backup,
                dry_run=args.dry_run,
            )
            if changed:
                n_modified += 1
            else:
                n_skipped += 1
        except Exception as exc:
            print(f'ERROR processing {nc_file}: {exc}')
            import traceback
            traceback.print_exc()
            n_skipped += 1

    print(f'\n{"="*70}')
    if args.dry_run:
        print(f'DRY RUN complete. {n_modified} file(s) would be modified, {n_skipped} skipped.')
    else:
        print(f'Done. {n_modified} file(s) modified, {n_skipped} skipped.')


if __name__ == '__main__':
    main()
