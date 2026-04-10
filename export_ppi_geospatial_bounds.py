#!/usr/bin/env python
"""
export_ppi_geospatial_bounds.py
--------------------------------
Extract geospatial bounds from all PPI L1 NetCDF files for a given day and
write them as a CSV.

Usage
-----
    python export_ppi_geospatial_bounds.py YYYYMMDD INPATH [--output OUTPUT.csv]

Arguments
---------
    YYYYMMDD   Date string, e.g. 20240209
    INPATH     Root input path; PPI files are expected under INPATH/YYYYMMDD/*ppi*.nc
    --output   Output CSV file path  (default: ppi_bounds_YYYYMMDD.csv)
"""
import argparse
import csv
import glob
import os
import re
import sys

try:
    import netCDF4 as nc4
except ImportError:
    sys.exit("Error: netCDF4 is required. Install it with: conda install netCDF4")


# ---------------------------------------------------------------------------
# Bounds parsing helpers
# ---------------------------------------------------------------------------

def _parse_geospatial_bounds(bounds_str):
    """Parse 'Bounding box: 50.95N -1.77E, 51.29N -1.44E' -> (lon_min, lon_max, lat_min, lat_max)."""
    m = re.search(
        r'([+-]?\d+\.?\d*)([NS])\s+([+-]?\d+\.?\d*)([EW])'
        r'\s*,\s*'
        r'([+-]?\d+\.?\d*)([NS])\s+([+-]?\d+\.?\d*)([EW])',
        bounds_str,
    )
    if m is None:
        return None
    lat1 = float(m.group(1)) * (-1 if m.group(2) == 'S' else 1)
    lon1 = float(m.group(3)) * (-1 if m.group(4) == 'W' else 1)
    lat2 = float(m.group(5)) * (-1 if m.group(6) == 'S' else 1)
    lon2 = float(m.group(7)) * (-1 if m.group(8) == 'W' else 1)
    return min(lon1, lon2), max(lon1, lon2), min(lat1, lat2), max(lat1, lat2)


def extract_bounds(ncfile):
    """Return (lon_min, lon_max, lat_min, lat_max, source) or None."""
    with nc4.Dataset(ncfile) as ds:
        raw = getattr(ds, 'geospatial_bounds', '')
        bounds = _parse_geospatial_bounds(raw)
        if bounds is not None:
            return bounds + ('geospatial_bounds',)

        lon_min = getattr(ds, 'geospatial_lon_min', None)
        lon_max = getattr(ds, 'geospatial_lon_max', None)
        lat_min = getattr(ds, 'geospatial_lat_min', None)
        lat_max = getattr(ds, 'geospatial_lat_max', None)
        if None not in (lon_min, lon_max, lat_min, lat_max):
            return float(lon_min), float(lon_max), float(lat_min), float(lat_max), 'geospatial_lon/lat_min/max'

    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Export geospatial bounds for all PPI files on a given day.'
    )
    parser.add_argument('datestr', help='Date string YYYYMMDD')
    parser.add_argument('inpath', help='Root input path (PPI files under INPATH/YYYYMMDD/)')
    parser.add_argument('--output', default=None, help='Output CSV file (default: ppi_bounds_YYYYMMDD.csv)')
    args = parser.parse_args()

    inpath_date = os.path.join(args.inpath, args.datestr)
    if not os.path.isdir(inpath_date):
        sys.exit(f"Error: directory not found: {inpath_date}")

    ppi_files = sorted(glob.glob(os.path.join(inpath_date, '*ppi*.nc')))
    if not ppi_files:
        sys.exit(f"No PPI files found in {inpath_date}")

    print(f"Found {len(ppi_files)} PPI file(s) in {inpath_date}")

    output_csv = args.output or f'ppi_bounds_{args.datestr}.csv'

    rows = []
    day_bounds = None

    for ncfile in ppi_files:
        basename = os.path.basename(ncfile)
        try:
            result = extract_bounds(ncfile)
        except Exception as exc:
            print(f"  WARNING: {basename} -> could not open ({exc})")
            rows.append({'filename': basename, 'lon_min': '', 'lon_max': '',
                         'lat_min': '', 'lat_max': '', 'source': f'error: {exc}'})
            continue

        if result is None:
            print(f"  {basename} -> no bounds found")
            rows.append({'filename': basename, 'lon_min': '', 'lon_max': '',
                         'lat_min': '', 'lat_max': '', 'source': 'unavailable'})
            continue

        lon_min, lon_max, lat_min, lat_max, source = result
        print(f"  {basename} ({source})")
        print(f"    lon [{lon_min:.4f}, {lon_max:.4f}]  lat [{lat_min:.4f}, {lat_max:.4f}]")

        rows.append({
            'filename': basename,
            'lon_min': f'{lon_min:.6f}',
            'lon_max': f'{lon_max:.6f}',
            'lat_min': f'{lat_min:.6f}',
            'lat_max': f'{lat_max:.6f}',
            'source': source,
        })

        # Accumulate day-wide union
        if day_bounds is None:
            day_bounds = [lon_min, lon_max, lat_min, lat_max]
        else:
            day_bounds[0] = min(day_bounds[0], lon_min)
            day_bounds[1] = max(day_bounds[1], lon_max)
            day_bounds[2] = min(day_bounds[2], lat_min)
            day_bounds[3] = max(day_bounds[3], lat_max)

    # Append a summary row with the day-wide union
    if day_bounds is not None:
        rows.append({
            'filename': f'DAY_UNION ({args.datestr})',
            'lon_min': f'{day_bounds[0]:.6f}',
            'lon_max': f'{day_bounds[1]:.6f}',
            'lat_min': f'{day_bounds[2]:.6f}',
            'lat_max': f'{day_bounds[3]:.6f}',
            'source': 'union',
        })
        print(f"\nDay-wide union:")
        print(f"  lon [{day_bounds[0]:.4f}, {day_bounds[1]:.4f}]  lat [{day_bounds[2]:.4f}, {day_bounds[3]:.4f}]")

    fieldnames = ['filename', 'lon_min', 'lon_max', 'lat_min', 'lat_max', 'source']
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWritten to {output_csv}")


if __name__ == '__main__':
    main()
