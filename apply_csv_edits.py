#!/usr/bin/env python3
"""
Apply edited CSV files to NetCDF files for a given date.

The script discovers which NetCDF files to update by scanning the CSV directory
for files ending in '_sweeps.csv' or '_antenna_transition.csv'. The corresponding
NetCDF filename is derived by stripping that suffix and appending '.nc'. This
means any scan type (vpt, man, ppi, ...) can be edited.

Usage:
    python apply_csv_edits.py -d YYYYMMDD -i /path/to/data [-c csv_dir] [--no-backup]

Example:
    python apply_csv_edits.py -d 20171216 -i /gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/picasso/L1_v1.0.0
    python apply_csv_edits.py -d 20250305 -i /path/to/data -c /path/to/csv_edits
    python apply_csv_edits.py -d 20250305 -i /path/to/data --no-backup
"""

import argparse
import sys
from kepler_utils import apply_csv_edits_for_date


def parse_command_line():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Apply edited CSV files to NetCDF files for a given date'
    )
    
    parser.add_argument(
        '-d', '--date',
        required=True,
        help='Date string in YYYYMMDD format'
    )
    
    parser.add_argument(
        '-i', '--input-dir',
        required=True,
        help='Base directory containing date subdirectories with NetCDF files'
    )
    
    parser.add_argument(
        '-c', '--csv-dir',
        default=None,
        help='Directory containing CSV files in YYYYMMDD subdirectories (default: same as input-dir)'
    )
    
    parser.add_argument(
        '--no-backup',
        action='store_true',
        help='Do not create backup files'
    )
    
    args = parser.parse_args()
    
    # Validate date format
    import datetime
    try:
        datetime.datetime.strptime(args.date, '%Y%m%d')
    except ValueError:
        print(f"ERROR: Invalid date format '{args.date}'. Use YYYYMMDD format.")
        sys.exit(1)
    
    return args


def main():
    """Main entry point."""
    args = parse_command_line()
    
    print("="*70)
    print("Apply CSV Edits to NetCDF Files")
    print("="*70)
    print(f"Date:        {args.date}")
    print(f"Data dir:    {args.input_dir}")
    print(f"CSV dir:     {args.csv_dir if args.csv_dir else args.input_dir + ' (same as data dir)'}")
    print(f"Backup:      {'No' if args.no_backup else 'Yes'}")
    print("="*70)
    print()
    
    try:
        apply_csv_edits_for_date(
            datestr=args.date,
            data_dir=args.input_dir,
            csv_dir=args.csv_dir,
            backup=not args.no_backup,
        )
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print("\nDone!")


if __name__ == '__main__':
    main()
