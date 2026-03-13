#!/usr/bin/env python3
"""
Export sweep metadata and antenna_transition flags from NetCDF files to CSV for a given date.

Usage:
    python export_sweep_metadata.py -d YYYYMMDD -i /path/to/data [-p pattern] [--overwrite]

Example:
    python export_sweep_metadata.py -d 20171213 -i /gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/picasso/L1_v1.0.0
    python export_sweep_metadata.py -d 20171213 -i /path/to/data -p '*_man_*.nc'
    python export_sweep_metadata.py -d 20171213 -i /path/to/data --overwrite
"""

import argparse
import sys
from kepler_utils import export_sweep_metadata_for_date


def parse_command_line():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Export sweep metadata and antenna_transition flags from NetCDF files to CSV'
    )
    
    parser.add_argument(
        '-d', '--date',
        required=True,
        help='Date string in YYYYMMDD format'
    )
    
    parser.add_argument(
        '-i', '--input-dir',
        required=True,
        help='Base directory containing date subdirectories'
    )
    
    parser.add_argument(
        '-p', '--pattern',
        default='*.nc',
        help='Filename pattern to match (default: *.nc)'
    )
    
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite existing CSV files'
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
    print("Export Sweep Metadata from NetCDF Files")
    print("="*70)
    print(f"Date:      {args.date}")
    print(f"Directory: {args.input_dir}")
    print(f"Pattern:   {args.pattern}")
    print(f"Overwrite: {'Yes' if args.overwrite else 'No'}")
    print("="*70)
    print()
    
    try:
        export_sweep_metadata_for_date(
            datestr=args.date,
            data_dir=args.input_dir,
            pattern=args.pattern,
            overwrite=args.overwrite
        )
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print("\nDone!")


if __name__ == '__main__':
    main()
