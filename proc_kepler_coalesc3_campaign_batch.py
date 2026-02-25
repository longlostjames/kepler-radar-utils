#!/usr/bin/env python3
"""
proc_kepler_coalesc3_campaign_batch.py

Batch process COALESC3 campaign data over a date range using campaign_processing module.

Usage:
    python proc_kepler_coalesc3_campaign_batch.py --start-date YYYYMMDD --end-date YYYYMMDD
    python proc_kepler_coalesc3_campaign_batch.py -d YYYYMMDD

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
"""

import sys
import os
import argparse
import datetime
from pathlib import Path

# Add the script directory to Python path
script_dir = Path(__file__).parent.absolute()
sys.path.insert(0, str(script_dir))

from campaign_processing import process_campaign_day, get_campaign_info

def setup_coalesc3_paths(outpath=None):
    """Set up file and directory paths for COALESC3 campaign."""

    # Base COALESC3 campaign paths
    base_inpath = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-mobile-ka-band-radar-1/data/campaign/coalesc3/mom'

    # Default output path if not specified
    if outpath is None:
        outpath = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/coalesc3/L1_v1.0.1'

    paths = {
        'inpath': base_inpath,
        'outpath': outpath,
        'yaml_project_file': str(script_dir / 'campaigns' / 'coalesc3_project.yml'),
        'yaml_instrument_file': str(script_dir / 'instrument_metadata.yml')
    }

    # Ensure output directory exists
    Path(paths['outpath']).mkdir(parents=True, exist_ok=True)

    return paths

def generate_date_list(start_date_str, end_date_str):
    """
    Generate a list of date strings between start and end dates (inclusive).

    Args:
        start_date_str: Start date in YYYYMMDD format
        end_date_str: End date in YYYYMMDD format

    Returns:
        List of date strings in YYYYMMDD format
    """
    start_date = datetime.datetime.strptime(start_date_str, '%Y%m%d')
    end_date = datetime.datetime.strptime(end_date_str, '%Y%m%d')

    date_list = []
    current_date = start_date
    while current_date <= end_date:
        date_list.append(current_date.strftime('%Y%m%d'))
        current_date += datetime.timedelta(days=1)

    return date_list

def check_input_data(inpath, datestr):
    """
    Check if input data exists for the given date.

    Args:
        inpath: Base input directory path
        datestr: Date string in YYYYMMDD format

    Returns:
        bool: True if data exists, False otherwise
    """
    date_path = Path(inpath) / datestr
    print(f"Looking for data in: {date_path}")

    if not date_path.exists():
        print(f"Date directory does not exist: {date_path}")
        return False

    mmclx_files = list(date_path.glob('*.mmclx')) + list(date_path.glob('*.mmclx.gz'))
    if not mmclx_files:
        print(f"No mmclx files found in: {date_path}")
        return False

    print(f"Found {len(mmclx_files)} mmclx files")
    return True

def main():
    """Main batch processing function."""
    parser = argparse.ArgumentParser(
        description="Batch process COALESC3 campaign radar data over a date range",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--start-date', required=False, help='Start date in YYYYMMDD format')
    parser.add_argument('--end-date', required=False, help='End date in YYYYMMDD format')
    parser.add_argument('-d', '--date', help='Process a single date in YYYYMMDD format')
    parser.add_argument('--gzip', action='store_true', default=True, help='Input files are gzip compressed (default: True)')
    parser.add_argument('--data-version', default='1.0.1', help='Data version string')
    parser.add_argument('--force', action='store_true',
                        help='Force processing even if output files exist')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be processed without actually processing')
    parser.add_argument('--skip-missing', action='store_true',
                        help='Skip dates with no input data instead of stopping')
    parser.add_argument('--single-sweep', action='store_true',
                        help='Create separate files for each sweep')
    parser.add_argument('--north-angle', type=float, default=287.0,
                        help='North angle correction for COALESC3 deployment (default: 287.0 degrees)')
    parser.add_argument('--outpath', type=str, default=None,
                        help='Output directory path (default: /gws/.../coalesc3/L1_v1.0.1)')
    parser.add_argument('--no-vpt', action='store_true', help='Disable vertical profiling')

    args = parser.parse_args()

    print(f"Arguments: {args}")

    # Set up paths
    try:
        paths = setup_coalesc3_paths(outpath=args.outpath)
        print(f"COALESC3 Campaign Processing")
        print(f"Base input path: {paths['inpath']}")
        print(f"Output path: {paths['outpath']}")
        print(f"Project YAML: {paths['yaml_project_file']}")
        print(f"Instrument YAML: {paths['yaml_instrument_file']}")
    except Exception as e:
        print(f"Error setting up paths: {e}")
        sys.exit(1)

    # Validate that required files exist
    for key, path in paths.items():
        if key.endswith('_file'):
            if not Path(path).exists():
                print(f"Warning: {key} not found at {path}")

    # Generate date list
    try:
        if args.date:
            date_list = [args.date]
            datetime.datetime.strptime(args.date, '%Y%m%d')
            print(f"Processing single date: {args.date}")
        elif args.start_date and args.end_date:
            date_list = generate_date_list(args.start_date, args.end_date)
            print(f"Processing {len(date_list)} dates from {args.start_date} to {args.end_date}")
        else:
            parser.error("Either --date (-d) or both --start-date and --end-date must be provided")
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Get campaign configuration
    try:
        campaign_info = get_campaign_info('coalesc3')
        print(f"Campaign info: {campaign_info}")
    except Exception as e:
        print(f"Warning: Could not get campaign info for COALESC3: {e}")
        campaign_info = {
            'tracking_tag': 'AMF_07092016101810',
            'location': 'wao',
            'revised_northangle': args.north_angle
        }

    if args.dry_run:
        print("\nDRY RUN: Would process the following dates:")
        for datestr in date_list:
            has_data = check_input_data(paths['inpath'], datestr)
            status = "HAS DATA" if has_data else "NO DATA"
            print(f"  {datestr} - {status}")
        print(f"\nConfiguration:")
        print(f"  Data version: {args.data_version}")
        print(f"  Single sweep mode: {args.single_sweep}")
        print(f"  North angle correction: {args.north_angle}°")
        print(f"  Gzip input: {args.gzip}")
        sys.exit(0)

    # Process each date
    total_processed = 0
    total_skipped = 0
    total_errors = 0

    overall_start_time = datetime.datetime.now()

    for i, datestr in enumerate(date_list, 1):
        print(f"\n{'='*60}")
        print(f"Processing COALESC3 date {i}/{len(date_list)}: {datestr}")
        print(f"{'='*60}")

        if not check_input_data(paths['inpath'], datestr):
            print(f"No input data found for {datestr}")
            if args.skip_missing:
                print("Skipping (--skip-missing enabled)")
                total_skipped += 1
                continue
            else:
                print("Stopping (use --skip-missing to skip missing dates)")
                break

        output_dir = Path(paths['outpath']) / datestr
        if output_dir.exists() and not args.force:
            existing_files = list(output_dir.glob('*.nc'))
            if existing_files:
                print(f"Output files already exist for {datestr}. Skipping (use --force to overwrite).")
                total_skipped += 1
                continue

        try:
            print(f"Starting COALESC3 processing for {datestr}...")
            start_time = datetime.datetime.now()

            process_campaign_day(
                campaign='coalesc3',
                datestr=datestr,
                inpath=paths['inpath'],
                outpath=paths['outpath'],
                yaml_project_file=str(paths['yaml_project_file']),
                yaml_instrument_file=str(paths['yaml_instrument_file']),
                gzip_flag=args.gzip,
                data_version=args.data_version,
                single_sweep=args.single_sweep,
                revised_northangle=args.north_angle,
                no_vpt=args.no_vpt
            )

            duration = datetime.datetime.now() - start_time
            print(f"✓ Successfully completed COALESC3 processing for {datestr}")
            print(f"Processing time: {duration}")
            total_processed += 1

        except Exception as e:
            print(f"✗ Error during processing of {datestr}: {e}")
            import traceback
            traceback.print_exc()
            total_errors += 1
            if not args.skip_missing:
                break

    overall_duration = datetime.datetime.now() - overall_start_time
    print(f"\n{'='*60}")
    print(f"COALESC3 BATCH PROCESSING SUMMARY")
    print(f"{'='*60}")
    print(f"Total processed: {total_processed}")
    print(f"Total skipped:   {total_skipped}")
    print(f"Total errors:    {total_errors}")
    print(f"Total time:      {overall_duration}")

if __name__ == '__main__':
    main()
