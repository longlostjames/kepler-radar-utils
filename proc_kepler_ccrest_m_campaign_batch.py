#!/usr/bin/env python3
"""
proc_kepler_ccrest_m_campaign_batch.py

Batch process CCREST-M campaign data over a date range using campaign_processing module.

Usage:
    python proc_kepler_ccrest_m_campaign_batch.py --start-date YYYYMMDD --end-date YYYYMMDD
    python proc_kepler_ccrest_m_campaign_batch.py -d YYYYMMDD
    python proc_kepler_ccrest_m_campaign_batch.py -d YYYYMMDD --latest

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

def setup_ccrest_m_paths(use_latest=False, outpath=None):
    """Set up file and directory paths for CCREST-M campaign."""

    # Base CCREST-M campaign paths
    base_inpath = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-mobile-ka-band-radar-1/data/campaign/ccrest-m/mom'

    if use_latest:
        print("Using 'latest' subdirectory for input data")

    # Default output path if not specified
    if outpath is None:
        outpath = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/ccrest-m/L1_v1.0.0'

    paths = {
        'inpath': base_inpath,
        'outpath': outpath,
        'yaml_project_file': str(script_dir / 'campaigns' / 'ccrest-m_project.yml'),
        'yaml_instrument_file': str(script_dir / 'instrument_metadata.yml'),
        'use_latest': use_latest
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

def check_input_data(inpath, datestr, use_latest=False):
    """
    Check if input data exists for the given date.

    Args:
        inpath: Base input directory path
        datestr: Date string in YYYYMMDD format
        use_latest: If True, look in datestr/latest subdirectory

    Returns:
        bool: True if data exists, False otherwise
    """
    if use_latest:
        date_path = Path(inpath) / datestr / 'latest'
        print(f"Looking for data in LATEST subdirectory: {date_path}")
    else:
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

def get_input_path(base_inpath, datestr, use_latest=False):
    """Get the appropriate input path for the given date."""
    if use_latest:
        return os.path.join(base_inpath, datestr, 'latest')
    else:
        return os.path.join(base_inpath, datestr)

def main():
    """Main batch processing function."""
    parser = argparse.ArgumentParser(
        description="Batch process CCREST-M campaign radar data over a date range",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--start-date', required=False, help='Start date in YYYYMMDD format')
    parser.add_argument('--end-date', required=False, help='End date in YYYYMMDD format')
    parser.add_argument('-d', '--date', help='Process a single date in YYYYMMDD format')
    parser.add_argument('--latest', action='store_true',
                        help='Process files from the "latest" subdirectory')
    parser.add_argument('--gzip', action='store_true', help='Input files are gzip compressed')
    parser.add_argument('--data-version', default='1.0.0', help='Data version string')
    parser.add_argument('--force', action='store_true',
                        help='Force processing even if output files exist')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be processed without actually processing')
    parser.add_argument('--skip-missing', action='store_true',
                        help='Skip dates with no input data instead of stopping')
    parser.add_argument('--single-sweep', action='store_true',
                        help='Create separate files for each sweep')
    parser.add_argument('--north-angle', type=float, default=None,
                        help='North angle correction (default: load from YAML config, typically 55.7 degrees for CCREST-M)')
    parser.add_argument('--outpath', type=str, default=None,
                        help='Output directory path (default: /gws/.../ccrest-m/L1_v1.0.0)')
    parser.add_argument('--no-vpt', action='store_true', help='Disable vertical profiling')
    parser.add_argument('--max-age', type=int, default=6,
                        help='Maximum age of input files in hours (default: 6 hours)')
    parser.add_argument('--stationary-threshold', type=float, default=0.5,
                        help='Per-ray |d_az| threshold (deg) for stationary ray classification')
    parser.add_argument('--stationary-std-threshold', type=float, default=1.0,
                        help='Maximum azimuth std (deg) for a stationary segment')
    parser.add_argument('--stationary-range-threshold', type=float, default=2.0,
                        help='Maximum azimuth range (deg) for a stationary segment')
    parser.add_argument('--stationary-total-change-threshold', type=float, default=3.0,
                        help='Maximum total unwrapped azimuth excursion (deg) for stationary segments')
    parser.add_argument('--stationary-scan-rate-threshold', type=float, default=0.2,
                        help='Median |scan_rate| threshold (deg/s) above which a segment is treated as scanning')
    parser.add_argument('--direction-rate-threshold', type=float, default=0.1,
                        help='Median per-ray d_az threshold (deg/ray) for scan direction')
    parser.add_argument('--pointing-window-size', type=int, default=10,
                        help='Rolling window size (rays) for pointing split detection')
    parser.add_argument('--pointing-az-std-threshold', type=float, default=1.0,
                        help='Rolling-window azimuth std threshold (deg) for pointing split')
    parser.add_argument('--min-pointing-rays', type=int, default=10,
                        help='Minimum contiguous rays to keep a pointing phase')
    parser.add_argument('--pointing-total-change-threshold', type=float, default=3.0,
                        help='Maximum total unwrapped azimuth excursion (deg) for pointing phases')
    parser.add_argument('--pointing-scan-rate-threshold', type=float, default=0.2,
                        help='Median |scan_rate| threshold (deg/s) above which pointing phases are reclassified as ppi')

    args = parser.parse_args()

    print(f"Arguments: {args}")

    # Set up paths
    try:
        paths = setup_ccrest_m_paths(use_latest=args.latest, outpath=args.outpath)
        print(f"CCREST-M Campaign Processing")
        if args.latest:
            print(f"Mode: Processing from 'latest' subdirectories")
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
        campaign_info = get_campaign_info('ccrest-m')
        print(f"Campaign info: {campaign_info}")
    except Exception as e:
        print(f"Warning: Could not get campaign info for CCREST-M: {e}")
        campaign_info = {
            'tracking_tag': 'AMOF_20230201132601',
            'location': 'cao',
            'revised_northangle': args.north_angle
        }

    if args.dry_run:
        print("\nDRY RUN: Would process the following dates:")
        for datestr in date_list:
            has_data = check_input_data(paths['inpath'], datestr, use_latest=args.latest)
            status = "HAS DATA" if has_data else "NO DATA"
            input_path = get_input_path(paths['inpath'], datestr, use_latest=args.latest)
            print(f"  {datestr} - {status} (from: {input_path})")
        print(f"\nConfiguration:")
        print(f"  Data version: {args.data_version}")
        print(f"  Single sweep mode: {args.single_sweep}")
        print(f"  North angle correction: {args.north_angle}°")
        print(f"  Gzip input: {args.gzip}")
        print(f"  Use latest subdirectory: {args.latest}")
        print(f"  Stationary threshold: {args.stationary_threshold} deg/ray")
        print(f"  Stationary std threshold: {args.stationary_std_threshold} deg")
        print(f"  Stationary range threshold: {args.stationary_range_threshold} deg")
        print(f"  Stationary total change threshold: {args.stationary_total_change_threshold} deg")
        print(f"  Stationary scan rate threshold: {args.stationary_scan_rate_threshold} deg/s")
        print(f"  Direction rate threshold: {args.direction_rate_threshold} deg/ray")
        print(f"  Pointing window size: {args.pointing_window_size} rays")
        print(f"  Pointing az std threshold: {args.pointing_az_std_threshold} deg")
        print(f"  Minimum pointing rays: {args.min_pointing_rays}")
        print(f"  Pointing total change threshold: {args.pointing_total_change_threshold} deg")
        print(f"  Pointing scan rate threshold: {args.pointing_scan_rate_threshold} deg/s")
        sys.exit(0)

    # Process each date
    total_processed = 0
    total_skipped = 0
    total_errors = 0

    overall_start_time = datetime.datetime.now()

    for i, datestr in enumerate(date_list, 1):
        print(f"\n{'='*60}")
        print(f"Processing CCREST-M date {i}/{len(date_list)}: {datestr}")
        print(f"{'='*60}")

        if not check_input_data(paths['inpath'], datestr, use_latest=args.latest):
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
            print(f"Starting CCREST-M processing for {datestr}...")
            start_time = datetime.datetime.now()

            process_campaign_day(
                campaign='ccrest-m',
                datestr=datestr,
                inpath=paths['inpath'],
                outpath=paths['outpath'],
                yaml_project_file=str(paths['yaml_project_file']),
                yaml_instrument_file=str(paths['yaml_instrument_file']),
                gzip_flag=args.gzip,
                data_version=args.data_version,
                single_sweep=args.single_sweep,
                revised_northangle=args.north_angle,
                no_vpt=args.no_vpt,
                max_age=args.max_age,
                stationary_threshold=args.stationary_threshold,
                stationary_std_threshold=args.stationary_std_threshold,
                stationary_range_threshold=args.stationary_range_threshold,
                stationary_total_change_threshold=args.stationary_total_change_threshold,
                stationary_scan_rate_threshold=args.stationary_scan_rate_threshold,
                direction_rate_threshold=args.direction_rate_threshold,
                pointing_window_size=args.pointing_window_size,
                pointing_az_std_threshold=args.pointing_az_std_threshold,
                min_pointing_rays=args.min_pointing_rays,
                pointing_total_change_threshold=args.pointing_total_change_threshold,
                pointing_scan_rate_threshold=args.pointing_scan_rate_threshold,
            )

            duration = datetime.datetime.now() - start_time
            print(f"✓ Successfully completed CCREST-M processing for {datestr}")
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
    print(f"CCREST-M BATCH PROCESSING SUMMARY")
    print(f"{'='*60}")
    print(f"Total processed: {total_processed}")
    print(f"Total skipped:   {total_skipped}")
    print(f"Total errors:    {total_errors}")
    print(f"Total time:      {overall_duration}")

if __name__ == '__main__':
    main()
