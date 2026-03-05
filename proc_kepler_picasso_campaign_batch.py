#!/usr/bin/env python3
"""
proc_kepler_picasso_campaign_batch.py

Batch process PICASSO campaign data over a date range using campaign_processing module.
Processes RHI, PPI, VPT, and MAN (manual aircraft tracking) scans.

Usage:
    python proc_kepler_picasso_campaign_batch.py --start-date YYYYMMDD --end-date YYYYMMDD
    python proc_kepler_picasso_campaign_batch.py -d YYYYMMDD

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

def setup_picasso_paths(outpath=None, data_version='1.0.0'):
    """Set up file and directory paths for PICASSO campaign."""
    
    # Base PICASSO campaign paths
    base_inpath = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-mobile-ka-band-radar-1/data/campaign/picasso/mom'
    
    # Default output path if not specified
    if outpath is None:
        outpath = f'/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/picasso/L1_v{data_version}'
    
    paths = {
        'inpath': base_inpath,
        'outpath': outpath,
        'yaml_project_file': str(script_dir / 'campaigns' / 'picasso_project.yml'),
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
    try:
        start_date = datetime.datetime.strptime(start_date_str, '%Y%m%d')
        end_date = datetime.datetime.strptime(end_date_str, '%Y%m%d')
    except ValueError as e:
        raise ValueError(f"Invalid date format. Use YYYYMMDD. Error: {e}")
    
    if start_date > end_date:
        raise ValueError("Start date must be before or equal to end date")
    
    date_list = []
    current_date = start_date
    
    while current_date <= end_date:
        date_list.append(current_date.strftime('%Y%m%d'))
        current_date += datetime.timedelta(days=1)
    
    return date_list

def check_input_data(inpath, datestr):
    """
    Check if input data exists for the specified date.
    
    Args:
        inpath: Input directory path
        datestr: Date string in YYYYMMDD format
        
    Returns:
        bool: True if data exists, False otherwise
    """
    date_path = Path(inpath) / datestr
    
    if not date_path.exists():
        print(f"Date directory does not exist: {date_path}")
        return False
    
    # Check for mmclx files (both compressed and uncompressed)
    mmclx_files = list(date_path.rglob('*.mmclx*'))
    
    return len(mmclx_files) > 0

def main():
    """Main batch processing function."""
    parser = argparse.ArgumentParser(
        description="Batch process PICASSO campaign radar data over a date range",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--start-date', required=False, help='Start date in YYYYMMDD format')
    parser.add_argument('--end-date', required=False, help='End date in YYYYMMDD format')
    parser.add_argument('-d', '--date', help='Process a single date in YYYYMMDD format')
    parser.add_argument('--no-gzip', action='store_false', dest='gzip', 
                        help='Input files are NOT gzip compressed (default: files ARE gzipped)')
    parser.add_argument('--data-version', default='1.0.0', help='Data version string')
    parser.add_argument('--force', action='store_true',
                        help='Force processing even if output files exist')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be processed without actually processing')
    parser.add_argument('--skip-missing', action='store_true',
                        help='Skip dates with no input data instead of stopping')
    parser.add_argument('--single-sweep', action='store_true', default=True,
                        help='Create separate files for each sweep (default: True)')
    parser.add_argument('--split-man-phases', action='store_true', default=True,
                        help='Split MAN (manual tracking) scans into separate sweeps per phase (upward/dwell/downward/turning) (default: True)')
    parser.add_argument('--north-angle', type=float, default=None,
                        help='North angle correction (default: load from YAML config, 97.1 degrees for PICASSO)')
    parser.add_argument('--outpath', type=str, default=None,
                        help='Output directory path (auto-generated based on data-version if not specified)')
    
    args = parser.parse_args()
    
    print(f"Arguments: {args}")
    
    # Set up paths
    try:
        paths = setup_picasso_paths(outpath=args.outpath, data_version=args.data_version)
        print(f"PICASSO Campaign Processing")
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
    
    # Check for data availability
    print("\nChecking data availability:")
    dates_with_data = []
    dates_without_data = []
    
    for datestr in date_list:
        if check_input_data(paths['inpath'], datestr):
            dates_with_data.append(datestr)
            print(f"  {datestr} - HAS DATA")
        else:
            dates_without_data.append(datestr)
            print(f"  {datestr} - NO DATA")
    
    print(f"\nSummary: {len(dates_with_data)} dates with data, {len(dates_without_data)} without")
    
    if args.dry_run:
        print("\nDry run mode - no processing performed")
        print("\nConfiguration:")
        print(f"  Data version: {args.data_version}")
        print(f"  North angle: {args.north_angle}°")
        print(f"  Single sweep mode: {args.single_sweep}")
        print(f"  Split MAN phases: {args.split_man_phases}")
        print(f"  Gzip input: {args.gzip}")
        sys.exit(0)
    
    # Process each date
    success_count = 0
    error_count = 0
    skipped_count = 0
    
    for datestr in dates_with_data:
        try:
            print(f"\n{'='*60}")
            print(f"Processing {datestr}")
            print(f"{'='*60}")
            
            # Call the campaign processor
            process_campaign_day(
                campaign='picasso',
                datestr=datestr,
                inpath=paths['inpath'],
                outpath=paths['outpath'],
                yaml_project_file=paths['yaml_project_file'],
                yaml_instrument_file=paths['yaml_instrument_file'],
                gzip_flag=args.gzip,
                data_version=args.data_version,
                single_sweep=args.single_sweep,
                revised_northangle=args.north_angle,
                split_man_phases=args.split_man_phases
            )
            
            success_count += 1
            print(f"Successfully processed {datestr}")
            
        except Exception as e:
            error_count += 1
            print(f"Error processing {datestr}: {e}")
            if not args.skip_missing:
                print("Stopping due to error (use --skip-missing to continue)")
                sys.exit(1)
    
    # Process summary
    print(f"\n{'='*60}")
    print("PROCESSING SUMMARY")
    print(f"{'='*60}")
    print(f"Total dates: {len(date_list)}")
    print(f"  With data: {len(dates_with_data)}")
    print(f"  Without data: {len(dates_without_data)} (skipped)")
    print(f"  Processed successfully: {success_count}")
    print(f"  Errors: {error_count}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
