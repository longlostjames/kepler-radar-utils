#!/usr/bin/env python3
"""
proc_kepler_cobalt_campaign_batch.py

Batch process COBALT campaign data over a date range using campaign_processing module.

Usage:
    python proc_kepler_cobalt_campaign_batch.py --start-date YYYYMMDD --end-date YYYYMMDD

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

def setup_paths():
    """Set up file and directory paths for COBALT campaign."""
    home_path = Path.home()
    
    paths = {
        'inpath': '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-mobile-ka-band-radar-1/data/campaign/20241210_cobalt/mom',
        'outpath': '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/cobalt/L1c',
        'yaml_project_file': '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-mobile-ka-band-radar-1/data/campaign/20241210_cobalt/yaml/cobalt_project.yml',
        'yaml_instrument_file': '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-mobile-ka-band-radar-1/data/campaign/20241210_cobalt/yaml/amof_instruments.yml'
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
    
    print(f"Looking for data in: {date_path}")
    
    if not date_path.exists():
        print(f"Date directory does not exist: {date_path}")
        return False
    
    # Check for mmclx files (both compressed and uncompressed)
    mmclx_files = list(date_path.glob('*.mmclx*'))
    print(f"Found {len(mmclx_files)} mmclx files in {date_path}")
    
    if mmclx_files:
        print(f"Sample files: {[f.name for f in mmclx_files[:3]]}")
    else:
        # Debug: show what files are actually there
        all_files = list(date_path.glob('*'))
        print(f"All files in directory: {[f.name for f in all_files[:10]]}")
    
    return len(mmclx_files) > 0

def main():
    """Main batch processing function."""
    parser = argparse.ArgumentParser(
        description="Batch process COBALT campaign radar data over a date range",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--start-date', 
        required=False,  # Change this from True to False
        help='Start date in YYYYMMDD format'
    )
    
    parser.add_argument(
        '--end-date', 
        required=False,  # Change this from True to False
        help='End date in YYYYMMDD format'
    )
    
    parser.add_argument(
        '--gzip', 
        action='store_true',
        help='Input files are gzip compressed'
    )
    
    parser.add_argument(
        '--data-version', 
        default='1.0.1',
        help='Data version string'
    )
    
    parser.add_argument(
        '--force', 
        action='store_true',
        help='Force processing even if output files exist'
    )
    
    parser.add_argument(
        '--dry-run', 
        action='store_true',
        help='Show what would be processed without actually processing'
    )
    
    parser.add_argument(
        '--skip-missing', 
        action='store_true',
        help='Skip dates with no input data instead of stopping'
    )
    
    parser.add_argument(
        '-d', '--date',
        help='Process a single date in YYYYMMDD format (alternative to --start-date/--end-date)'
    )

    parser.add_argument(
        '--single-sweep', 
        action='store_true',
        help='Create separate files for each sweep instead of multi-sweep files'
    )

    parser.add_argument(
        '--no-vpt', 
        action='store_true',
        help='Disable vertical profiling for the processing'
    )

    args = parser.parse_args()
    
    # Set up paths
    try:
        paths = setup_paths()
        print(f"Input path: {paths['inpath']}")
        print(f"Output path: {paths['outpath']}")
    except Exception as e:
        print(f"Error setting up paths: {e}")
        sys.exit(1)
    
    # Generate date list
    try:
        if args.date:
            # Single date processing
            date_list = [args.date]
            # Validate the date format
            datetime.datetime.strptime(args.date, '%Y%m%d')
            print(f"Processing single date: {args.date}")
        elif args.start_date and args.end_date:
            # Date range processing
            date_list = generate_date_list(args.start_date, args.end_date)
            print(f"Processing {len(date_list)} dates from {args.start_date} to {args.end_date}")
        else:
            parser.error("Either --date (-d) or both --start-date and --end-date must be provided")
        
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Get campaign configuration
    campaign_info = get_campaign_info('cobalt')
    print(f"Campaign info: {campaign_info}")
    
    if args.dry_run:
        print("DRY RUN: Would process the following dates:")
        for datestr in date_list:
            has_data = check_input_data(paths['inpath'], datestr)
            status = "HAS DATA" if has_data else "NO DATA"
            print(f"  {datestr} - {status}")
        sys.exit(0)
    
    # Process each date
    total_processed = 0
    total_skipped = 0
    total_errors = 0
    
    overall_start_time = datetime.datetime.now()
    
    for i, datestr in enumerate(date_list, 1):
        print(f"\n{'='*60}")
        print(f"Processing date {i}/{len(date_list)}: {datestr}")
        print(f"{'='*60}")
        
        # Check for input data
        if not check_input_data(paths['inpath'], datestr):
            print(f"No input data found for {datestr}")
            if args.skip_missing:
                print("Skipping (--skip-missing enabled)")
                total_skipped += 1
                continue
            else:
                print("Stopping (use --skip-missing to skip missing dates)")
                break
        
        # Check if output already exists (unless forced)
        output_dir = Path(paths['outpath']) / datestr
        if output_dir.exists() and not args.force:
            existing_files = list(output_dir.glob('*.nc'))
            if existing_files:
                print(f"Output files already exist for {datestr}. Skipping (use --force to overwrite).")
                total_skipped += 1
                continue
        
        # Process the date
        try:
            print(f"Starting COBALT processing for {datestr}...")
            start_time = datetime.datetime.now()
            
            process_campaign_day(
                campaign='cobalt',
                datestr=datestr,
                inpath=paths['inpath'],
                outpath=paths['outpath'],
                yaml_project_file=str(paths['yaml_project_file']),
                yaml_instrument_file=str(paths['yaml_instrument_file']),
                gzip_flag=args.gzip,
                data_version=args.data_version,
                single_sweep=args.single_sweep,  # Add comma here
                no_vpt=args.no_vpt  # Pass the new parameter
            )
            
            end_time = datetime.datetime.now()
            duration = end_time - start_time
            
            print(f"Successfully completed COBALT processing for {datestr}")
            print(f"Processing time: {duration}")
            total_processed += 1
            
        except Exception as e:
            print(f"Error during processing of {datestr}: {e}")
            import traceback
            traceback.print_exc()
            total_errors += 1
            
            if not args.skip_missing:
                print("Stopping due to error (use --skip-missing to continue)")
                break
    
    # Summary
    overall_end_time = datetime.datetime.now()
    total_duration = overall_end_time - overall_start_time
    
    print(f"\n{'='*60}")
    print("BATCH PROCESSING SUMMARY")
    print(f"{'='*60}")
    print(f"Date range: {args.start_date} to {args.end_date}")
    print(f"Total dates in range: {len(date_list)}")
    print(f"Successfully processed: {total_processed}")
    print(f"Skipped: {total_skipped}")
    print(f"Errors: {total_errors}")
    print(f"Total processing time: {total_duration}")
    print(f"{'='*60}")
    
    # Exit with error code if there were any errors
    if total_errors > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()