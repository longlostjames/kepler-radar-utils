#!/usr/bin/env python3
"""
proc_kepler_reading_general_campaign_batch.py

Batch process Reading General campaign data over a date range using campaign_processing module.

Raw data is expected under /mnt/keplerdata/mom_backup/YYYY/YYYYMMDD/

Usage:
    python proc_kepler_reading_general_campaign_batch.py -d YYYYMMDD
    python proc_kepler_reading_general_campaign_batch.py --start-date YYYYMMDD --end-date YYYYMMDD
    python proc_kepler_reading_general_campaign_batch.py -d YYYYMMDD --latest

Author: Chris Walden, Science and Technology Facilities Council (STFC)
        as part of UK Research and Innovation (UKRI)
"""

import sys
import os
import re
import shutil
import argparse
import datetime
from pathlib import Path

# Add the script directory to Python path
script_dir = Path(__file__).parent.absolute()
sys.path.insert(0, str(script_dir))

from campaign_processing import process_campaign_day, get_campaign_info


def fix_midnight_crossings_for_date(outpath: str, datestr: str,
                                    force: bool = False) -> None:
    """Check VPT output files for midnight crossings and fix any that are found.

    Also checks the previous day's output directory for orphaned standalone
    spillover files (written when the next-day directory did not yet exist) and
    merges them into the current day's VPT file.

    For each VPT NetCDF file in <outpath>/<datestr>/ the function:
      1. Detects whether rays spill past midnight into the following day.
      2. Trims the source file to its own calendar day (backup saved as *_original.nc).
      3. Prepends the spillover as a new sweep 0 in the next-day file (or writes a
         standalone spillover file if no next-day file exists yet).

    Args:
        outpath: Base output directory (e.g. /data/processing/kepler/reading-general)
        datestr: Date string in YYYYMMDD format
        force:   Overwrite existing backup files if True
    """
    import netCDF4 as nc4
    import fix_midnight_crossing_vpt as fmcv

    date_outdir = Path(outpath) / datestr
    if not date_outdir.exists():
        print(f"  Midnight fix: output directory not found: {date_outdir}")
        return

    # ------------------------------------------------------------------
    # Step 0: merge any orphaned standalone spillover from the previous day.
    # Primary source: prev-day output directory.
    # Fallback source: spillover/ backup directory (used on reruns when the
    # standalone file in prev-day/ has already been consumed).
    # ------------------------------------------------------------------
    prev_datestr = (datetime.datetime.strptime(datestr, '%Y%m%d')
                    - datetime.timedelta(days=1)).strftime('%Y%m%d')
    prev_outdir = Path(outpath) / prev_datestr
    spillover_bak_dir = Path(outpath) / 'spillover'

    orphan_spillovers = []
    orphan_source = None

    if prev_outdir.exists():
        orphan_spillovers = sorted([
            f for f in prev_outdir.glob(f'*{datestr}*vpt*.nc')
            if '_original' not in f.name
        ])
        if orphan_spillovers:
            orphan_source = 'live'

    # Fallback: check the spillover backup directory when the live file is gone.
    if not orphan_spillovers and spillover_bak_dir.exists():
        orphan_spillovers = sorted([
            f for f in spillover_bak_dir.glob(f'*{datestr}*vpt*.nc')
            if '_original' not in f.name
        ])
        if orphan_spillovers:
            orphan_source = 'backup'
            print(f"  Midnight fix: no live spillover in {prev_datestr}/ — "
                  f"using backup from spillover/")

    if orphan_spillovers:
        today_vpt_files = sorted([
            f for f in date_outdir.glob('*vpt*.nc')
            if '_original' not in f.name
        ])
        if today_vpt_files:
            print(f"  Midnight fix: found {len(orphan_spillovers)} orphaned spillover "
                  f"file(s) from {prev_datestr} — merging into today's VPT...")
            for spillover_file in orphan_spillovers:
                # Match by filename similarity (strip the date token) or
                # fall back to the first (and usually only) today's VPT.
                spill_stem = re.sub(r'\d{8}-\d{6}|\d{8}', '', spillover_file.name)
                matched = next(
                    (f for f in today_vpt_files
                     if re.sub(r'\d{8}-\d{6}|\d{8}', '', f.name) == spill_stem),
                    today_vpt_files[0]
                )
                print(f"    {spillover_file.name} → {matched.name}")
                spill_date = datetime.datetime.strptime(prev_datestr, '%Y%m%d').date()
                if fmcv._spillover_already_present(str(matched), spill_date):
                    print(f"    Already merged — skipping")
                    continue
                if orphan_source == 'live':
                    # Back up before consuming the live file.
                    spillover_bak_dir.mkdir(exist_ok=True)
                    shutil.copy2(spillover_file, spillover_bak_dir / spillover_file.name)
                    print(f"    Backed up spillover to: spillover/{spillover_file.name}")
                    fmcv.merge_spillover_into_day(str(spillover_file), str(matched),
                                                  force=force)
                else:
                    # Sourced from backup: merge from a temporary copy so that
                    # merge_spillover_into_day's os.unlink() removes the copy,
                    # not the backup itself — keeping it intact for future reruns.
                    tmp_spill = spillover_file.with_suffix('.mergetmp.nc')
                    shutil.copy2(spillover_file, tmp_spill)
                    try:
                        fmcv.merge_spillover_into_day(str(tmp_spill), str(matched),
                                                      force=force)
                    except Exception:
                        tmp_spill.unlink(missing_ok=True)
                        raise
        else:
            print(f"  Midnight fix: orphaned spillover(s) from {prev_datestr} found "
                  f"but no today VPT to merge into yet — will retry next run.")

    # ------------------------------------------------------------------
    # Step 1: check today's VPT files for midnight crossings
    # ------------------------------------------------------------------
    vpt_files = sorted([
        f for f in date_outdir.glob('*vpt*.nc')
        if '_original' not in f.name
    ])
    if not vpt_files:
        print(f"  Midnight fix: no VPT files found in {date_outdir}")
        return

    print(f"  Midnight fix: checking {len(vpt_files)} VPT file(s)...")
    for vpt_file in vpt_files:
        src_path = str(vpt_file)
        print(f"    {vpt_file.name} ...", end=' ', flush=True)
        try:
            with nc4.Dataset(src_path, 'r') as ds:
                split = fmcv._find_split(ds.variables['time'])
        except ValueError:
            print("no crossing.")
            continue

        print(f"crossing at ray {split}.")
        src_date = fmcv._src_date(src_path)
        next_day_file = fmcv._detect_next_day_file(src_path, src_date)

        with nc4.Dataset(src_path, 'r') as ds:
            t = ds.variables['time']
            times = nc4.num2date(t[split:split + 1], t.units,
                                 only_use_cftime_datetimes=False)
        first_spill_time = datetime.datetime(
            times[0].year, times[0].month, times[0].day,
            times[0].hour, times[0].minute, times[0].second
        )
        out_path = fmcv._output_path_for_spillover(src_path, first_spill_time, next_day_file)
        bak = src_path.replace('.nc', '_original.nc')

        fmcv.trim_source(src_path, split, force=force)

        if next_day_file and fmcv._spillover_already_present(next_day_file, src_date):
            print(f"    Next-day file already starts with rays from {src_date} — "
                  f"skipping prepend (reprocessed day?)")
        else:
            # Always back up the spillover rays to spillover/ so reruns can
            # find them even after the live standalone file has been consumed.
            # Write a standalone spillover file directly into spillover/.
            spillover_bak_dir = Path(outpath) / 'spillover'
            spillover_bak_dir.mkdir(exist_ok=True)
            spill_bak_path = str(spillover_bak_dir / Path(out_path).name)
            if not Path(spill_bak_path).exists() or force:
                fmcv.prepend_spillover(bak, split, None, spill_bak_path,
                                       no_backup=True, force=force)
                print(f"    Backed up spillover to: spillover/{Path(out_path).name}")

            fmcv.prepend_spillover(bak, split, next_day_file, out_path, force=force)
            if next_day_file:
                print(f"    Spillover prepended to: {Path(out_path).name}")
            else:
                print(f"    Spillover written standalone: {Path(out_path).name}")
                print(f"    (will be merged into next-day VPT on the next processing run)")

        # Remove any stale next-day VPT files that have the right date token but
        # a different (earlier) start-time token than out_path.  These are left
        # over from a run that predated the rename-on-prepend logic.
        if out_path and Path(out_path).exists():
            out_dir = Path(out_path).parent
            date_token = Path(out_path).name  # e.g. '...20260403-000000...'
            # Extract the YYYYMMDD portion from out_path to match siblings
            _m = re.search(r'(\d{8})-\d{6}', Path(out_path).name)
            if _m:
                _day = _m.group(1)
                for _stale in out_dir.glob(f'*{_day}*vpt*.nc'):
                    if '_original' in _stale.name:
                        continue
                    if _stale == Path(out_path):
                        continue
                    _stale.unlink()
                    print(f"    Removed stale file: {_stale.name}")

    # ------------------------------------------------------------------
    # Step 2: move any *_original.nc files from the date directory to
    #         the central original/ directory.
    # ------------------------------------------------------------------
    original_dir = Path(outpath) / 'original'
    original_dir.mkdir(exist_ok=True)
    for orig_file in sorted(date_outdir.glob('*_original.nc')):
        dest = original_dir / orig_file.name
        if dest.exists() and not force:
            print(f"  Original backup already in original/: {orig_file.name} — skipping move")
        else:
            shutil.move(str(orig_file), dest)
            print(f"  Moved to original/: {orig_file.name}")


def setup_reading_general_paths(use_latest=True):
    """Set up file and directory paths for Reading General campaign."""

    # Base input path - raw data is under YYYY/YYYYMMDD below this
    base_inpath = '/mnt/keplerdata/mom_backup'

    if use_latest:
        print("Using 'latest' subdirectory for input data")

    paths = {
        'inpath': base_inpath,
        'outpath': '/data/processing/kepler/reading-general',
        'yaml_project_file': str(script_dir / 'campaigns' / 'reading-general_project.yml'),
        'yaml_instrument_file': str(script_dir / 'instrument_metadata.yml'),
        'use_latest': use_latest
    }

    # Ensure output directory exists
    Path(paths['outpath']).mkdir(parents=True, exist_ok=True)

    return paths


def get_date_input_path(base_inpath, datestr, use_latest=False):
    """
    Get the input path for the given date.

    Raw data is stored as <base_inpath>/YYYY/YYYYMMDD/
    with an optional 'latest' subdirectory.

    Args:
        base_inpath: Base input directory path
        datestr: Date string in YYYYMMDD format
        use_latest: If True, use datestr/latest subdirectory

    Returns:
        str: Full input path for the date
    """
    year = datestr[:4]
    date_path = os.path.join(base_inpath, year, datestr)
    if use_latest:
        return os.path.join(date_path, 'latest')
    return date_path


def check_input_data(base_inpath, datestr, use_latest=False):
    """
    Check if input data exists for the specified date.

    Args:
        base_inpath: Base input directory path
        datestr: Date string in YYYYMMDD format
        use_latest: If True, look in the 'latest' subdirectory

    Returns:
        bool: True if mmclx data files exist, False otherwise
    """
    date_path = Path(get_date_input_path(base_inpath, datestr, use_latest))
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
        all_files = list(date_path.glob('*'))
        print(f"All files in directory: {[f.name for f in all_files[:10]]}")

    return len(mmclx_files) > 0


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


def main():
    """Main batch processing function."""
    parser = argparse.ArgumentParser(
        description="Batch process Reading General radar data over a date range",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--start-date',
        required=False,
        help='Start date in YYYYMMDD format'
    )

    parser.add_argument(
        '--end-date',
        required=False,
        help='End date in YYYYMMDD format'
    )

    parser.add_argument(
        '-d', '--date',
        help='Process a single date in YYYYMMDD format (alternative to --start-date/--end-date)'
    )

    parser.add_argument(
        '--latest',
        action='store_true',
        help='Process files from the "latest" subdirectory within each date directory'
    )

    parser.add_argument(
        '--gzip',
        action='store_true',
        help='Input files are gzip compressed'
    )

    parser.add_argument(
        '--data-version',
        default='1.0.0',
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
        '--single-sweep',
        action='store_true',
        help='Create separate files for each sweep instead of multi-sweep files'
    )

    parser.add_argument(
        '--north-angle',
        type=float,
        default=None,
        help='North angle correction in degrees (default: load from YAML config)'
    )

    parser.add_argument(
        '--no-vpt',
        action='store_true',
        help='Disable vertical profiling for the processing'
    )

    parser.add_argument(
        '--max-age',
        type=int,
        default=6,
        help='Maximum age of input files in hours (default: 6 hours)'
    )

    parser.add_argument(
        '--no-midnight-fix',
        action='store_true',
        help='Disable automatic midnight crossing fix for VPT files after processing'
    )

    args = parser.parse_args()

    print(f"Arguments: {args}")

    # Set up paths
    try:
        paths = setup_reading_general_paths(use_latest=args.latest)
        print("Reading General Campaign Processing")
        if args.latest:
            print("Mode: Processing from 'latest' subdirectories")
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
                print("You may need to create this YAML file or update the path")

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
        campaign_info = get_campaign_info('reading-general')
        print(f"Campaign info: {campaign_info}")
    except Exception as e:
        print(f"Warning: Could not get campaign info for reading-general: {e}")
        campaign_info = {
            'tracking_tag': 'reading-general',
            'location': 'cao',
            'revised_northangle': args.north_angle
        }

    if args.dry_run:
        print("\nDRY RUN: Would process the following dates:")
        for datestr in date_list:
            has_data = check_input_data(paths['inpath'], datestr, use_latest=args.latest)
            status = "HAS DATA" if has_data else "NO DATA"
            input_path = get_date_input_path(paths['inpath'], datestr, use_latest=args.latest)
            print(f"  {datestr} - {status} (from: {input_path})")
        print(f"\nConfiguration:")
        print(f"  Data version: {args.data_version}")
        print(f"  Single sweep mode: {args.single_sweep}")
        print(f"  North angle correction: {args.north_angle}°")
        print(f"  Gzip input: {args.gzip}")
        print(f"  Use latest subdirectory: {args.latest}")
        sys.exit(0)

    # Process each date
    total_processed = 0
    total_skipped = 0
    total_errors = 0

    overall_start_time = datetime.datetime.now()

    for i, datestr in enumerate(date_list, 1):
        print(f"\n{'='*60}")
        print(f"Processing Reading General date {i}/{len(date_list)}: {datestr}")
        if args.latest:
            print(f"Using 'latest' subdirectory for {datestr}")
        print(f"{'='*60}")

        # Check for input data
        if not check_input_data(paths['inpath'], datestr, use_latest=args.latest):
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
            print(f"Starting Reading General processing for {datestr}...")
            start_time = datetime.datetime.now()

            date_input_path = get_date_input_path(paths['inpath'], datestr, use_latest=args.latest)
            print(f"Processing from: {date_input_path}")

            process_campaign_day(
                campaign='reading-general',
                datestr=datestr,
                inpath=date_input_path,
                outpath=paths['outpath'],
                yaml_project_file=str(paths['yaml_project_file']),
                yaml_instrument_file=str(paths['yaml_instrument_file']),
                gzip_flag=args.gzip,
                data_version=args.data_version,
                single_sweep=args.single_sweep,
                revised_northangle=args.north_angle,
                no_vpt=args.no_vpt,
                max_age=args.max_age,
                logp_dir='/mnt/keplerdata/Logp',
                force=args.force
            )

            if not args.no_midnight_fix:
                try:
                    fix_midnight_crossings_for_date(
                        paths['outpath'], datestr, force=args.force
                    )
                except Exception as mc_exc:
                    print(f"Warning: midnight crossing fix failed for {datestr}: {mc_exc}")

            end_time = datetime.datetime.now()
            duration = end_time - start_time

            print(f"Successfully completed Reading General processing for {datestr}")
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
    print("Reading General Processing Summary")
    print(f"{'='*60}")
    print(f"Total dates processed: {total_processed}")
    print(f"Total dates skipped:   {total_skipped}")
    print(f"Total errors:          {total_errors}")
    print(f"Total time:            {total_duration}")
    print(f"{'='*60}")

    sys.exit(0 if total_errors == 0 else 1)


if __name__ == '__main__':
    main()
