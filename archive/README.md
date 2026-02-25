# Archive Directory

This directory contains deprecated and old code that has been moved from the main repository to keep it organized.

## Directory Structure

### `old_scripts/`
Deprecated Python scripts that have been superseded by `campaign_processing.py` and the updated `kepler_utils.py`:

- **`kepler_utils_old.py`** - Old version of utility functions with hardcoded metadata
- **`proc_kepler2ncas_cobalt.py`** - Old COBALT processing script (replaced by `proc_kepler_cobalt_campaign_batch.py`)
- **`proc_kepler2ncas_cobalt_rhi.py`** - Old COBALT RHI processing script
- **`proc_kepler2ncas_kasbex.py`** - Old KASBEX processing script (replaced by `proc_kepler_kasbex_campaign_batch.py`)
- **`proc_kepler_cobalt_campaign_batch_old.py`** - Earlier version of campaign batch processing
- **`add_processing_software_metadata.py`** - Utility to retroactively add metadata (no longer needed as new processing includes it)
- **`instrument_metadata.yml`** - Old metadata file (superseded by campaign-specific YAML files in `campaigns/`)

### `old_bash/`
Bash scripts that used the old Python scripts:

- **`kepler_cobalt_batch*.sh`** - Old batch processing scripts for COBALT
- **`kepler_cobalt_cron*.sh`** - Old cron job scripts  
- **`kepler_proc_cobalt*.sh`** - Old processing launcher scripts
- **`kepler_cobalt_lotus2_[1-5].sh`** - Old LOTUS2 array job scripts
- **`kepler_quicklooks_cobalt[0-8].sh`** - Old numbered quicklook scripts

### `backups/`
Timestamped backup files (`.backup.YYYYMMDD_HHMMSS`) created during development:

- Various backup copies of quicklook scripts

### `logs/`
Old log files from previous processing runs:

- `clearout.log` - Log from file cleanup operations
- `kasbex_cron.log` - KASBEX cron job logs
- `proc.log` - General processing logs
- `quicklook.log` - Quicklook generation logs

## Current Active Code

The main repository now uses:

- **`campaign_processing.py`** - Unified campaign processing module
- **`kepler_utils.py`** - Updated utility functions with YAML-based metadata
- **`proc_kepler_*_campaign_batch.py`** - Modern campaign-specific processing scripts
- **Campaign YAML files** in `campaigns/` directory for metadata management

## Notes

- Editor backup files (`*~`, `#*#`) have been deleted as they were temporary files
- These archived files are kept for reference and potential recovery
- Do not use archived scripts for new processing - use the current campaign processing system

---
*Archived: February 25, 2026*
