# COALESC3 Quicklooks Generation

This directory contains scripts for generating quicklook plots from COALESC3 campaign radar data.

## Campaign Information

- **Campaign**: COALESC3 (Combined Observations of the Atmospheric Boundary Layer to study the Evolution of StratoCumulus)
- **Location**: Weybourne Atmospheric Observatory (WAO), Norfolk, UK (52.95025°N, 1.12170°E)
- **Dates**: March - July 2017 (20170308 - 20170705)
- **Instrument**: NCAS Mobile Ka-band Radar (Kepler)
- **PI**: Simon Osborne (Met Office)

## Scripts

### 1. `make_coalesc3_quicklooks.py`

Python script to generate quicklook plots for a single day of COALESC3 data.

**Generates**:
- VPT (Vertical Pointing Time) plots showing DBZ, VEL, WIDTH, LDR
- VAD (Velocity Azimuth Display) plots showing wind direction and speed profiles

**Usage**:
```bash
# Process a single date
python make_coalesc3_quicklooks.py -d 20170410

# Specify custom paths
python make_coalesc3_quicklooks.py -d 20170410 \
    -i /path/to/input \
    -o /path/to/output

# Boundary layer mode (limits plots to 4km height)
python make_coalesc3_quicklooks.py -d 20170410 -b
```

**Arguments**:
- `-d, --date`: Date string in YYYYMMDD format (required)
- `-i, --inpath`: Input directory containing CF-Radial files (optional)
- `-o, --outpath`: Output directory for quicklook images (optional)
- `-b`: Boundary layer flag - limits plots to 4km height (optional)

### 2. `kepler_coalesc3_quicklooks.sh`

SLURM batch script for processing a single date on JASMIN.

**Usage**:
```bash
# Process a specific date
sbatch kepler_coalesc3_quicklooks.sh 20170410

# Or set as environment variable
DATESTR=20170410 sbatch kepler_coalesc3_quicklooks.sh

# Override paths
INPATH=/custom/input OUTPATH=/custom/output sbatch kepler_coalesc3_quicklooks.sh 20170410
```

### 3. `kepler_coalesc3_quicklooks_batch.sh`

SLURM array job script for batch processing multiple dates.

**Usage**:
```bash
# Process entire campaign (120 days from 2017-03-08 to 2017-07-05)
sbatch --array=1-120 kepler_coalesc3_quicklooks_batch.sh

# Process specific date range (e.g., March only - days 1-24)
sbatch --array=1-24 kepler_coalesc3_quicklooks_batch.sh

# Process with custom date range
START_DATE=20170401 END_DATE=20170430 sbatch --array=1-30 kepler_coalesc3_quicklooks_batch.sh

# Process non-contiguous dates (e.g., days 1, 5, 10-20)
sbatch --array=1,5,10-20 kepler_coalesc3_quicklooks_batch.sh
```

**Environment Variables**:
- `START_DATE`: Campaign start date (default: 20170308)
- `END_DATE`: Campaign end date (default: 20170705)
- `INPATH`: Input data path (default: `/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/coalesc3/L1_v1.0.1`)
- `OUTPATH`: Output quicklooks path (default: `$INPATH/quicklooks`)

## Output Structure

Quicklook plots are saved in subdirectories:

```
quicklooks/
├── vpt/
│   └── ncas-mobile-ka-band-radar-1_wao_coalesc3_YYYYMMDD_vpt_l1_VERSION.png
└── vad/
    └── ncas-mobile-ka-band-radar-1_wao_coalesc3_YYYYMMDD_vad_l1_VERSION.png
```

## Requirements

- Python 3.11+ (cao_3_11 conda environment)
- PyART
- NetCDF4
- NumPy
- Matplotlib
- cftime
- cmocean

## Examples

### Process a single test date
```bash
python make_coalesc3_quicklooks.py -d 20170410
```

### Process all March 2017 dates (24 days)
```bash
sbatch --array=1-24 kepler_coalesc3_quicklooks_batch.sh
```

### Process entire campaign
```bash
# Campaign duration: 2017-03-08 to 2017-07-05 = 120 days
sbatch --array=1-120 kepler_coalesc3_quicklooks_batch.sh
```

### Process with boundary layer focus
```bash
# Edit make_coalesc3_quicklooks.py call in the shell script to add -b flag
python make_coalesc3_quicklooks.py -d 20170410 -b
```

## Notes

- Processing one day typically takes ~5-30 minutes depending on data volume
- Memory requirement: ~64GB for standard processing
- SLURM logs are saved to `slurm_logs/` directory
- If a date has no processed data files, the script will skip it gracefully
- VPT plots include data from the last sweep of the previous day for continuity
