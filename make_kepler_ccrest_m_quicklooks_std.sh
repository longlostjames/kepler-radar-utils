#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=ccrest-m-quicklooks
#SBATCH -o slurm_logs/quicklooks_%A_%a.out
#SBATCH -e slurm_logs/quicklooks_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=128G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# CCREST-M Quicklooks Generation - SLURM Array Job
# This script will be updated by submit_kepler_ccrest_m_quicklooks.sh

# Configuration - values can be overridden via sbatch --export
START_DATE="20240219"  # YYYYMMDD format
END_DATE="20240219"      # YYYYMMDD format

# Optional flags
BOUNDARY_LAYER_FLAG="${BOUNDARY_LAYER_FLAG:-}"       # Set to -b to limit plots to 4 km height
SKIP_ALL_TRANSITION_FLAG="${SKIP_ALL_TRANSITION_FLAG:-}"  # Set to --skip-all-transition to skip 100% transition sweeps
VPT_ONLY_FLAG="${VPT_ONLY_FLAG:-}"             # Set to --vpt-only to generate only VPT plots
PPI_ONLY_FLAG="${PPI_ONLY_FLAG:-}"             # Set to --ppi-only to generate only PPI/PPI-map plots
POINTING_ONLY_FLAG="${POINTING_ONLY_FLAG:-}"      # Set to --pointing-only to generate only pointing plots
PPI_MAP_DAY_EXTENT_FLAG="${PPI_MAP_DAY_EXTENT_FLAG:-}"  # Set to --ppi-map-day-extent for day-wide PPI-map axis limits
PPI_MAP_COLORBAR_SHRINK="${PPI_MAP_COLORBAR_SHRINK:-}"  # Set numeric value (e.g. 0.82) for PPI-map colorbar size

# Input/Output paths (optional - leave empty to use defaults)
INPUT_PATH="${INPUT_PATH:-}"   # Override default input path if needed
OUTPUT_PATH="${OUTPUT_PATH:-}"  # Override default output path if needed

# ============================================================================
# Script execution - DO NOT MODIFY BELOW THIS LINE
# ============================================================================

# Load conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Calculate the target date based on array task ID
start_epoch=$(date -d "$START_DATE" +%s)
target_epoch=$((start_epoch + SLURM_ARRAY_TASK_ID * 86400))
target_date=$(date -d "@$target_epoch" +%Y%m%d)

# Validate date range
end_epoch=$(date -d "$END_DATE" +%s)
if [ $target_epoch -gt $end_epoch ]; then
    echo "Error: Array task ID $SLURM_ARRAY_TASK_ID results in date $target_date beyond end date $END_DATE"
    exit 1
fi

# Build command line arguments
cmd_args="-d $target_date"

if [ -n "$BOUNDARY_LAYER_FLAG" ]; then
    cmd_args="$cmd_args $BOUNDARY_LAYER_FLAG"
fi

if [ -n "$SKIP_ALL_TRANSITION_FLAG" ]; then
    cmd_args="$cmd_args $SKIP_ALL_TRANSITION_FLAG"
fi

if [ -n "$VPT_ONLY_FLAG" ]; then
    cmd_args="$cmd_args $VPT_ONLY_FLAG"
fi

if [ -n "$PPI_ONLY_FLAG" ]; then
    cmd_args="$cmd_args $PPI_ONLY_FLAG"
fi

if [ -n "$POINTING_ONLY_FLAG" ]; then
    cmd_args="$cmd_args $POINTING_ONLY_FLAG"
fi

if [ -n "$PPI_MAP_DAY_EXTENT_FLAG" ]; then
    cmd_args="$cmd_args $PPI_MAP_DAY_EXTENT_FLAG"
fi

if [ -n "$PPI_MAP_COLORBAR_SHRINK" ]; then
    cmd_args="$cmd_args --ppi-map-colorbar-shrink $PPI_MAP_COLORBAR_SHRINK"
fi

if [ -n "$INPUT_PATH" ]; then
    cmd_args="$cmd_args -i $INPUT_PATH"
fi

if [ -n "$OUTPUT_PATH" ]; then
    cmd_args="$cmd_args -o $OUTPUT_PATH"
fi

# Print job information
echo "============================================================"
echo "CCREST-M Quicklooks Generation"
echo "============================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start Date: $START_DATE"
echo "End Date: $END_DATE"
echo "Processing Date: $target_date"
if [ -n "$VPT_ONLY_FLAG" ]; then
    echo "Plot mode: VPT only"
fi
if [ -n "$PPI_ONLY_FLAG" ]; then
    echo "Plot mode: PPI + PPI-map only"
fi
if [ -n "$POINTING_ONLY_FLAG" ]; then
    echo "Plot mode: Pointing only"
fi
if [ -n "$PPI_MAP_DAY_EXTENT_FLAG" ]; then
    echo "PPI-map extent mode: day outermost bounds"
fi
if [ -n "$PPI_MAP_COLORBAR_SHRINK" ]; then
    echo "PPI-map colorbar shrink: $PPI_MAP_COLORBAR_SHRINK"
fi
echo "Command: make_ccrest_quicklooks.py $cmd_args"
echo "============================================================"

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Create log directory if it doesn't exist
mkdir -p slurm_logs

# Execute the quicklooks script
time python $SCRIPT_DIR/make_ccrest_quicklooks.py $cmd_args

exit_code=$?

echo "============================================================"
echo "Job completed with exit code: $exit_code"
echo "Processing time: $(date)"
echo "============================================================"

exit $exit_code
