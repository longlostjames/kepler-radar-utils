#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=kasbex-processing
#SBATCH -o slurm_logs/kasbex_%A_%a.out
#SBATCH -e slurm_logs/kasbex_%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# KASBEX Campaign Processing - SLURM Array Job
# Usage: sbatch --array=0-N submit_kasbex_processing.sh
# where N = number of days - 1

# Configuration - THESE WILL BE UPDATED BY THE HELPER SCRIPT
START_DATE="20250802"  # YYYYMMDD format
END_DATE="20250809"    # YYYYMMDD format

# Processing options
GZIP_FLAG="--gzip"           # Set to "--gzip" if input files are compressed
DATA_VERSION="1.0.0"   # Data version string
SINGLE_SWEEP_FLAG="--single-sweep"   # Set to "--single-sweep" for individual sweep files
NORTH_ANGLE="55.9"     # North angle correction (degrees)
FORCE_FLAG=""          # Set to "--force" to overwrite existing files
SKIP_MISSING_FLAG="--skip-missing"  # Skip dates with no data

# ============================================================================
# Script execution - DO NOT MODIFY BELOW THIS LINE
# ============================================================================

# Load conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Calculate the target date based on array task ID
start_epoch=$(date -d "$START_DATE" +%s)
target_epoch=$((start_epoch + SLURM_ARRAY_TASK_ID * 86400))  # 86400 = seconds per day
target_date=$(date -d "@$target_epoch" +%Y%m%d)

# Validate date range
end_epoch=$(date -d "$END_DATE" +%s)
if [ $target_epoch -gt $end_epoch ]; then
    echo "Error: Array task ID $SLURM_ARRAY_TASK_ID results in date $target_date beyond end date $END_DATE"
    exit 1
fi

# Build command line arguments
cmd_args="-d $target_date"

if [ -n "$GZIP_FLAG" ]; then
    cmd_args="$cmd_args $GZIP_FLAG"
fi

if [ -n "$DATA_VERSION" ]; then
    cmd_args="$cmd_args --data-version $DATA_VERSION"
fi

if [ -n "$SINGLE_SWEEP_FLAG" ]; then
    cmd_args="$cmd_args $SINGLE_SWEEP_FLAG"
fi

if [ -n "$NORTH_ANGLE" ]; then
    cmd_args="$cmd_args --north-angle $NORTH_ANGLE"
fi

if [ -n "$FORCE_FLAG" ]; then
    cmd_args="$cmd_args $FORCE_FLAG"
fi

if [ -n "$SKIP_MISSING_FLAG" ]; then
    cmd_args="$cmd_args $SKIP_MISSING_FLAG"
fi

# Print job information
echo "============================================================"
echo "KASBEX Campaign Processing"
echo "============================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start Date: $START_DATE"
echo "End Date: $END_DATE"
echo "Processing Date: $target_date"
echo "Command: proc_kepler_kasbex_campaign_batch.py $cmd_args"
echo "Configuration:"
echo "  Data version: $DATA_VERSION"
echo "  North angle: ${NORTH_ANGLE}°"
echo "  Gzip input: $([ -n "$GZIP_FLAG" ] && echo "Yes" || echo "No")"
echo "  Single sweep: $([ -n "$SINGLE_SWEEP_FLAG" ] && echo "Yes" || echo "No")"
echo "  Force overwrite: $([ -n "$FORCE_FLAG" ] && echo "Yes" || echo "No")"
echo "============================================================"

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Create log directory if it doesn't exist
mkdir -p slurm_logs

# Execute the KASBEX processing script
time python $SCRIPT_DIR/proc_kepler_kasbex_campaign_batch.py $cmd_args

exit_code=$?

echo "============================================================"
echo "Job completed with exit code: $exit_code"
echo "Processing date: $target_date"
echo "End time: $(date)"
echo "============================================================"

exit $exit_code
