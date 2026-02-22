#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=cobalt-quicklooks
#SBATCH -o slurm_logs/quicklooks_%A_%a.out
#SBATCH -e slurm_logs/quicklooks_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# COBALT Quicklooks Generation - SLURM Array Job
# This script will be updated by submit_cobalt_quicklooks.sh

# Configuration - THESE WILL BE UPDATED BY THE HELPER SCRIPT
START_DATE="20241210"  # YYYYMMDD format
END_DATE="20251201"    # YYYYMMDD format

# Optional flags
BOUNDARY_LAYER_FLAG=""  # Set to "-b" to limit plots to 4km height
LATEST_FLAG=""          # Set to "-l" for latest data processing

# Input/Output paths (optional - leave empty to use defaults)
INPUT_PATH=""   # Override default input path if needed
OUTPUT_PATH=""  # Override default output path if needed

# ============================================================================
# Script execution - DO NOT MODIFY BELOW THIS LINE
# ============================================================================

# Load conda environment
source $HOME/anaconda3/etc/profile.d/conda.sh
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

if [ -n "$BOUNDARY_LAYER_FLAG" ]; then
    cmd_args="$cmd_args $BOUNDARY_LAYER_FLAG"
fi

if [ -n "$LATEST_FLAG" ]; then
    cmd_args="$cmd_args $LATEST_FLAG"
fi

if [ -n "$INPUT_PATH" ]; then
    cmd_args="$cmd_args -i $INPUT_PATH"
fi

if [ -n "$OUTPUT_PATH" ]; then
    cmd_args="$cmd_args -o $OUTPUT_PATH"
fi

# Print job information
echo "============================================================"
echo "COBALT Quicklooks Generation"
echo "============================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start Date: $START_DATE"
echo "End Date: $END_DATE"
echo "Processing Date: $target_date"
echo "Command: make_kepler_cobalt_quicklooks.py $cmd_args"
echo "============================================================"

# Create log directory if it doesn't exist
mkdir -p slurm_logs

# Execute the quicklooks script
time python /home/users/cjwalden/git/kepler-radar-utils-cobalt/make_cobalt_quicklooks.py $cmd_args

exit_code=$?

echo "============================================================"
echo "Job completed with exit code: $exit_code"
echo "Processing time: $(date)"
echo "============================================================"

exit $exit_code
