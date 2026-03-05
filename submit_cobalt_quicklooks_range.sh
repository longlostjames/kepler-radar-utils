#!/bin/bash
# Helper script to automatically calculate array range and submit

if [ $# -ne 2 ]; then
    echo "Usage: $0 START_DATE END_DATE"
    echo "Example: $0 20250110 20250131"
    exit 1
fi

START_DATE=$1
END_DATE=$2

# Name of the actual SLURM submit script
SUBMIT_SCRIPT="make_kepler_cobalt_quicklooks.sh"

# Check if the submit script exists
if [ ! -f "$SUBMIT_SCRIPT" ]; then
    echo "Error: Submit script '$SUBMIT_SCRIPT' not found in current directory"
    echo "Available files:"
    ls -la *.sh
    exit 1
fi

# Update the submit script with new dates
sed -i "s/START_DATE=\".*\"/START_DATE=\"$START_DATE\"/" "$SUBMIT_SCRIPT"
sed -i "s/END_DATE=\".*\"/END_DATE=\"$END_DATE\"/" "$SUBMIT_SCRIPT"

# Verify the sed commands worked
if ! grep -q "START_DATE=\"$START_DATE\"" "$SUBMIT_SCRIPT"; then
    echo "Warning: Failed to update START_DATE in $SUBMIT_SCRIPT"
    echo "Please check the file format and try again"
    exit 1
fi

if ! grep -q "END_DATE=\"$END_DATE\"" "$SUBMIT_SCRIPT"; then
    echo "Warning: Failed to update END_DATE in $SUBMIT_SCRIPT"
    echo "Please check the file format and try again"
    exit 1
fi

# Calculate number of days
start_epoch=$(date -d "$START_DATE" +%s)
end_epoch=$(date -d "$END_DATE" +%s)
num_days=$(( ( end_epoch - start_epoch ) / 86400 ))

if [ $num_days -lt 0 ]; then
    echo "Error: End date must be after start date"
    exit 1
fi

echo "Updated $SUBMIT_SCRIPT with:"
echo "  START_DATE=$START_DATE"
echo "  END_DATE=$END_DATE"
echo ""
echo "Submitting quicklooks job for $START_DATE to $END_DATE ($((num_days + 1)) days)"
echo "Array range: 0-$num_days"

# Submit the job
sbatch --array=0-$num_days "$SUBMIT_SCRIPT"

if [ $? -eq 0 ]; then
    echo "Job submitted successfully!"
else
    echo "Error: Job submission failed"
    exit 1
fi
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
# Usage: sbatch --array=0-N make_kepler_cobalt_quicklooks.sh
# where N = number of days - 1

# Configuration - THESE WILL BE UPDATED BY THE HELPER SCRIPT
START_DATE="20241210"  # YYYYMMDD format
END_DATE="20250728"    # YYYYMMDD format

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
echo "Command: make_cobalt_quicklooks.py $cmd_args"
echo "============================================================"

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Create log directory if it doesn't exist
mkdir -p slurm_logs

# Execute the quicklooks script
time python $SCRIPT_DIR/make_cobalt_quicklooks.py $cmd_args

exit_code=$?

echo "============================================================"
echo "Job completed with exit code: $exit_code"
echo "Processing time: $(date)"
echo "============================================================"

exit $exit_code