#!/bin/bash
# Helper script to automatically calculate array range and submit

if [ $# -ne 2 ]; then
    echo "Usage: $0 START_DATE END_DATE"
    echo "Example: $0 20241210 20251201"
    exit 1
fi

START_DATE=$1
END_DATE=$2

# Name of the actual SLURM submit script
SUBMIT_SCRIPT="make_kepler_cobalt_quicklooks.sh"

echo "============================================================"
echo "COBALT Quicklooks Submission Helper"
echo "============================================================"

# Check if the submit script exists
if [ ! -f "$SUBMIT_SCRIPT" ]; then
    echo "Error: Submit script '$SUBMIT_SCRIPT' not found in current directory"
    echo "Current directory: $(pwd)"
    echo "Available .sh files:"
    ls -la *.sh 2>/dev/null || echo "No .sh files found"
    exit 1
fi

# Validate dates
if ! date -d "$START_DATE" >/dev/null 2>&1; then
    echo "Error: Invalid start date format: $START_DATE"
    echo "Use YYYYMMDD format (e.g., 20241210)"
    exit 1
fi

if ! date -d "$END_DATE" >/dev/null 2>&1; then
    echo "Error: Invalid end date format: $END_DATE"
    echo "Use YYYYMMDD format (e.g., 20250728)"
    exit 1
fi

# Calculate number of days
start_epoch=$(date -d "$START_DATE" +%s)
end_epoch=$(date -d "$END_DATE" +%s)
num_days=$(( ( end_epoch - start_epoch ) / 86400 ))

if [ $num_days -lt 0 ]; then
    echo "Error: End date must be after start date"
    echo "Start: $START_DATE"
    echo "End: $END_DATE"
    exit 1
fi

echo "Date range: $START_DATE to $END_DATE"
echo "Number of days: $((num_days + 1))"
echo "Array range: 0-$num_days"
echo ""

# Create backup of submit script
cp "$SUBMIT_SCRIPT" "${SUBMIT_SCRIPT}.backup.$(date +%Y%m%d_%H%M%S)"

# Update the submit script with new dates
echo "Updating $SUBMIT_SCRIPT with new dates..."
sed -i "s/START_DATE=\".*\"/START_DATE=\"$START_DATE\"/" "$SUBMIT_SCRIPT"
sed -i "s/END_DATE=\".*\"/END_DATE=\"$END_DATE\"/" "$SUBMIT_SCRIPT"

# Verify the updates worked
if grep -q "START_DATE=\"$START_DATE\"" "$SUBMIT_SCRIPT" && grep -q "END_DATE=\"$END_DATE\"" "$SUBMIT_SCRIPT"; then
    echo "✓ Successfully updated dates in $SUBMIT_SCRIPT"
else
    echo "✗ Failed to update dates in $SUBMIT_SCRIPT"
    echo "Please check the file format and try again"
    exit 1
fi

# Create log directory
mkdir -p slurm_logs

# Submit the job
echo ""
echo "Submitting job..."
echo "Command: sbatch --array=0-$num_days $SUBMIT_SCRIPT"
echo ""

sbatch_output=$(sbatch --array=0-$num_days "$SUBMIT_SCRIPT" 2>&1)
sbatch_exit=$?

if [ $sbatch_exit -eq 0 ]; then
    echo "✓ Job submitted successfully!"
    echo "$sbatch_output"
    job_id=$(echo "$sbatch_output" | grep -o '[0-9]*')
    echo ""
    echo "Job monitoring commands:"
    echo "  squeue -u \$USER"
    echo "  squeue -j $job_id"
    echo "  scancel $job_id    # to cancel all array tasks"
    echo ""
    echo "Log files will be in: slurm_logs/"
else
    echo "✗ Job submission failed!"
    echo "Error output: $sbatch_output"
    exit 1
fi