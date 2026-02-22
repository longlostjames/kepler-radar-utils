#!/bin/bash
# Helper script to automatically calculate array range and submit KASBEX processing

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 START_DATE END_DATE [OPTIONS]"
    echo "Example: $0 20230501 20230531"
    echo "         $0 20230501 20230531 --single-sweep"
    echo ""
    echo "Options:"
    echo "  --single-sweep    Create separate files for each sweep"
    echo "  --gzip           Input files are gzip compressed"
    echo "  --force          Overwrite existing output files"
    echo ""
    exit 1
fi

START_DATE=$1
END_DATE=$2
OPTIONS=${3:-""}

# Name of the actual SLURM submit script
SUBMIT_SCRIPT="submit_kasbex_processing.sh"

echo "============================================================"
echo "KASBEX Processing Submission Helper"
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
    echo "Use YYYYMMDD format (e.g., 20230501)"
    exit 1
fi

if ! date -d "$END_DATE" >/dev/null 2>&1; then
    echo "Error: Invalid end date format: $END_DATE"
    echo "Use YYYYMMDD format (e.g., 20230531)"
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
echo "Options: $OPTIONS"
echo ""

# Create backup of submit script
cp "$SUBMIT_SCRIPT" "${SUBMIT_SCRIPT}.backup.$(date +%Y%m%d_%H%M%S)"

# Update the submit script with new dates
echo "Updating $SUBMIT_SCRIPT with new dates..."
sed -i "s/START_DATE=\".*\"/START_DATE=\"$START_DATE\"/" "$SUBMIT_SCRIPT"
sed -i "s/END_DATE=\".*\"/END_DATE=\"$END_DATE\"/" "$SUBMIT_SCRIPT"

# Handle options
if [[ "$OPTIONS" == *"--single-sweep"* ]]; then
    sed -i "s/SINGLE_SWEEP_FLAG=\".*\"/SINGLE_SWEEP_FLAG=\"--single-sweep\"/" "$SUBMIT_SCRIPT"
    echo "Enabled single-sweep mode"
fi

if [[ "$OPTIONS" == *"--gzip"* ]]; then
    sed -i "s/GZIP_FLAG=\".*\"/GZIP_FLAG=\"--gzip\"/" "$SUBMIT_SCRIPT"
    echo "Enabled gzip input mode"
fi

if [[ "$OPTIONS" == *"--force"* ]]; then
    sed -i "s/FORCE_FLAG=\".*\"/FORCE_FLAG=\"--force\"/" "$SUBMIT_SCRIPT"
    echo "Enabled force overwrite mode"
fi

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

# Check available data before submitting
echo ""
echo "Checking for available KASBEX data..."
python -c "
import datetime
from pathlib import Path

start_date = datetime.datetime.strptime('$START_DATE', '%Y%m%d')
end_date = datetime.datetime.strptime('$END_DATE', '%Y%m%d')
inpath = Path('/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-mobile-ka-band-radar-1/data/campaign/kasbex/mom')

total_days = 0
data_days = 0
current_date = start_date

while current_date <= end_date:
    datestr = current_date.strftime('%Y%m%d')
    date_path = inpath / datestr
    
    total_days += 1
    if date_path.exists() and list(date_path.glob('*.mmclx*')):
        data_days += 1
        if data_days <= 5:  # Show first 5 dates with data
            print(f'  ✓ {datestr} - Data available')
    elif data_days <= 5:
        print(f'  ✗ {datestr} - No data')
        
    current_date += datetime.timedelta(days=1)

print(f'\\nData availability: {data_days}/{total_days} days have data')
if data_days < total_days:
    print('Note: --skip-missing is enabled to handle missing data')
"

# Submit the job
echo ""
echo "Submitting KASBEX processing job..."
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
    echo "Output files will be in: /gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/kasbex/L1c/"
else
    echo "✗ Job submission failed!"
    echo "Error output: $sbatch_output"
    exit 1
fi