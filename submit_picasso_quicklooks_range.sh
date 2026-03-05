#!/bin/bash
#
# Helper script to submit PICASSO quicklooks for a specific date range
# Usage: ./submit_picasso_quicklooks_range.sh START_DATE END_DATE
# Example: ./submit_picasso_quicklooks_range.sh 20171212 20171231
#

if [ $# -ne 2 ]; then
    echo "Usage: $0 START_DATE END_DATE"
    echo "Example: $0 20171212 20171231"
    echo ""
    echo "PICASSO Campaign: 2017-12-12 to 2019-05-23"
    exit 1
fi

START_DATE=$1
END_DATE=$2

# Validate date format
if ! [[ $START_DATE =~ ^[0-9]{8}$ ]] || ! [[ $END_DATE =~ ^[0-9]{8}$ ]]; then
    echo "Error: Dates must be in YYYYMMDD format"
    exit 1
fi

# Check if the array job script exists
ARRAY_SCRIPT="kepler_picasso_quicklooks_lotus2.sh"
if [ ! -f "$ARRAY_SCRIPT" ]; then
    echo "Error: Array job script '$ARRAY_SCRIPT' not found in current directory"
    exit 1
fi

# Calculate number of days between dates
num_days=$(python -c "
import datetime
import sys
try:
    start = datetime.datetime.strptime('${START_DATE}', '%Y%m%d')
    end = datetime.datetime.strptime('${END_DATE}', '%Y%m%d')
    delta = (end - start).days + 1
    if delta < 1:
        print('Error: End date must be on or after start date', file=sys.stderr)
        sys.exit(1)
    print(delta)
except ValueError as e:
    print(f'Error: Invalid date format - {e}', file=sys.stderr)
    sys.exit(1)
")

if [ $? -ne 0 ]; then
    exit 1
fi

echo "PICASSO Quicklook Submission"
echo "============================="
echo "Start date:  $START_DATE"
echo "End date:    $END_DATE"
echo "Total days:  $num_days"
echo "Array range: 1-$num_days"
echo ""

# Confirm submission
read -p "Submit $num_days quicklook jobs? [y/N]: " confirm
if [[ ! $confirm == [yY] ]]; then
    echo "Cancelled"
    exit 0
fi

# Create slurm_logs directory if it doesn't exist
mkdir -p slurm_logs

# Submit the job with environment variables and array range
echo "Submitting to SLURM..."
START_DATE=$START_DATE END_DATE=$END_DATE sbatch --array=1-$num_days "$ARRAY_SCRIPT"

if [ $? -eq 0 ]; then
    echo ""
    echo "Job submitted successfully!"
    echo ""
    echo "Monitor with: squeue -u $USER"
    echo "Logs: slurm_logs/picasso_quicklooks_*.out/err"
    echo "Output: /gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/picasso/L1_v1.0.0/quicklooks/"
else
    echo ""
    echo "Error: Job submission failed"
    exit 1
fi
