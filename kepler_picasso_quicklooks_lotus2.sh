#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=picasso-quicklooks
#SBATCH -o slurm_logs/picasso_quicklooks_%A_%a.out
#SBATCH -e slurm_logs/picasso_quicklooks_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# PICASSO Quicklook Generation - SLURM Array Job
#
# This script processes a range of dates using SLURM array jobs.
# Each array task generates quicklooks for one day.

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_SCRIPT="$SCRIPT_DIR/make_picasso_quicklooks.py"

# Set date range (can be overridden via environment variables)
# PICASSO campaign: December 2017 to May 2019
START_DATE=${START_DATE:-20171212}
END_DATE=${END_DATE:-20190523}

# Calculate the date for this array task
DATESTR=$(python -c "
import datetime
start = datetime.datetime.strptime('${START_DATE}', '%Y%m%d')
target = start + datetime.timedelta(days=${SLURM_ARRAY_TASK_ID}-1)
end = datetime.datetime.strptime('${END_DATE}', '%Y%m%d')
if target <= end:
    print(target.strftime('%Y%m%d'))
else:
    exit(1)
")

# Exit gracefully if beyond date range
if [ $? -ne 0 ]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} beyond date range. Exiting gracefully."
    exit 0
fi

echo "============================================================"
echo "PICASSO Quicklook Generation"
echo "============================================================"
echo "Processing date: ${DATESTR}"
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Python script: ${PYTHON_SCRIPT}"
echo "============================================================"

# Run the quicklook script
# Uses default paths from make_picasso_quicklooks.py
time python $PYTHON_SCRIPT -d ${DATESTR}

exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "============================================================"
    echo "Completed quicklooks for ${DATESTR}"
    echo "============================================================"
else
    echo "============================================================"
    echo "Quicklooks failed for ${DATESTR} with exit code $exit_code"
    echo "============================================================"
fi

exit $exit_code
