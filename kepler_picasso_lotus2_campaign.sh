#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=picasso-campaign
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_SCRIPT="$SCRIPT_DIR/proc_kepler_picasso_campaign_batch.py"

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

echo "Processing date: ${DATESTR}"
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"

# Process the date
# Note: --outpath is auto-generated from --data-version unless explicitly set
# Gzip is enabled by default (files are .mmclx.gz)
# --split-man-phases enables phase detection for MAN scans
time python $PYTHON_SCRIPT -d ${DATESTR} --skip-missing --split-man-phases --data-version 1.0.0

echo "Completed processing for ${DATESTR}"
