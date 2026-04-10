#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=cobalt-campaign
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --mem=256G
#SBATCH --account=ncas_radar
#SBATCH --qos=high

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_SCRIPT="$SCRIPT_DIR/proc_kepler_cobalt_campaign_batch.py"

# Set date range and product version (can be overridden via environment variables)
START_DATE=${START_DATE:-20241210}
END_DATE=${END_DATE:-20251201}
DATA_VERSION=${DATA_VERSION:-1.0.3}

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
time python $PYTHON_SCRIPT -d ${DATESTR} --skip-missing --single-sweep --gzip --data-version ${DATA_VERSION}

echo "Completed processing for ${DATESTR}"
