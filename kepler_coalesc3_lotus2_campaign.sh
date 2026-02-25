#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=coalesc3-processing
#SBATCH -o slurm_logs/coalesc3_%A_%a.out
#SBATCH -e slurm_logs/coalesc3_%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path
SCRIPT_DIR="/home/users/cjwalden/git/kepler-radar-utils-cobalt"
PYTHON_SCRIPT="$SCRIPT_DIR/proc_kepler_coalesc3_campaign_batch.py"

# Set date range and output path (can be overridden via environment variables)
START_DATE=${START_DATE:-20170308}
END_DATE=${END_DATE:-20170705}
OUTPATH=${OUTPATH:-/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/coalesc3/L1_v1.0.1}

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
time python $PYTHON_SCRIPT -d ${DATESTR} --skip-missing --single-sweep --gzip --data-version 1.0.1 --outpath ${OUTPATH}

echo "Completed processing for ${DATESTR}"
