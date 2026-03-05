#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=coalesc3_quicklooks_batch
#SBATCH -o slurm_logs/coalesc3_ql_%A_%a.out
#SBATCH -e slurm_logs/coalesc3_ql_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_SCRIPT="$SCRIPT_DIR/make_coalesc3_quicklooks.py"

# Set date range (can be overridden via environment variables)
# COALESC3 campaign ran from March to July 2017
START_DATE=${START_DATE:-20170308}
END_DATE=${END_DATE:-20170705}
INPATH=${INPATH:-/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/coalesc3/L1_v1.0.1}
OUTPATH=${OUTPATH:-$INPATH/quicklooks}

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

echo "Processing COALESC3 quicklooks for date: ${DATESTR}"
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Input path: ${INPATH}"
echo "Output path: ${OUTPATH}"

# Generate quicklooks
time python $PYTHON_SCRIPT -d ${DATESTR} -i ${INPATH} -o ${OUTPATH}

echo "Completed quicklook generation for ${DATESTR}"
