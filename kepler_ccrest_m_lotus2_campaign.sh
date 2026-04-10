#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=ccrest-m-campaign
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --mem=128G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_SCRIPT="$SCRIPT_DIR/proc_kepler_ccrest_m_campaign_batch.py"

# Set date range and output path (can be overridden via environment variables)
START_DATE=${START_DATE:-20240201}
END_DATE=${END_DATE:-20240416}
OUTPATH=${OUTPATH:-/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/ccrest-m/L1_v1.0.1}
DATA_VERSION=${DATA_VERSION:-1.0.1}

# Optional PPI/pointing classifier tuning (override via environment variables)
STATIONARY_THRESHOLD=${STATIONARY_THRESHOLD:-0.5}
STATIONARY_STD_THRESHOLD=${STATIONARY_STD_THRESHOLD:-1.0}
STATIONARY_RANGE_THRESHOLD=${STATIONARY_RANGE_THRESHOLD:-2.0}
STATIONARY_TOTAL_CHANGE_THRESHOLD=${STATIONARY_TOTAL_CHANGE_THRESHOLD:-3.0}
STATIONARY_SCAN_RATE_THRESHOLD=${STATIONARY_SCAN_RATE_THRESHOLD:-0.2}
DIRECTION_RATE_THRESHOLD=${DIRECTION_RATE_THRESHOLD:-0.1}
POINTING_WINDOW_SIZE=${POINTING_WINDOW_SIZE:-10}
POINTING_AZ_STD_THRESHOLD=${POINTING_AZ_STD_THRESHOLD:-1.0}
MIN_POINTING_RAYS=${MIN_POINTING_RAYS:-10}
POINTING_TOTAL_CHANGE_THRESHOLD=${POINTING_TOTAL_CHANGE_THRESHOLD:-3.0}
POINTING_SCAN_RATE_THRESHOLD=${POINTING_SCAN_RATE_THRESHOLD:-0.2}

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
echo "Data version: ${DATA_VERSION}"
echo "Classifier tuning: stationary_total_change=${STATIONARY_TOTAL_CHANGE_THRESHOLD}, pointing_total_change=${POINTING_TOTAL_CHANGE_THRESHOLD}"
echo "Classifier tuning: stationary_scan_rate=${STATIONARY_SCAN_RATE_THRESHOLD}, pointing_scan_rate=${POINTING_SCAN_RATE_THRESHOLD}"

# Process the date
time python $PYTHON_SCRIPT \
    -d ${DATESTR} \
    --skip-missing \
    --single-sweep \
    --gzip \
    --data-version ${DATA_VERSION} \
    --outpath ${OUTPATH} \
    --stationary-threshold ${STATIONARY_THRESHOLD} \
    --stationary-std-threshold ${STATIONARY_STD_THRESHOLD} \
    --stationary-range-threshold ${STATIONARY_RANGE_THRESHOLD} \
    --stationary-total-change-threshold ${STATIONARY_TOTAL_CHANGE_THRESHOLD} \
    --stationary-scan-rate-threshold ${STATIONARY_SCAN_RATE_THRESHOLD} \
    --direction-rate-threshold ${DIRECTION_RATE_THRESHOLD} \
    --pointing-window-size ${POINTING_WINDOW_SIZE} \
    --pointing-az-std-threshold ${POINTING_AZ_STD_THRESHOLD} \
    --min-pointing-rays ${MIN_POINTING_RAYS} \
    --pointing-total-change-threshold ${POINTING_TOTAL_CHANGE_THRESHOLD} \
    --pointing-scan-rate-threshold ${POINTING_SCAN_RATE_THRESHOLD}

echo "Completed processing for ${DATESTR}"
