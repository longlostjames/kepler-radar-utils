#!/bin/bash
#SBATCH --partition=test
#SBATCH --job-name=ccrest-m-test
#SBATCH -o slurm_logs/test_ccrest_m_%j.out
#SBATCH -e slurm_logs/test_ccrest_m_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --account=ncas_radar
#SBATCH --qos=interactive

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_SCRIPT="$SCRIPT_DIR/proc_kepler_ccrest_m_campaign_batch.py"

# Test date can be overridden by first positional argument
DATESTR=${1:-20240215}

# Output path can be overridden via environment variable
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

echo "=========================================="
echo "CCREST-M Campaign Test Processing"
echo "=========================================="
echo "Test date: ${DATESTR}"
echo "Output path: ${OUTPATH}"
echo "Data version: ${DATA_VERSION}"
echo "Script: ${PYTHON_SCRIPT}"
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo ""

# Process the test date
# Match campaign defaults used in kepler_ccrest_m_lotus2_campaign.sh
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

echo ""
echo "Completed test processing for ${DATESTR}"
