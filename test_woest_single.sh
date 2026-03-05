#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=woest-test
#SBATCH -o slurm_logs/test_%j.out
#SBATCH -e slurm_logs/test_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_SCRIPT="$SCRIPT_DIR/proc_kepler_woest_campaign_batch.py"

# Test with a single date (use a date from WOEST campaign)
# Update this to an actual date when you know the campaign dates
DATESTR=${1:-20230615}

echo "Testing WOEST processing for date: ${DATESTR}"
echo "SLURM Job ID: ${SLURM_JOB_ID}"

# Process the date
time python $PYTHON_SCRIPT -d ${DATESTR} --skip-missing --gzip --data-version 1.0.1

echo "Completed test processing for ${DATESTR}"
