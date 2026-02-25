#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=coalesc3-test
#SBATCH -o slurm_logs/test_%j.out
#SBATCH -e slurm_logs/test_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path
SCRIPT_DIR="/home/users/cjwalden/git/kepler-radar-utils-cobalt"
PYTHON_SCRIPT="$SCRIPT_DIR/proc_kepler_coalesc3_campaign_batch.py"

# Test with a single date
DATESTR=20170310

echo "Testing COALESC3 processing for date: ${DATESTR}"
echo "SLURM Job ID: ${SLURM_JOB_ID}"

# Process the date
time python $PYTHON_SCRIPT -d ${DATESTR} --skip-missing --single-sweep --gzip --data-version 1.0.1

echo "Completed test processing for ${DATESTR}"
