#!/bin/bash
#SBATCH --partition=test
#SBATCH --job-name=picasso-test
#SBATCH -o slurm_logs/test_picasso_%j.out
#SBATCH -e slurm_logs/test_picasso_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --account=ncas_radar
#SBATCH --qos=interactive

# Test script for processing a single PICASSO date
# Update the test date to a date with actual data

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_SCRIPT="$SCRIPT_DIR/proc_kepler_picasso_campaign_batch.py"

# Test date - Known date with MAN scan data for testing
TEST_DATE=20171214

echo "=========================================="
echo "PICASSO Campaign Test Processing"
echo "=========================================="
echo "Test date: ${TEST_DATE}"
echo "Script: ${PYTHON_SCRIPT}"
echo ""

# Process the test date with phase splitting for MAN scans
# Gzip is enabled by default (files are .mmclx.gz)
# --split-man-phases enables phase detection for MAN scans
time python $PYTHON_SCRIPT -d ${TEST_DATE} --split-man-phases --data-version 1.0.0

echo ""
echo "Test processing completed for ${TEST_DATE}"
