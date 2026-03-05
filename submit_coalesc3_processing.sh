#!/bin/bash
#
# Script to submit COALESC3 processing jobs to SLURM
# Run this from sci-vm-05 or another JASMIN sci VM
#

echo "COALESC3 Campaign Data Processing Submission"
echo "============================================="
echo ""

# Check we're in the right directory (self-locating)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -d "$SCRIPT_DIR" ]; then
    echo "Error: Script directory not found: $SCRIPT_DIR"
    exit 1
fi

cd $SCRIPT_DIR

# Ensure slurm_logs directory exists
mkdir -p slurm_logs

# Campaign date range: 2017-03-08 to 2017-07-05 = 120 days
START_DATE=20170308
END_DATE=20170705

echo "Campaign dates: $START_DATE to $END_DATE"
echo "Total days: 120"
echo ""

# Ask user what to do
echo "Choose submission option:"
echo "  1) Test single date (20170310) on test partition"
echo "  2) Process first 10 dates (testing)"
echo "  3) Process entire campaign (120 dates)"
echo "  4) Custom array range"
echo ""
read -p "Enter choice [1-4]: " choice

case $choice in
    1)
        echo "Submitting test job for single date..."
        sbatch test_coalesc3_single.sh
        ;;
    2)
        echo "Submitting array job for first 10 dates..."
        sbatch --array=1-10 kepler_coalesc3_lotus2_campaign.sh
        ;;
    3)
        echo "Submitting array job for entire campaign (120 dates)..."
        echo "This will submit 120 array tasks"
        read -p "Are you sure? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            sbatch --array=1-120 kepler_coalesc3_lotus2_campaign.sh
        else
            echo "Cancelled"
        fi
        ;;
    4)
        read -p "Enter array range (e.g., 1-30 or 1,5,10-20): " array_range
        echo "Submitting array job with range: $array_range"
        sbatch --array=$array_range kepler_coalesc3_lotus2_campaign.sh
        ;;
    *)
        echo "Invalid choice"
        exit 1
        ;;
esac

echo ""
echo "Job submitted. Check status with: squeue -u $USER"
echo "Monitor logs in: $SCRIPT_DIR/slurm_logs/"
