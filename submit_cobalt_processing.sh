#!/bin/bash
#
# Script to submit COBALT processing jobs to SLURM
# Run this from sci-vm-05 or another JASMIN sci VM
#

echo "COBALT Campaign Data Processing Submission"
echo "==========================================="
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

# Campaign date range - COBALT: December 2024 to December 2025
START_DATE=20241210
END_DATE=20251201

# Product version - used to define the output path (L1_v<DATA_VERSION>)
DATA_VERSION=${DATA_VERSION:-1.0.3}
read -p "Product version [${DATA_VERSION}]: " input_version
DATA_VERSION=${input_version:-$DATA_VERSION}
echo "Using product version: $DATA_VERSION"

# Calculate total days
TOTAL_DAYS=$(python -c "
import datetime
start = datetime.datetime.strptime('${START_DATE}', '%Y%m%d')
end = datetime.datetime.strptime('${END_DATE}', '%Y%m%d')
print((end - start).days + 1)
")

echo "Campaign dates: $START_DATE to $END_DATE"
echo "Total days: $TOTAL_DAYS"
echo ""

# Ask user what to do
echo "Choose submission option:"
echo "  1) Test single date on test partition"
echo "  2) Process first 10 dates (testing)"
echo "  3) Process entire campaign ($TOTAL_DAYS dates)"
echo "  4) Custom array range"
echo "  5) Process custom date range"
echo ""
read -p "Enter choice [1-5]: " choice

case $choice in
    1)
        echo "Submitting test job for single date..."
        if [ -f "test_cobalt_single.sh" ]; then
            sbatch test_cobalt_single.sh
        else
            echo "test_cobalt_single.sh not found. Create it first."
            exit 1
        fi
        ;;
    2)
        echo "Submitting array job for first 10 dates..."
        DATA_VERSION=$DATA_VERSION sbatch --array=1-10 kepler_cobalt_lotus2_campaign.sh
        ;;
    3)
        echo "Submitting array job for entire campaign ($TOTAL_DAYS dates)..."
        echo "This will submit $TOTAL_DAYS array tasks"
        read -p "Are you sure? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            DATA_VERSION=$DATA_VERSION sbatch --array=1-$TOTAL_DAYS kepler_cobalt_lotus2_campaign.sh
        else
            echo "Cancelled"
        fi
        ;;
    4)
        read -p "Enter array range (e.g., 1-30 or 1,5,10-20): " array_range
        echo "Submitting array job with range: $array_range"
        DATA_VERSION=$DATA_VERSION sbatch --array=$array_range kepler_cobalt_lotus2_campaign.sh
        ;;
    5)
        read -p "Enter start date (YYYYMMDD): " custom_start
        read -p "Enter end date (YYYYMMDD): " custom_end

        # Calculate days for custom range
        custom_days=$(python -c "
import datetime
start = datetime.datetime.strptime('${custom_start}', '%Y%m%d')
end = datetime.datetime.strptime('${custom_end}', '%Y%m%d')
print((end - start).days + 1)
")

        echo "Processing $custom_days days from $custom_start to $custom_end"
        read -p "Submit? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            START_DATE=$custom_start END_DATE=$custom_end DATA_VERSION=$DATA_VERSION sbatch --array=1-$custom_days kepler_cobalt_lotus2_campaign.sh
        else
            echo "Cancelled"
        fi
        ;;
    *)
        echo "Invalid choice"
        exit 1
        ;;
esac

echo ""
echo "Job submitted. Check status with: squeue -u $USER"
echo "Monitor logs in: $SCRIPT_DIR/slurm_logs/"
