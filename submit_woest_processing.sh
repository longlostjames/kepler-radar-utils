#!/bin/bash
#
# Script to submit WOEST processing jobs to SLURM
# Run this from sci-vm-05 or another JASMIN sci VM
#

echo "WOEST Campaign Data Processing Submission"
echo "=========================================="
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

# Campaign date range: Update these to match actual WOEST campaign dates
# Using placeholder dates: 2023-06-01 to 2023-08-31 = 92 days
START_DATE=20230601
END_DATE=20230831

echo "Campaign dates: $START_DATE to $END_DATE"
echo "Total days: 92"
echo ""
echo "NOTE: Update START_DATE and END_DATE in this script to match actual WOEST campaign dates"
echo ""

# Ask user what to do
echo "Choose submission option:"
echo "  1) Test single date (20230615) on test partition"
echo "  2) Process first 10 dates (testing)"
echo "  3) Process entire campaign (92 dates)"
echo "  4) Custom array range"
echo "  5) Custom date range"
echo ""
read -p "Enter choice [1-5]: " choice

case $choice in
    1)
        echo "Submitting test job for single date..."
        sbatch test_woest_single.sh
        ;;
    2)
        echo "Submitting array job for first 10 dates..."
        START_DATE=$START_DATE END_DATE=$END_DATE sbatch --array=1-10 kepler_woest_lotus2_campaign.sh
        ;;
    3)
        echo "Submitting array job for entire campaign (92 dates)..."
        echo "This will submit 92 array tasks"
        read -p "Are you sure? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            START_DATE=$START_DATE END_DATE=$END_DATE sbatch --array=1-92 kepler_woest_lotus2_campaign.sh
        else
            echo "Cancelled"
        fi
        ;;
    4)
        read -p "Enter array range (e.g., 1-30 or 1,5,10-20): " array_range
        echo "Submitting array job with range: $array_range"
        START_DATE=$START_DATE END_DATE=$END_DATE sbatch --array=$array_range kepler_woest_lotus2_campaign.sh
        ;;
    5)
        read -p "Enter start date (YYYYMMDD): " custom_start
        read -p "Enter end date (YYYYMMDD): " custom_end
        
        # Calculate number of days
        DAYS=$(python -c "
import datetime
start = datetime.datetime.strptime('${custom_start}', '%Y%m%d')
end = datetime.datetime.strptime('${custom_end}', '%Y%m%d')
days = (end - start).days + 1
print(days)
")
        
        echo "Processing date range: $custom_start to $custom_end ($DAYS days)"
        read -p "Submit array job with $DAYS tasks? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            START_DATE=$custom_start END_DATE=$custom_end sbatch --array=1-$DAYS kepler_woest_lotus2_campaign.sh
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
