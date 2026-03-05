#!/bin/bash
#
# Script to submit PICASSO quicklook generation jobs to SLURM
# Run this from sci-vm-05 or another JASMIN sci VM
#

echo "PICASSO Campaign Quicklook Generation Submission"
echo "================================================"
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

# Campaign date range - PICASSO: December 2017 to May 2019
START_DATE=20171212
END_DATE=20190523

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
echo "  1) Test single date (20171214) on short-serial partition"
echo "  2) Generate quicklooks for first 10 dates (testing)"
echo "  3) Generate quicklooks for entire campaign ($TOTAL_DAYS dates)"
echo "  4) Custom array range"
echo "  5) Custom date range"
echo ""
read -p "Enter choice [1-5]: " choice

case $choice in
    1)
        echo "Submitting test job for single date (20171214)..."
        DATE=20171214 sbatch make_picasso_quicklooks.sh
        ;;
    2)
        echo "Submitting array job for first 10 dates..."
        sbatch --array=1-10 kepler_picasso_quicklooks_lotus2.sh
        ;;
    3)
        echo "Generating quicklooks for entire campaign ($TOTAL_DAYS dates)..."
        echo "This will submit $TOTAL_DAYS array tasks"
        read -p "Are you sure? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            sbatch --array=1-$TOTAL_DAYS kepler_picasso_quicklooks_lotus2.sh
        else
            echo "Cancelled"
        fi
        ;;
    4)
        read -p "Enter array range (e.g., 1-30 or 1,5,10-20): " array_range
        echo "Submitting array job with range: $array_range"
        sbatch --array=$array_range kepler_picasso_quicklooks_lotus2.sh
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
        
        echo "Generating quicklooks for $custom_days days from $custom_start to $custom_end"
        read -p "Submit? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            START_DATE=$custom_start END_DATE=$custom_end sbatch --array=1-$custom_days kepler_picasso_quicklooks_lotus2.sh
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
echo "Job(s) submitted. Check status with: squeue -u $USER"
echo "Monitor logs in: $SCRIPT_DIR/slurm_logs/"
echo ""
echo "Quicklook output location:"
echo "/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/picasso/L1_v1.0.0/quicklooks/"
