#!/bin/bash
#
# Script to submit CCREST-M processing jobs to SLURM
# Run this from sci-vm-05 or another JASMIN sci VM
#

echo "CCREST-M Campaign Data Processing Submission"
echo "============================================"
echo ""

# Check we're in the right directory (self-locating)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -d "$SCRIPT_DIR" ]; then
    echo "Error: Script directory not found: $SCRIPT_DIR"
    exit 1
fi

cd "$SCRIPT_DIR"

# Ensure slurm_logs directory exists
mkdir -p slurm_logs

# Campaign date range - CCREST-M: February 2024 to April 2024
START_DATE=20240203
END_DATE=20240408

# Output path can be overridden for reprocessing/versioning
OUTPATH=${OUTPATH:-/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/ccrest-m/L1_v1.0.1}
STATIONARY_SCAN_RATE_THRESHOLD=${STATIONARY_SCAN_RATE_THRESHOLD:-0.2}
POINTING_SCAN_RATE_THRESHOLD=${POINTING_SCAN_RATE_THRESHOLD:-0.2}
read -p "Output path [${OUTPATH}]: " input_outpath
OUTPATH=${input_outpath:-$OUTPATH}
echo "Using output path: $OUTPATH"
echo "Using stationary scan-rate threshold: $STATIONARY_SCAN_RATE_THRESHOLD"
echo "Using pointing scan-rate threshold: $POINTING_SCAN_RATE_THRESHOLD"

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
        if [ -f "test_ccrest_m_single.sh" ]; then
            STATIONARY_SCAN_RATE_THRESHOLD=$STATIONARY_SCAN_RATE_THRESHOLD \
            POINTING_SCAN_RATE_THRESHOLD=$POINTING_SCAN_RATE_THRESHOLD \
            sbatch test_ccrest_m_single.sh
        else
            echo "test_ccrest_m_single.sh not found. Create it first."
            exit 1
        fi
        ;;
    2)
        echo "Submitting array job for first 10 dates..."
        OUTPATH=$OUTPATH \
        STATIONARY_SCAN_RATE_THRESHOLD=$STATIONARY_SCAN_RATE_THRESHOLD \
        POINTING_SCAN_RATE_THRESHOLD=$POINTING_SCAN_RATE_THRESHOLD \
        sbatch --array=1-10 kepler_ccrest_m_lotus2_campaign.sh
        ;;
    3)
        echo "Submitting array job for entire campaign ($TOTAL_DAYS dates)..."
        echo "This will submit $TOTAL_DAYS array tasks"
        read -p "Are you sure? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            OUTPATH=$OUTPATH \
            STATIONARY_SCAN_RATE_THRESHOLD=$STATIONARY_SCAN_RATE_THRESHOLD \
            POINTING_SCAN_RATE_THRESHOLD=$POINTING_SCAN_RATE_THRESHOLD \
            sbatch --array=1-$TOTAL_DAYS kepler_ccrest_m_lotus2_campaign.sh
        else
            echo "Cancelled"
        fi
        ;;
    4)
        read -p "Enter array range (e.g., 1-30 or 1,5,10-20): " array_range
        echo "Submitting array job with range: $array_range"
        OUTPATH=$OUTPATH \
        STATIONARY_SCAN_RATE_THRESHOLD=$STATIONARY_SCAN_RATE_THRESHOLD \
        POINTING_SCAN_RATE_THRESHOLD=$POINTING_SCAN_RATE_THRESHOLD \
        sbatch --array=$array_range kepler_ccrest_m_lotus2_campaign.sh
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
            START_DATE=$custom_start END_DATE=$custom_end OUTPATH=$OUTPATH \
            STATIONARY_SCAN_RATE_THRESHOLD=$STATIONARY_SCAN_RATE_THRESHOLD \
            POINTING_SCAN_RATE_THRESHOLD=$POINTING_SCAN_RATE_THRESHOLD \
            sbatch --array=1-$custom_days kepler_ccrest_m_lotus2_campaign.sh
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
