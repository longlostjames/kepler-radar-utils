#!/bin/bash
#
# submit_reading_general_processing.sh
#
# Interactive script to submit Reading General processing jobs.
# Run this from the kepler-radar-utils directory.
#

echo "Reading General Campaign Data Processing Submission"
echo "===================================================="
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -d "$SCRIPT_DIR" ]; then
    echo "Error: Script directory not found: $SCRIPT_DIR"
    exit 1
fi

cd $SCRIPT_DIR

# Load conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Campaign date range defaults
START_DATE=${START_DATE:-20260401}
END_DATE=${END_DATE:-20261231}

# Product version - used to define the output path (L1_v<DATA_VERSION>)
DATA_VERSION=${DATA_VERSION:-1.0.0}
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
echo "Choose processing option:"
echo "  1) Process single date (local)"
echo "  2) Process entire campaign [$TOTAL_DAYS dates] (local)"
echo "  3) Process custom date range (local)"
echo "  4) Dry run (show what would be processed)"
echo ""
read -p "Enter choice [1-4]: " choice

PYTHON_SCRIPT="$SCRIPT_DIR/proc_kepler_reading_general_campaign_batch.py"
BASE_ARGS="--data-version $DATA_VERSION --skip-missing --single-sweep --latest"

case $choice in
    1)
        read -p "Enter date to process (YYYYMMDD): " single_date
        echo "Processing single date: $single_date"
        python $PYTHON_SCRIPT -d $single_date $BASE_ARGS
        ;;
    2)
        echo "Processing entire campaign ($TOTAL_DAYS dates) from $START_DATE to $END_DATE"
        read -p "Are you sure? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            python $PYTHON_SCRIPT --start-date $START_DATE --end-date $END_DATE $BASE_ARGS
        else
            echo "Cancelled"
        fi
        ;;
    3)
        read -p "Enter start date (YYYYMMDD): " custom_start
        read -p "Enter end date (YYYYMMDD): " custom_end
        custom_days=$(python -c "
import datetime
start = datetime.datetime.strptime('${custom_start}', '%Y%m%d')
end = datetime.datetime.strptime('${custom_end}', '%Y%m%d')
print((end - start).days + 1)
")
        echo "Processing $custom_days days from $custom_start to $custom_end"
        read -p "Submit? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            python $PYTHON_SCRIPT --start-date $custom_start --end-date $custom_end $BASE_ARGS
        else
            echo "Cancelled"
        fi
        ;;
    4)
        read -p "Enter start date for dry run (YYYYMMDD) [${START_DATE}]: " dr_start
        read -p "Enter end date for dry run (YYYYMMDD) [${END_DATE}]: " dr_end
        dr_start=${dr_start:-$START_DATE}
        dr_end=${dr_end:-$END_DATE}
        python $PYTHON_SCRIPT --start-date $dr_start --end-date $dr_end $BASE_ARGS --dry-run
        ;;
    *)
        echo "Invalid choice"
        exit 1
        ;;
esac
