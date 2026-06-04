#!/bin/bash
#
# submit_reading_general_quicklooks.sh
#
# Interactive script to generate Reading General quicklooks locally.
# Run this from the kepler-radar-utils directory.
#

echo "Reading General Quicklooks Submission"
echo "======================================"
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $SCRIPT_DIR

# Load conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Defaults
DEFAULT_INPATH="/data/processing/kepler/reading-general"
DEFAULT_OUTPATH="/data/processing/kepler/reading-general/quicklooks"
START_DATE=${START_DATE:-20260401}
END_DATE=${END_DATE:-20261231}

PYTHON_SCRIPT="$SCRIPT_DIR/make_reading_general_quicklooks.py"
BASE_ARGS="-i $DEFAULT_INPATH -o $DEFAULT_OUTPATH"

echo "Choose option:"
echo "  1) Single date (VPT + VAD + PPI)"
echo "  2) Date range (VPT + VAD + PPI)"
echo "  3) Custom date range (VPT + VAD + PPI)"
echo "  4) Boundary layer mode (single date, 4 km height limit)"
echo "  5) Single date — VPT only"
echo "  6) Single date — PPI only"
echo "  7) Single date — PPI with day-outermost map extent"
echo "  8) Single date — VAD only"
echo "  9) Date range — VAD only"
echo ""
read -p "Enter choice [1-9]: " choice

case $choice in
    1)
        read -p "Enter date (YYYYMMDD): " single_date
        echo "Generating quicklooks for $single_date..."
        python $PYTHON_SCRIPT -d $single_date $BASE_ARGS
        ;;
    2)
        echo "Processing $START_DATE to $END_DATE"
        read -p "Are you sure? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            current=$START_DATE
            while [[ $current -le $END_DATE ]]; do
                echo "--- $current ---"
                python $PYTHON_SCRIPT -d $current $BASE_ARGS
                current=$(date -d "$current + 1 day" +%Y%m%d)
            done
        else
            echo "Cancelled"
        fi
        ;;
    3)
        read -p "Enter start date (YYYYMMDD): " custom_start
        read -p "Enter end date (YYYYMMDD): " custom_end
        echo "Processing $custom_start to $custom_end"
        read -p "Submit? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            current=$custom_start
            while [[ $current -le $custom_end ]]; do
                echo "--- $current ---"
                python $PYTHON_SCRIPT -d $current $BASE_ARGS
                current=$(date -d "$current + 1 day" +%Y%m%d)
            done
        else
            echo "Cancelled"
        fi
        ;;
    4)
        read -p "Enter date (YYYYMMDD): " single_date
        echo "Generating BL quicklooks for $single_date..."
        python $PYTHON_SCRIPT -d $single_date $BASE_ARGS -b
        ;;
    5)
        read -p "Enter date (YYYYMMDD): " single_date
        echo "Generating VPT-only quicklooks for $single_date..."
        python $PYTHON_SCRIPT -d $single_date $BASE_ARGS --vpt-only
        ;;
    6)
        read -p "Enter date (YYYYMMDD): " single_date
        echo "Generating PPI-only quicklooks for $single_date..."
        python $PYTHON_SCRIPT -d $single_date $BASE_ARGS --ppi-only
        ;;
    7)
        read -p "Enter date (YYYYMMDD): " single_date
        echo "Generating PPI-only quicklooks (day-outermost extent) for $single_date..."
        python $PYTHON_SCRIPT -d $single_date $BASE_ARGS --ppi-only --ppi-map-day-extent
        ;;
    8)
        read -p "Enter date (YYYYMMDD): " single_date
        echo "Generating VAD-only quicklooks for $single_date..."
        python $PYTHON_SCRIPT -d $single_date $BASE_ARGS --vad-only
        ;;
    9)
        read -p "Enter start date (YYYYMMDD): " custom_start
        read -p "Enter end date (YYYYMMDD): " custom_end
        echo "Processing VAD quicklooks $custom_start to $custom_end"
        read -p "Submit? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            current=$custom_start
            while [[ $current -le $custom_end ]]; do
                echo "--- $current ---"
                python $PYTHON_SCRIPT -d $current $BASE_ARGS --vad-only
                current=$(date -d "$current + 1 day" +%Y%m%d)
            done
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
echo "Regenerating quicklooks index..."
python "$SCRIPT_DIR/generate_quicklooks_index.py" -o "$DEFAULT_OUTPATH"

echo ""
echo "Done."
