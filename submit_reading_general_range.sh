#!/bin/bash
# submit_reading_general_range.sh
# Helper script to process a date range for the Reading General campaign.
#
# Usage: ./submit_reading_general_range.sh [START_DATE] [END_DATE]
#   e.g. ./submit_reading_general_range.sh 20260401 20260430

START_DATE=${1:-20260401}
END_DATE=${2:-20261231}

# Calculate number of days
DAYS=$(python -c "
import datetime
start = datetime.datetime.strptime('${START_DATE}', '%Y%m%d')
end = datetime.datetime.strptime('${END_DATE}', '%Y%m%d')
days = (end - start).days + 1
print(days)
")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

echo "Reading General processing: ${START_DATE} to ${END_DATE} (${DAYS} days)"

python $SCRIPT_DIR/proc_kepler_reading_general_campaign_batch.py \
    --start-date ${START_DATE} \
    --end-date ${END_DATE} \
    --skip-missing \
    --single-sweep \
    --latest \
    --data-version 1.0.0
