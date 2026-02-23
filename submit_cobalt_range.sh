#!/bin/bash
# submit_cobalt_range.sh - Helper script to submit date range processing

START_DATE=${1:-20241210}
END_DATE=${2:-20251201}

# Calculate number of days
DAYS=$(python -c "
import datetime
start = datetime.datetime.strptime('${START_DATE}', '%Y%m%d')
end = datetime.datetime.strptime('${END_DATE}', '%Y%m%d')
days = (end - start).days + 1
print(days)
")

echo "Processing date range: ${START_DATE} to ${END_DATE} (${DAYS} days)"
echo "Submitting SLURM array job with ${DAYS} tasks..."

# Submit with correct array size
START_DATE=${START_DATE} END_DATE=${END_DATE} sbatch --array=1-${DAYS} kepler_cobalt_lotus2_campaign.sh