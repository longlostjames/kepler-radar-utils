#!/bin/bash
# submit_ccrest_m_range.sh - Helper script to submit CCREST-M date range processing

START_DATE=${1:-20240201}
END_DATE=${2:-20240416}
STATIONARY_SCAN_RATE_THRESHOLD=${STATIONARY_SCAN_RATE_THRESHOLD:-0.2}
POINTING_SCAN_RATE_THRESHOLD=${POINTING_SCAN_RATE_THRESHOLD:-0.2}

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
echo "Using stationary scan-rate threshold: ${STATIONARY_SCAN_RATE_THRESHOLD}"
echo "Using pointing scan-rate threshold: ${POINTING_SCAN_RATE_THRESHOLD}"

# Submit with correct array size
START_DATE=${START_DATE} END_DATE=${END_DATE} \
STATIONARY_SCAN_RATE_THRESHOLD=${STATIONARY_SCAN_RATE_THRESHOLD} \
POINTING_SCAN_RATE_THRESHOLD=${POINTING_SCAN_RATE_THRESHOLD} \
sbatch --array=1-${DAYS} kepler_ccrest_m_lotus2_campaign.sh
