#!/bin/bash
#
# cron_reading_general_processing.sh
#
# Intended to be run daily at 06:00 UTC via cron:
#   0 6 * * * /home/cwalden/git/kepler-radar-utils/cron_reading_general_processing.sh >> /home/cwalden/logs/reading_general_processing.log 2>&1
#
# Processes data for two days ago and overwrites any existing output.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate cao_3_11

TARGET_DATE=$(date -u -d "2 days ago" +%Y%m%d)
DATA_VERSION="1.0.0"

echo "=========================================="
echo "Reading General processing cron"
echo "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "Target date: $TARGET_DATE"
echo "=========================================="

INPATH="/data/processing/kepler/reading-general"
OUTPATH="/data/processing/kepler/reading-general/quicklooks"

python "$SCRIPT_DIR/proc_kepler_reading_general_campaign_batch.py" \
    -d "$TARGET_DATE" \
    --data-version "$DATA_VERSION" \
    --skip-missing \
    --single-sweep \
    --force \
    --gzip

echo "Processing complete: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo ""
echo "Generating quicklooks for $TARGET_DATE..."

python "$SCRIPT_DIR/make_reading_general_quicklooks.py" \
    -d "$TARGET_DATE" \
    -i "$INPATH" \
    -o "$OUTPATH"

echo "Regenerating quicklooks index..."
python "$SCRIPT_DIR/generate_quicklooks_index.py" -o "$OUTPATH"

echo "Quicklooks complete: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
