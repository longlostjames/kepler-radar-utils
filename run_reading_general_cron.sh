#!/bin/bash
#
# run_reading_general_cron.sh
#
# Non-interactive script for cron: processes today's data from the 'latest'
# subdirectory, generates VPT quicklooks, and regenerates the index.html.
#
# Intended to be run every 30 minutes via cron, e.g.:
#   */30 * * * * /home/cwalden/git/kepler-radar-utils/run_reading_general_cron.sh >> /data/processing/kepler/reading-general/logs/cron.log 2>&1
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATESTR=$(date -u +%Y%m%d)
DATA_VERSION=${DATA_VERSION:-1.0.0}

INPATH="/data/processing/kepler/reading-general"
OUTPATH="$INPATH/quicklooks"
LOG_DIR="$INPATH/logs"
mkdir -p "$LOG_DIR"

echo "========================================"
echo "Reading General cron run: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "Processing date: $DATESTR"
echo "========================================"

# Load conda
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate cao_3_11

# ── 1. Process NetCDF files ───────────────────────────────────────────────────
echo ""
echo "--- Processing ---"
python "$SCRIPT_DIR/proc_kepler_reading_general_campaign_batch.py" \
    -d "$DATESTR" \
    --data-version "$DATA_VERSION" \
    --skip-missing \
    --single-sweep \
    --latest \
    --force

# ── 2. Generate VPT quicklook ─────────────────────────────────────────────────
echo ""
echo "--- Quicklooks ---"
python "$SCRIPT_DIR/make_reading_general_quicklooks.py" \
    -d "$DATESTR" \
    -i "$INPATH" \
    -o "$OUTPATH"

# ── 3. Regenerate index.html ──────────────────────────────────────────────────
echo ""
echo "--- Regenerating index ---"
python "$SCRIPT_DIR/generate_quicklooks_index.py" -o "$OUTPATH"

echo ""
echo "Done: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
