#!/bin/bash
# run_kepler_kasbex_proc_with_cleanup.sh
# Enhanced wrapper script that includes file cleanup

# Path to your conda installation
CONDA_BASE="$HOME/miniforge3"
SCRIPT_DIR="$HOME/git/kepler-radar-utils-cobalt"

# Get today's date or use the provided date
if [ -n "$1" ]; then
    TODAY="$1"
else
    TODAY=$(date +%Y%m%d)
fi

# Optional: Set max-age in hours (e.g., 6 hours)
MAX_AGE_HOURS=5

# Function to log messages
log_message() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] $level: $message"
}

log_message "INFO" "Starting KASBEX processing for $TODAY"

# Activate conda
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate cao_3_11

# Step 1: Clean up any existing numbered files before processing
log_message "INFO" "Step 1: Cleaning up existing numbered files"
bash "$SCRIPT_DIR/cleanup_kepler_kasbex_files.sh" -d "$TODAY" --verbose

if [ $? -ne 0 ]; then
    log_message "WARN" "File cleanup had issues, but continuing with processing"
fi

# Step 2: Run the processing
log_message "INFO" "Step 2: Running KASBEX data processing"
python "$SCRIPT_DIR/proc_kepler_kasbex_campaign_batch.py" \
    -d "$TODAY" \
    --gzip \
    --skip-missing \
    --single-sweep \
    --latest \
    --force \
    --no-vpt \
    --max-age "$MAX_AGE_HOURS"

if [ $? -ne 0 ]; then
    log_message "ERROR" "Data processing failed"
    conda deactivate
    exit 1
fi

# Step 3: Clean up numbered files created by --force option
log_message "INFO" "Step 3: Cleaning up numbered files created by processing"
bash "$SCRIPT_DIR/cleanup_kepler_kasbex_files.sh" -d "$TODAY" --verbose

if [ $? -ne 0 ]; then
    log_message "WARN" "Post-processing cleanup had issues"
fi

# Step 4: Generate quicklooks
log_message "INFO" "Step 4: Generating quicklooks"
python "$SCRIPT_DIR/make_kepler_kasbex_quicklooks.py" \
    -d "$TODAY" \
    --max-age "$MAX_AGE_HOURS"  # Add this line

if [ $? -ne 0 ]; then
    log_message "ERROR" "Quicklook generation failed"
    conda deactivate
    exit 1
fi

# Deactivate conda environment  
conda deactivate

log_message "INFO" "KASBEX processing completed successfully for $TODAY"

# End of script