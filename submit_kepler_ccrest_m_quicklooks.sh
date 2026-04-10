#!/bin/bash
# Helper script to automatically calculate array range and submit CCREST-M quicklooks jobs

if [ $# -lt 2 ]; then
    echo "Usage: $0 START_DATE END_DATE [-b] [--skip-all-transition] [--vpt-only] [--ppi-only] [--pointing-only] [--ppi-map-day-extent] [--ppi-map-colorbar-shrink VALUE]"
    echo "Example: $0 20240201 20240416"
    echo "         $0 20240201 20240416 -b"
    echo "         $0 20240201 20240416 --skip-all-transition"
    echo "         $0 20240201 20240416 -b --skip-all-transition"
    echo "         $0 20240201 20240416 --vpt-only"
    echo "         $0 20240201 20240416 --ppi-only"
    echo "         $0 20240201 20240416 --ppi-only --ppi-map-day-extent"
    echo "         $0 20240201 20240416 --ppi-only --ppi-map-colorbar-shrink 0.82"
    exit 1
fi

START_DATE=$1
END_DATE=$2
shift 2

BOUNDARY_LAYER_FLAG=""
SKIP_ALL_TRANSITION_FLAG=""
VPT_ONLY_FLAG=""
PPI_ONLY_FLAG=""
POINTING_ONLY_FLAG=""
PPI_MAP_DAY_EXTENT_FLAG=""
PPI_MAP_COLORBAR_SHRINK=""

# Parse optional flags
for arg in "$@"; do
    case "$arg" in
        -b)
            BOUNDARY_LAYER_FLAG="-b"
            ;;
        --skip-all-transition)
            SKIP_ALL_TRANSITION_FLAG="--skip-all-transition"
            ;;
        --vpt-only)
            VPT_ONLY_FLAG="--vpt-only"
            ;;
        --ppi-only)
            PPI_ONLY_FLAG="--ppi-only"
            ;;
        --pointing-only)
            POINTING_ONLY_FLAG="--pointing-only"
            ;;
        --ppi-map-day-extent)
            PPI_MAP_DAY_EXTENT_FLAG="--ppi-map-day-extent"
            ;;
        --ppi-map-colorbar-shrink)
            PPI_MAP_COLORBAR_SHRINK_NEXT=1
            ;;
        *)
            if [ "${PPI_MAP_COLORBAR_SHRINK_NEXT:-0}" -eq 1 ]; then
                PPI_MAP_COLORBAR_SHRINK="$arg"
                PPI_MAP_COLORBAR_SHRINK_NEXT=0
            else
                echo "Error: Unknown option '$arg'"
                echo "Supported options: -b, --skip-all-transition, --vpt-only, --ppi-only, --pointing-only, --ppi-map-day-extent, --ppi-map-colorbar-shrink VALUE"
                exit 1
            fi
            ;;
    esac
done

if [ "${PPI_MAP_COLORBAR_SHRINK_NEXT:-0}" -eq 1 ]; then
    echo "Error: --ppi-map-colorbar-shrink requires a numeric value"
    exit 1
fi

exclusive_count=$(( (${#VPT_ONLY_FLAG} > 0 ? 1 : 0) + (${#PPI_ONLY_FLAG} > 0 ? 1 : 0) + (${#POINTING_ONLY_FLAG} > 0 ? 1 : 0) ))
if [ "$exclusive_count" -gt 1 ]; then
    echo "Error: --vpt-only, --ppi-only and --pointing-only are mutually exclusive"
    exit 1
fi

# Name of the actual SLURM submit script
SUBMIT_SCRIPT="make_kepler_ccrest_m_quicklooks.sh"

echo "============================================================"
echo "CCREST-M Quicklooks Submission Helper"
echo "============================================================"

# Check if the submit script exists
if [ ! -f "$SUBMIT_SCRIPT" ]; then
    echo "Error: Submit script '$SUBMIT_SCRIPT' not found in current directory"
    echo "Current directory: $(pwd)"
    echo "Available .sh files:"
    ls -la *.sh 2>/dev/null || echo "No .sh files found"
    exit 1
fi

# Validate dates
if ! date -d "$START_DATE" >/dev/null 2>&1; then
    echo "Error: Invalid start date format: $START_DATE"
    echo "Use YYYYMMDD format (e.g., 20240201)"
    exit 1
fi

if ! date -d "$END_DATE" >/dev/null 2>&1; then
    echo "Error: Invalid end date format: $END_DATE"
    echo "Use YYYYMMDD format (e.g., 20240416)"
    exit 1
fi

# Calculate number of days
start_epoch=$(date -d "$START_DATE" +%s)
end_epoch=$(date -d "$END_DATE" +%s)
num_days=$(( ( end_epoch - start_epoch ) / 86400 ))

if [ $num_days -lt 0 ]; then
    echo "Error: End date must be after start date"
    echo "Start: $START_DATE"
    echo "End: $END_DATE"
    exit 1
fi

echo "Date range: $START_DATE to $END_DATE"
echo "Number of days: $((num_days + 1))"
echo "Array range: 0-$num_days"
if [ -n "$BOUNDARY_LAYER_FLAG" ]; then
    echo "Boundary layer mode: enabled"
fi
if [ -n "$SKIP_ALL_TRANSITION_FLAG" ]; then
    echo "Skip all-transition sweeps: enabled"
fi
if [ -n "$VPT_ONLY_FLAG" ]; then
    echo "Plot mode: VPT only"
fi
if [ -n "$PPI_ONLY_FLAG" ]; then
    echo "Plot mode: PPI + PPI-map only"
fi
if [ -n "$POINTING_ONLY_FLAG" ]; then
    echo "Plot mode: Pointing only"
fi
if [ -n "$PPI_MAP_DAY_EXTENT_FLAG" ]; then
    echo "PPI-map extent mode: day outermost bounds"
fi
if [ -n "$PPI_MAP_COLORBAR_SHRINK" ]; then
    echo "PPI-map colorbar shrink: $PPI_MAP_COLORBAR_SHRINK"
fi
echo ""

# Create backup of submit script
cp "$SUBMIT_SCRIPT" "${SUBMIT_SCRIPT}.backup.$(date +%Y%m%d_%H%M%S)"

# Update the submit script with new dates and flags
echo "Updating $SUBMIT_SCRIPT with new configuration..."
sed -i "s/START_DATE=\".*\"/START_DATE=\"$START_DATE\"/" "$SUBMIT_SCRIPT"
sed -i "s/END_DATE=\".*\"/END_DATE=\"$END_DATE\"/" "$SUBMIT_SCRIPT"
# Flag variables are propagated to jobs via sbatch --export (see below)

# Verify the updates worked
if grep -q "START_DATE=\"$START_DATE\"" "$SUBMIT_SCRIPT" && grep -q "END_DATE=\"$END_DATE\"" "$SUBMIT_SCRIPT"; then
    echo "Successfully updated configuration in $SUBMIT_SCRIPT"
else
    echo "Failed to update configuration in $SUBMIT_SCRIPT"
    echo "Please check the file format and try again"
    exit 1
fi

# Create log directory
mkdir -p slurm_logs

# Submit the job
echo ""
echo "Submitting job..."
echo "Command: sbatch --export=ALL,START_DATE=$START_DATE,END_DATE=$END_DATE,BOUNDARY_LAYER_FLAG=$BOUNDARY_LAYER_FLAG,SKIP_ALL_TRANSITION_FLAG=$SKIP_ALL_TRANSITION_FLAG,VPT_ONLY_FLAG=$VPT_ONLY_FLAG,PPI_ONLY_FLAG=$PPI_ONLY_FLAG,POINTING_ONLY_FLAG=$POINTING_ONLY_FLAG,PPI_MAP_DAY_EXTENT_FLAG=$PPI_MAP_DAY_EXTENT_FLAG,PPI_MAP_COLORBAR_SHRINK=$PPI_MAP_COLORBAR_SHRINK --array=0-$num_days $SUBMIT_SCRIPT"
echo ""

sbatch_output=$(sbatch \
    --export=ALL,START_DATE="$START_DATE",END_DATE="$END_DATE",BOUNDARY_LAYER_FLAG="$BOUNDARY_LAYER_FLAG",SKIP_ALL_TRANSITION_FLAG="$SKIP_ALL_TRANSITION_FLAG",VPT_ONLY_FLAG="$VPT_ONLY_FLAG",PPI_ONLY_FLAG="$PPI_ONLY_FLAG",POINTING_ONLY_FLAG="$POINTING_ONLY_FLAG",PPI_MAP_DAY_EXTENT_FLAG="$PPI_MAP_DAY_EXTENT_FLAG",PPI_MAP_COLORBAR_SHRINK="$PPI_MAP_COLORBAR_SHRINK" \
    --array=0-$num_days "$SUBMIT_SCRIPT" 2>&1)
sbatch_exit=$?

if [ $sbatch_exit -eq 0 ]; then
    echo "Job submitted successfully"
    echo "$sbatch_output"
    job_id=$(echo "$sbatch_output" | grep -o '[0-9]*')
    echo ""
    echo "Job monitoring commands:"
    echo "  squeue -u \$USER"
    echo "  squeue -j $job_id"
    echo "  scancel $job_id    # to cancel all array tasks"
    echo ""
    echo "Log files will be in: slurm_logs/"
else
    echo "Job submission failed"
    echo "Error output: $sbatch_output"
    exit 1
fi
