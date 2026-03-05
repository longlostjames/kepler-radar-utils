#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=picasso_quicklook_single
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard
#SBATCH --output=slurm_logs/picasso_quicklooks_%j.out
#SBATCH --error=slurm_logs/picasso_quicklooks_%j.err

# PICASSO Quicklook Generation Script
#
# Usage: 
#   DATE=YYYYMMDD sbatch make_picasso_quicklooks.sh
#   or
#   bash make_picasso_quicklooks.sh YYYYMMDD
#
# Example:
#   DATE=20171214 sbatch make_picasso_quicklooks.sh

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Date from environment variable or command line
if [ -n "$DATE" ]; then
    # Use environment variable (SLURM submission)
    DATESTR=$DATE
elif [ $# -ge 1 ]; then
    # Use command line argument (direct bash execution)
    DATESTR=$1
else
    echo "Error: Date required as environment variable DATE or first argument (YYYYMMDD)"
    echo "Usage: DATE=YYYYMMDD sbatch $0"
    echo "   or: bash $0 YYYYMMDD"
    exit 1
fi

# Optional boundary layer flag
BL_FLAG=""
if [ $# -ge 2 ] && [ "$2" == "-b" ]; then
    BL_FLAG="-b"
    echo "Boundary layer mode enabled"
fi

# Default paths
INPATH="/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/picasso/L1_v1.0.0"
OUTPATH="$INPATH/quicklooks"

echo "============================================================"
echo "PICASSO Quicklook Generation"
echo "============================================================"
echo "Date: $DATESTR"
echo "Input path: $INPATH"
echo "Output path: $OUTPATH"
echo "Script: $SCRIPT_DIR/make_picasso_quicklooks.py"
echo "============================================================"

# Check if input directory exists
if [ ! -d "$INPATH/$DATESTR" ]; then
    echo "Error: Input directory does not exist: $INPATH/$DATESTR"
    exit 1
fi

# Run the quicklook script
echo "Running quicklook script..."
time python "$SCRIPT_DIR/make_picasso_quicklooks.py" -d $DATESTR $BL_FLAG

exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "============================================================"
    echo "PICASSO quicklooks completed successfully for ${DATESTR}"
    echo "============================================================"
else
    echo "============================================================"
    echo "PICASSO quicklooks failed for ${DATESTR} with exit code $exit_code"
    echo "============================================================"
fi

exit $exit_code
