#!/bin/bash 
#SBATCH --partition=high-mem 
#SBATCH --job-name=coalesc3_quicklooks
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

# Activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up paths
SCRIPT_DIR="/home/users/cjwalden/git/kepler-radar-utils-cobalt"
PYTHON_SCRIPT="$SCRIPT_DIR/make_coalesc3_quicklooks.py"
INPATH=${INPATH:-/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/coalesc3/L1_v1.0.1}
OUTPATH=${OUTPATH:-$INPATH/quicklooks}

# Get date from command line or use environment variable
DATESTR=${1:-${DATESTR}}

if [ -z "$DATESTR" ]; then
    echo "Error: No date specified. Usage: $0 YYYYMMDD"
    echo "Or set DATESTR environment variable"
    exit 1
fi

echo "COALESC3 Quicklooks Generation"
echo "Date: ${DATESTR}"
echo "Input path: ${INPATH}"
echo "Output path: ${OUTPATH}"
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "------------------------------"

# Run the quicklook generation script
time python $PYTHON_SCRIPT -d ${DATESTR} -i ${INPATH} -o ${OUTPATH}

echo "Completed quicklook generation for ${DATESTR}"
