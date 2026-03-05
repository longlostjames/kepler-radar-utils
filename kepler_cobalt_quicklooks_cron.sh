#!/bin/bash 
#SBATCH --partition=short-serial
#SBATCH --job-name=cobalt-quicklooks
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --mem=128G

source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

# Set up script path (use SLURM_SUBMIT_DIR for SLURM jobs)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

time $SCRIPT_DIR/make_cobalt_quicklooks_latest.py 



