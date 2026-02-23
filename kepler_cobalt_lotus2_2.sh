#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=cobalt-2
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --array=1-9
#SBATCH --mem=128G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

time /home/users/cjwalden/git/kepler-radar-utils-cobalt/proc_kepler2ncas_cobalt.py -d 2025060${SLURM_ARRAY_TASK_ID}


