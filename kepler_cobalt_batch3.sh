#!/bin/bash 
#SBATCH --partition=short-serial
#SBATCH --job-name=cobalt-3
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --array=10-18
#SBATCH --mem=128G

source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate cao_3_11

time /home/users/cjwalden/git/kepler-radar-utils-cobalt/proc_kepler2ncas_cobalt.py -d 202502${SLURM_ARRAY_TASK_ID}



