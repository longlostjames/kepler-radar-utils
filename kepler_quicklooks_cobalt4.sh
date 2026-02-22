#!/bin/bash 
#SBATCH --partition=standard
#SBATCH --job-name=cobalt-1
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --array=10-28
#SBATCH --mem=128G
#SBATCH --account=ncas_radar
#SBATCH --qos=standard

source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate cao_3_11

time /home/users/cjwalden/git/kepler-radar-utils-cobalt/make_cobalt_quicklooks.py -d 202502${SLURM_ARRAY_TASK_ID}


