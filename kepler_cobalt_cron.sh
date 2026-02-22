#!/bin/bash 
#SBATCH --partition=short-serial
#SBATCH --job-name=cobalt-proc
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=128G

source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate cao_3_11

time /home/users/cjwalden/git/kepler-radar-utils-cobalt/proc_kepler2ncas_cobalt.py -l




