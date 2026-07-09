#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G20.5_b60
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

python -u clustering.py wp 20.5 1 200 14 log --L max _ --b 60 --z_bins 0.0 1.0 2.0 3.0 4.6 
