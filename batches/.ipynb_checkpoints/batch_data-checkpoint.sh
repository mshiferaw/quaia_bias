#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=quaia_debug
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts

python data.py