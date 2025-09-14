#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=data_all_25x
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

for G in '20.5' '20.0'
do
    
    # compute clustering
    python data_all_25x.py ${G}

    done
done