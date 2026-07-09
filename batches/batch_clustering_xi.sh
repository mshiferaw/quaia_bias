#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=clustering_xi
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=31:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

# python -u clustering.py wp 20.5 1 200 14 log #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wtheta 20.5 0.1 25 10 log #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wp 20.5 30 180 30 linear minmax 0.5 --z_bins 0.0,1.0,2.0,3.0,4.6 #1 200 14 log
python -u clustering.py xi 20.5 1 200 14 log max 0.5 --z_bins 0.0,1.0,2.0,3.0,4.6 #1 200 14 log
