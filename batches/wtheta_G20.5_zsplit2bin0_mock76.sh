#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=wtheta_G20.5_zsplit2bin0_mock76
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

python mocks_loop.py 20.5 0 76
