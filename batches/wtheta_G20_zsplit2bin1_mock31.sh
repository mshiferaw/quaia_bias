#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=wtheta_G20_zsplit2bin1_mock31
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

python mocks_loop.py 20 1 31
