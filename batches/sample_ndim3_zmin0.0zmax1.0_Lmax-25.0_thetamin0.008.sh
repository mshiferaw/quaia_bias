#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=pool_ndim3_zmin0.0zmax1.0_Lmax-25.0_thetamin0.008
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo -e "Starting batch job\n"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

python -u sample.py 3 0.0 1.0 -25.0 --plot --thetamin 0.008
