#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=pool_ndim3_zmin2.0zmax3.0_Lmax-27.0_thetamin0.007
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo -e "Starting batch job\n"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

python -u sample.py 3 2.0 3.0 -27.0 --plot --thetamin 0.007
