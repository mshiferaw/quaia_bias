#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G20.5_b30_zsplit3_Lsplit3_wtheta
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

# python -u clustering.py wtheta 20.5 0.001 10 30 log --L_method maxsplit _MC_ --b 30 Planck --n_zbins  --n_Lbins 
# python -u clustering.py wtheta 20.5 0.001 10 30 log --L_method split _MC_ --b 30 Planck --n_zbins  --n_Lbins 
python -u clustering.py wtheta 20.5 0.001 10 30 log --L_method split _MC_ --b 30 Planck --n_zbins 3 --n_Lbins 3
