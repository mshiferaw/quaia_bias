#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G20.5_b20_zsplit2_Lmaxsplit2_wtheta_jackknife
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

# python -u jackknife.py wtheta 20.5 0.001 10 30 log --L_method maxsplit 50 --b 20 Planck --n_zbins  --n_Lbins 
# python -u jackknife.py wtheta 20.5 0.001 10 30 log --L_method split 50 --b 20 Planck --n_zbins  --n_Lbins 
python -u jackknife.py wtheta 20.5 0.001 10 30 log --L_method maxsplit 50 --b 20 Planck --n_zbins 2 --n_Lbins 2
