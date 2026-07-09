#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G20.5_b60_wtheta_jackknife
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

# python -u jackknife.py wp 20.5 1 200 14 log --L max 24 --b 60 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 logp
# python -u jackknife.py wp 20.5 0.15 240 20 log --L max 25 --b 60 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u jackknife.py wp 20.5 0.15 240 20 log --L max 25 --b 60 Quaia --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u jackknife.py wp 20.5 30 180 30 linear --L max 50 --b 60 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# python -u jackknife.py wp 20.5 30 180 30 linear --L minmax 50 --b 60 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# python -u jackknife.py wtheta 20.5 0.1 25 10 log --L max 14 --b 60 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
python -u jackknife.py wtheta 20.5 0.001 10 20 log --L max 50 --b 60 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 

