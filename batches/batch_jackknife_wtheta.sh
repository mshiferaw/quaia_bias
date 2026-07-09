#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=jackknife_wtheta_max
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

# python -u jackknife.py wtheta 20.5 0.1 25 10 log max 14 0.5 --z_bins 0.0,1.0,2.0,3.0,4.6 #20 30 180 30 linear #1 200 14 log
# python -u jackknife.py wp 20.5 30 180 30 linear minmax 50 #1 200 14 log
# python -u jackknife.py wtheta 20.5 0.1 25 10 log --L max 14 --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u jackknife.py wtheta 20.5 0.1 25 10 log --L max 14 --mask 50 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u jackknife.py wtheta 20.5 0.001 10 20 log --L max 50 --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u jackknife.py wtheta 20.5 0.001 10 20 log --L minmax 50 --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# python -u jackknife.py wtheta 20.5 0.001 10 30 log --L_method minmax 50 --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# python -u jackknife.py wtheta 20.5 0.001 10 30 log --L_method max 50 --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
python -u jackknife.py wtheta 20.5 0.001 10 24 log --L_method max 50 --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 