#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=clustering_wtheta_minmax
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

# python -u clustering.py wtheta 20.5 0.001 5 27 log 0.5 --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wtheta 20.5 0.001 5 27 log --mask 0 _ --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wtheta 20.5 0.001 5 27 log _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wtheta 20.5 0.001 5 27 log --mask 50 _MC_ Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wtheta 20.5 0.1 25 10 log --L max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wtheta 20.5 0.1 25 10 log --L max --mask 50 _MC_ Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wtheta 20.5 0.001 10 20 log --L max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wtheta 20.5 0.001 10 20 log --L minmax _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6
# python -u clustering.py wtheta 20.5 0.001 10 30 log --L_method max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log
python -u clustering.py wtheta 20.5 0.001 10 24 log --L_method max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #20 30 180 30 linear #1 200 14 log