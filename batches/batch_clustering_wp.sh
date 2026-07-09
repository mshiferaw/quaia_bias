#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=clustering_wp_minmax_log
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

# python -u clustering.py wp 20.5 1 200 14 log #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wtheta 20.5 0.1 25 10 log #20 30 180 30 linear #1 200 14 log
# python -u clustering.py wp 20.5 30 180 30 linear max 0.5 --z_bins 0.0,1.0,2.0,3.0 #1 200 14 log
# python -u clustering.py wp 20.5 30 180 30 linear max 0.8 --z_bins 3.0,4.6 #1 200 14 log
# python -u clustering.py wp 20.5 30 180 30 linear --L max _ --b 30 --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 1 200 14 log --L max _ --b 0 --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 1 200 14 log --L max _ --b 10 --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 1 200 14 log --L max _ --b 20 --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 1 200 14 log --L max _ --b 30 --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 1 200 14 log --L max _ --b 40 --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 1 200 14 log --L max _ --b 50 --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 1 200 14 log --L max _ --b 60 --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 0.15 240 20 log --L max --mask 50 _MC_ Quaia --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 0.15 240 20 log --L max _MC_ --b 30 Quaia --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.0 0.15 240 20 log --L max --mask 50 _MC_ Quaia --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.0 0.15 240 20 log --L max _MC_ --b 30 Quaia --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.0 0.15 240 20 log --L max --mask 50 _MC_ Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.0 0.15 240 20 log --L max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.0 30 180 30 linear --L max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.5 30 180 30 linear --L minmax --mask 50 _MC_ Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# python -u clustering.py wp 20.0 30 180 30 linear --L minmax _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
python -u clustering.py wp 20.5 1 200 14 log --L minmax _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 