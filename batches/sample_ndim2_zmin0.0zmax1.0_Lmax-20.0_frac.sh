#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=pool_ndim2_zmin0.0zmax1.0_Lmax-20.0_frac
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo -e "Starting batch job\n"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

python -u sampler.py 2 0.0 1.0 -20.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 1
python -u sampler.py 2 0.0 1.0 -20.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 1.4
python -u sampler.py 2 0.0 1.0 -20.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 1.8
python -u sampler.py 2 0.0 1.0 -20.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 2.1
python -u sampler.py 2 0.0 1.0 -20.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 2.5
