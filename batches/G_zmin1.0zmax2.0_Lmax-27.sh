#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=G_zmin1.0zmax2.0_Lmax-27
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

python selection_function_map.py ../data/quaia_G20.5__zmin1.0zmax2.0_Lmax-27.fits ../data/maps/selection_function_NSIDE64_G20.5__zmin1.0zmax2.0_Lmax-27.fits -p ../data/quaia_G20.5.fits
python generate_random.py 20.5 25 _zmin1.0zmax2.0_Lmax-27 
