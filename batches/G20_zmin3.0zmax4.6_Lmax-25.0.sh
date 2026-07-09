#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=G20_zmin3.0zmax4.6_Lmax-25.0
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

python selection_function_map.py ../data/quaia_G20_zmin3.0zmax4.6_Lmax-25.0.fits ../data/maps/selection_function_NSIDE64_G20_zmin3.0zmax4.6_Lmax-25.0.fits -p ../data/quaia_G20.fits
python generate_random.py 20 25 _zmin3.0zmax4.6_Lmax-25.0 
