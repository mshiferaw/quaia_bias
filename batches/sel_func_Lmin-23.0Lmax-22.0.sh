#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=sel_func_G20.0_Lmin-23.0Lmax-22.0
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env
python selection_function_map.py ../data/quaia_G20.0_Lmin-23.0Lmax-22.0.fits ../data/maps/selection_function_NSIDE64_G20.0_Lmin-23.0Lmax-22.0.fits -p ../data/quaia_G20.0.fits
