#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=sel_func_zsplit2bin0.fits
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

# generate selfunc
python selection_function_map.py ../data/quaia_G20.0_zsplit2bin0.fits ../data/maps/selection_function_NSIDE64_G20.0_zsplit2bin0.fits -p ../data/quaia_G20.0.fits

# generate randoms
python generate_random.py 20.0 25 _zsplit2bin0

# compute clustering
python data_25x.py 20.0 0
