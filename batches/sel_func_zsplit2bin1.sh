#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=sel_func_G20_zsplit2bin1
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
conda activate quaia-env

# # generate selfunc
# python selection_function_map.py ../data/quaia_G20_zsplit2bin1.fits ../data/maps/selection_function_NSIDE64_G20_zsplit2bin1.fits -p ../data/quaia_G20.fits

# # generate randoms
# python generate_random.py 20 10 _zsplit2bin1

# compute clustering
python data_debug.py 20 1
