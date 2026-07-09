#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G_zmin2.0zmax3.0_Lmaxsplit2bin0
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

python selection_function_map.py ../data/quaia_G20.5_zmin2.0zmax3.0_Lmaxsplit2bin0.fits ../data/maps/selection_function_NSIDE64_G20.5_zmin2.0zmax3.0_Lmaxsplit2bin0.fits -p ../data/quaia_G20.5.fits -m dust stars m10 mcs unwise unwisescan mcsunwise zodi1.25 zodi3.4 zodi4.6 --overwrite_ypred
python generate_random.py 20.5 25 _zmin2.0zmax3.0_Lmaxsplit2bin0 
