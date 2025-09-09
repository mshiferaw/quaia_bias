#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=selfunc_Lsplit
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts

# for Lbin in '0' '1' '2'
# do 
#     for G in '20.0' '20.5'
#     do
#         python selection_function_map.py ../data/quaia_G${G}_Lsplit3bin${Lbin}.fits ../data/maps/selection_function_NSIDE64_G${G}_Lsplit3bin${Lbin}.fits -p ../data/quaia_G${G}.fits
#     done
# done

for G in '20.0' '20.5'
do
    python selection_function_map.py ../data/quaia_G${G}_Lsplit3bin2.fits ../data/maps/selection_function_NSIDE64_G${G}_Lsplit3bin2.fits -p ../data/quaia_G${G}.fits
done