#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=selfunc_LminLmax
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=96:00:00

#conda init 
conda activate quaia-env
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts

for G in '20.0' '20.5'
do 

    python selection_function_map.py ../data/quaia_G${G}_Lmin-89.0Lmax-28.0.fits ../data/maps/selection_function_NSIDE64_G${G}_Lmin-89.0Lmax-28.0.fits -p ../data/quaia_G${G}.fits
    
    for Lmin in '-32.0' '-28.0' '-24.0'
    do
        Lmax=$(bc <<< "$Lmin + 4.0")
        python selection_function_map.py ../data/quaia_G${G}_Lmin${Lmin}Lmax${Lmax}.fits ../data/maps/selection_function_NSIDE64_G${G}_Lmin${Lmin}Lmax${Lmax}.fits -p ../data/quaia_G${G}.fits
    done
done
