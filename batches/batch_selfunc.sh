#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=selfunc_G20.5_Lsplit4bin0
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts

# for zmin in '0.0' '1.0' '2.0' '3.0'
# do 
#     zmax=$(bc <<< "$zmin + 1.0")
#     python selection_function_map.py ../data/quaia_G20.0_zmin${zmin}zmax${zmax}.fits ../data/maps/selection_function_NSIDE64_G20.0_zmin${zmin}zmax${zmax}.fits -p ../data/quaia_G20.0.fits
# done

# for Lbin in '0' '1' '2'
# do 
#     for G in '20.0' '20.5'
#     do
#         python selection_function_map.py ../data/quaia_G${G}_Lsplit3bin${Lbin}.fits ../data/maps/selection_function_NSIDE64_G${G}_Lsplit3bin${Lbin}.fits -p ../data/quaia_G${G}.fits
#     done
# done

# for zbin in '0' '1'
# do 
#     python selection_function_map.py ../data/quaia_G20.5_zsplit2bin${zbin}.fits ../data/maps/selection_function_NSIDE64_G20.5_zsplit2bin${zbin}.fits -p ../data/quaia_G20.5.fits
# done

python selection_function_map.py ../data/quaia_G20.5_Lsplit4bin0.fits ../data/maps/selection_function_NSIDE64_G20.5_Lsplit4bin0.fits -p ../data/quaia_G20.5.fits