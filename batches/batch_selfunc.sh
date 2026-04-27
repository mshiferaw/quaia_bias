#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=selfunc_zminzmax
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

# python selection_function_map.py ../data/quaia_G20.5_Lsplit4bin0.fits ../data/maps/selection_function_NSIDE64_G20.5_Lsplit4bin0.fits -p ../data/quaia_G20.5.fits

# python selection_function_map.py ../data/quaia_G20.0_Lmin-29.0Lmax-26.0.fits ../data/maps/selection_function_NSIDE64_G20.0_Lmin-29.0Lmax-26.0.fits -p ../data/quaia_G20.0.fits
# python selection_function_map.py ../data/quaia_G20.5_zmin3.0zmax4.6.fits ../data/maps/selection_function_NSIDE64_G20.5_zmin3.0zmax4.6.fits -p ../data/quaia_G20.5.fits

# python selection_function_map.py ../data/quaia_G20.5_zmin0.0zmax1.0_Lbolmin46.0.fits ../data/maps/selection_function_NSIDE64_G20.5_zmin0.0zmax1.0_Lbolmin46.0.fits -p ../data/quaia_G20.5.fits

# python selection_function_map.py ../data/quaia_G20.5_zmin0.0zmax1.0_Lbolmin46.5.fits ../data/maps/selection_function_NSIDE64_G20.5_zmin0.0zmax1.0_Lbolmin46.5.fits -p ../data/quaia_G20.5.fits

# for Lbolmin in '44.5' '45.0' '45.5' '46.0' '46.5'
# do 
#     python make_catalogs.py 20.5 --z_bins 0.0 0.1875 0.375 0.5625 0.75 0.9375 1.125 1.3125 1.5  --L_bolmin ${Lbolmin}
#     for zmin in '0.0000' '0.1875' '0.3750' '0.5625' '0.7500' '0.9375' '1.1250' '1.3125'
#     do 
#         # zmax=$(bc <<< "$zmin + 0.1875")
#         zmax=$(printf '%.4f' $(bc <<< "$zmin + 0.1875"))
#         python selection_function_map.py ../data/quaia_G20.5_zmin${zmin}zmax${zmax}_Lbolmin${Lbolmin}.fits ../data/maps/selection_function_NSIDE64_G20.5_zmin${zmin}zmax${zmax}_Lbolmin${Lbolmin}.fits -p ../data/quaia_G20.5.fits
#         python generate_random.py 20.5 25 _zmin${zmin}zmax${zmax}_Lbolmin${Lbolmin} 
#     done
# done

python selection_function_map.py ../data/quaia_G20.0_zmin3.0zmax4.6.fits ../data/maps/selection_function_NSIDE64_G20.0_zmin3.0zmax4.6.fits -p ../data/quaia_G20.0.fits
python generate_random.py 20.0 25 _zmin3.0zmax4.6