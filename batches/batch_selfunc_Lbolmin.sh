#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=selfunc_Lbolmin
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts

# python selection_function_map.py ../data/quaia_G20.5_zmin0.0zmax1.0_Lbolmin45.5.fits ../data/maps/selection_function_NSIDE64_G20.5_zmin0.0zmax1.0_Lbolmin45.5.fits -p ../data/quaia_G20.5.fits

python selection_function_map.py ../data/quaia_G20.5_zmin0.0zmax1.0_Lbolmin45.0.fits ../data/maps/selection_function_NSIDE64_G20.5_zmin0.0zmax1.0_Lbolmin45.0.fits -p ../data/quaia_G20.5.fits
python generate_random.py 20.5 25 _zmin0.0zmax1.0_Lbolmin45.0
# python selection_function_map.py ../data/quaia_G20.5_zmin0.0zmax1.0_Lbolmin44.5.fits ../data/maps/selection_function_NSIDE64_G20.5_zmin0.0zmax1.0_Lbolmin44.5.fits -p ../data/quaia_G20.5.fits