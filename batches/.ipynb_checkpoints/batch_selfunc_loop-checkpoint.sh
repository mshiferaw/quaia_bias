#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=selfunc_zminzmax_Lbolmin
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts

for Lbolmin in '44.5' '45.0' '45.5' '46.0' '46.5'
do 
    python make_catalogs.py 20.5 --z_bins 0.0 0.1875 0.375 0.5625 0.75 0.9375 1.125 1.3125 1.5  --L_bolmin ${Lbolmin}
    for zmin in '0.0000' '0.1875' '0.3750' '0.5625' '0.7500' '0.9375' '1.1250' '1.3125'
    do 
        zmax=$(printf '%.4f' $(bc <<< "$zmin + 0.1875"))
        cat > G${G}_zmin${zmin}zmax${zmax}_Lbolmin${Lbolmin}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=G${G}_zmin${zmin}zmax${zmax}_Lbolmin${Lbolmin}
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

python selection_function_map.py ../data/quaia_G20.5_zmin${zmin}zmax${zmax}_Lbolmin${Lbolmin}.fits ../data/maps/selection_function_NSIDE64_G20.5_zmin${zmin}zmax${zmax}_Lbolmin${Lbolmin}.fits -p ../data/quaia_G20.5.fits
python generate_random.py 20.5 25 _zmin${zmin}zmax${zmax}_Lbolmin${Lbolmin} 
EOF

            G${G}_zmin${zmin}zmax${zmax}_Lbolmin${Lbolmin}.sh
        done
    done
done