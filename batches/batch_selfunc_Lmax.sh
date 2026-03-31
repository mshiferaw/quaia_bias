#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=selfunc_zminzmax_Lmax
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env
zmin=(0.0 1.0 2.0 3.0)
zmax=(1.0 2.0 3.0 4.6)
G=20.5

for Lmax in '-20.0' '-25.0' '-26.0' '-27.0'
do 
    cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
    python make_catalogs.py ${G} --z_bins 0.0 1.0 2.0 3.0 4.6 --L_max ${Lmax}
    for j in "${!zmin[@]}" 
    do 
        fname="_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax}"
        cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/batches
        cat > G${G}${fname}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=G${G}_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax}
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

python selection_function_map.py ../data/quaia_G${G}${fname}.fits ../data/maps/selection_function_NSIDE64_G${G}${fname}.fits -p ../data/quaia_G${G}.fits
python generate_random.py ${G} 25 ${fname} 
EOF

            sbatch G${G}${fname}.sh
        done
    done
done