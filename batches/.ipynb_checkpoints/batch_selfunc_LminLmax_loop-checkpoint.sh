#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=selfunc_zminzmax_LminLmax
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

Lmax=(-20.0 -24.0 -25.0 -26.0 -27.0)
Lmin=(-24.0 -25.0 -26.0 -27.0 -32.0)
# for Lmax in '-20.0' '-24.0' '-25.0' '-26.0' '-27.0'
for i in "${!Lmin[@]}"
do 
    # cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
    # python make_catalogs.py ${G} --z_bins 0.0 1.0 2.0 3.0 4.6 --L_max ${Lmax}
    # Lmin="${Lmin_vals[$i]}"
    # Lmax="${Lmax_vals[$i]}"
    for j in "${!zmin[@]}" 
    do 
        fname="_zmin${zmin[j]}zmax${zmax[j]}_Lmin${Lmin[i]}Lmax${Lmax[i]}"
        cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/batches
        cat > G${G}${fname}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G${G}${fname}
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

python selection_function_map.py ../data/quaia_G${G}${fname}.fits ../data/maps/selection_function_NSIDE64_G${G}${fname}.fits -p ../data/quaia_G${G}.fits -m dust stars m10 mcs unwise unwisescan mcsunwise zodi1.25 zodi3.4 zodi4.6 --overwrite_ypred
python generate_random.py ${G} 25 ${fname} 
EOF

            sbatch G${G}${fname}.sh
        done
    done
done