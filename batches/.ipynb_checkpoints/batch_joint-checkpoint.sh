#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=joint
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
conda activate quaia-env
Lmin=(-32.0 -27.0 -26.0 -25.0)
Lmax=(-27.0 -26.0 -25.0 -20.0)
zmin=(0.0 1.0 2.0 3.0)
zmax=(1.0 2.0 3.0 4.6)

for G in 20.5 20.0
do

    cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts

    # make catalogs
    python make_catalogs.py ${G} --L_bins -32 -27 -26 -25 -20 --z_bins 0 1 2 3 4.6
    
    cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/batches

    for i in "${!Lmin[@]}" 
    do        
        for j in "${!zmin[@]}" 
        do
            fname="_zmin${zmin[j]}zmax${zmax[j]}_Lmin${Lmin[i]}Lmax${Lmax[i]}"
            fn_gcat_zbin="../data/quaia_G${G}${fname}.fits"
            fn_selfunc="../data/maps/selection_function_NSIDE64_G${G}${fname}.fits"
            cat > sel_func_G${G}${fname}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=sel_func_G${G}${fname}
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

# generate selfunc
python selection_function_map.py $fn_gcat_zbin $fn_selfunc -p ../data/quaia_G${G}.fits

# generate randoms
python generate_random.py ${G} 25 ${fname}
EOF

            sbatch sel_func_G${G}${fname}.sh
        done
    done
done