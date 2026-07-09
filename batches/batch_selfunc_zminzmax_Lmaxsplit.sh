#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=selfunc_zminzmax_Lmaxsplit
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env
zmin=(0.0 1.0 2.0 3.0)
zmax=(1.0 2.0 3.0 4.6)
# n_zbins=(4 4 4 3 3 2)
n_Lbins=(4 3 2)
G=20.5

for j in "${!n_Lbins[@]}" 
do 
    # for threshold in 'True' 'False'
    # do
    cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
    python make_catalogs.py ${G} --z_bins 0.0 1.0 2.0 3.0 4.6 --n_Lbins ${n_Lbins[j]} --no-threshold 
    python make_catalogs.py ${G} --z_bins 0.0 1.0 2.0 3.0 4.6 --n_Lbins ${n_Lbins[j]} --threshold 
    
        for i in "${!zmin[@]}" 
        do 
            for ((lb=0; lb<${n_Lbins[j]}; lb++))
            do
                for L in "_Lmaxsplit${n_Lbins[j]}" "_Lsplit${n_Lbins[j]}"
                do
                # zb is the z-bin index, lb is the L-bin index
                # if [ "${threshold}" == 'True' ]; then
                    # fname="_zmin${zmin[i]}zmax${zmax[i]}_Lmaxsplit${n_Lbins[j]}bin${lb}"
                # else
                    fname="_zmin${zmin[i]}zmax${zmax[i]}${L}bin${lb}"
                # fi
                    cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/batches
                    cat > G${G}${fname}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G${fname}
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
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
done