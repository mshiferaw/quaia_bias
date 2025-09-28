#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=LminLmax
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env

for G in '20.5' '20.0'
do
    for Lmin in '-32.0' '-25.0' 
    do
        Lmax=$(bc <<< "$Lmin + 5.0")
        fn_gcat_zbin="../data/quaia_G${G}_Lmin${Lmin}Lmax${Lmax}.fits"
        fn_selfunc="../data/maps/selection_function_NSIDE64_G${G}_Lmin${Lmin}Lmax${Lmax}.fits"
        
        cat > sel_func_G${G}_Lmin${Lmin}Lmax${Lmax}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=sel_func_G${G}_Lmin${Lmin}Lmax${Lmax}
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

# make catalogs
python make_catalogs.py ${G} --L_bins ${Lmin} ${Lmax}

# generate selfunc
python selection_function_map.py $fn_gcat_zbin $fn_selfunc -p ../data/quaia_G${G}.fits

# generate randoms
python generate_random.py ${G} 25 _Lmin${Lmin}Lmax${Lmax}
EOF
    
        sbatch sel_func_G${G}_Lmin${Lmin}Lmax${Lmax}.sh
    done
done