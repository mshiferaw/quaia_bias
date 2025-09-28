#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=data_25x_debug
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
conda activate quaia-env

G=20
# for G in '20.5' '20.0'
# do
for zbin in '0' '1'
do

    fn_gcat_zbin="../data/quaia_G${G}_zsplit2bin${zbin}.fits"
    fn_selfunc="../data/maps/selection_function_NSIDE64_G${G}_zsplit2bin${zbin}.fits"

    cat > sel_func_zsplit2bin${zbin}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=sel_func_G${G}_zsplit2bin${zbin}
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
conda activate quaia-env

# # generate selfunc
# python selection_function_map.py $fn_gcat_zbin $fn_selfunc -p ../data/quaia_G${G}.fits

# # generate randoms
# python generate_random.py ${G} 10 _zsplit2bin${zbin}

# compute clustering
python data_debug.py ${G} ${zbin}
EOF

    sbatch sel_func_zsplit2bin${zbin}.sh
done
# done