#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=mocks_loop_wp
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
conda activate quaia-env

G=20
# for G in '20.5' '20'
# do
for zbin in '0' '1'
do
    for ((i=0;i<=99;i++))
    do
        for b in '30' '40'
        do
        
            # cat > wtheta_G${G}_zsplit2bin${zbin}_mock${i}.sh << EOF
            cat > wp_G${G}_zsplit2bin${zbin}_mock${i}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=wp_G${G}_zsplit2bin${zbin}_mock_debug${i}
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

python mocks_loop_debug.py ${G} ${zbin} ${i} ${b}
EOF

    # sbatch wtheta_G${G}_zsplit2bin${zbin}_mock${i}.sh
            sbatch wp_G${G}_zsplit2bin${zbin}_mock${i}.sh
        done
    done
done
# done
#done