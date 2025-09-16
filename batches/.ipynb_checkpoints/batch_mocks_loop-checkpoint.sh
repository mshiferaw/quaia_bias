#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=mocks_loop
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env

#for G in '20' '20.5'
#do
for zbin in '0' '1'
do
    for ((i=0;i<=99;i++))
    do
        # cat > wtheta_G${G}_zsplit2bin${zbin}_mock${i}.sh << EOF
        cat > wtheta_G20_zsplit2bin${zbin}_mock${i}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=wtheta_G20_zsplit2bin${zbin}_mock${i}
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

python mocks_loop.py 20 ${zbin} ${i}
EOF

        # sbatch wtheta_G${G}_zsplit2bin${zbin}_mock${i}.sh
        sbatch wtheta_G20_zsplit2bin${zbin}_mock${i}.sh
    done
done
#done