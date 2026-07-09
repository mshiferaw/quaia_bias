#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=clustering_loop_linear
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/batches
conda activate quaia-env
G=20.5

for b in '0' '10' '20' '30' '40' '50' '60'
do 
    fname="_b${b}"
    cat > G${G}${fname}_linear.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G${G}${fname}_linear
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

python -u clustering.py wp 20.5 30 180 30 linear --L max _MC_ --b ${b} --z_bins 0.0 1.0 2.0 3.0 4.6 
EOF

    sbatch G${G}${fname}_linear.sh
done