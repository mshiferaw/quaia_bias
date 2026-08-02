#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=sample_loop
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/batches
conda activate pyccl-env

for ndim in '2' '3'
do
    for i in {0..4}
    do 
        for j in {0..4}
        do
            fname="_ndim${ndim}_i${i}_j${j}"
            cat > sample${fname}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=sample${fname}
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

python -u sample.py ${ndim} ${i} ${j} --plot
EOF

            sbatch sample${fname}.sh
        done
    done
done
