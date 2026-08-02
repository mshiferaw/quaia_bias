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

zmin=(0.0 1.0 2.0 3.0)
zmax=(1.0 2.0 3.0 4.6)
Lmax=(-20.0 -25.0 -26.0 -27.0)

for ndim in '2' '3'
do
    # for i in {0..4}
    for i in "${!Lmax[@]}"
    do 
        # for j in {0..4}
        for j in "${!zmin[@]}" 
        do
            # fname="_ndim${ndim}_i${i}_j${j}"
            fname="_ndim${ndim}_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax[i]}"
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

# python -u sample.py ${ndim} ${i} ${j} --plot
python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot
EOF

            sbatch sample${fname}.sh
        done
    done
done
