#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=random_LminLmax
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts

for G in '20.0' '20.5'
do
    for Lmin in '-32' '-29' '-26' '-23'
    do
        Lmax=$(bc <<< "$Lmin + 3.0")
        python generate_random.py ${G} 25 _Lmin${Lmin}Lmax${Lmax}
    done
done
