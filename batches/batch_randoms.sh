#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=random_25
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate quaia-env
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts

# for G in '20.0' '20.5'
# do
#     for zmin in '0.0' '1.0' '2.0' '3.0'
#     do 
#         zmax=$(bc <<< "$zmin + 1.0")
#         # python generate_random.py ${G} 100 _zmin${zmin}zmax${zmax}
#         python generate_random.py ${G} 25 _zmin${zmin}zmax${zmax}
#     done
# done

# for G in '20.0' '20.5'
# do
#     for facrand in '25' '50' '100'
#     do 
#         for allsky in '' '_allsky'
#         do
#             python generate_random.py ${G} ${facrand} '' ${allsky}
#         done
#     done
# done

# for G in '20.0' '20.5'
# do
#     for Lmin in '-32.0' '-31.0' '-30.0' '-29.0' '-28.0' '-27.0' '-26.0' '-25.0' '-24.0' '-23.0' '-22.0' '-21.0' '-20.0'
#     do
#         Lmax=$(bc <<< "$Lmin + 1.0")
#         python generate_random.py ${G} 25 _Lmin${Lmin}Lmax${Lmax}
#     done
# done

for G in '20.0' '20.5'
do
    python generate_random.py ${G} 25 _Lmin-32.0Lmax-31.0
done